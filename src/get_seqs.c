/*
 * =====================================================================================
 *
 *       Filename:  get_seqs.c
 *
 *    Description:  given a list of regions, remove them from the original sequence file and generate a new one
 *
 *        Version:  1.0
 *        Created:  14/03/2019 08:28:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <zlib.h>

#include "sdict.h"
#include "kseq.h"
#include "kvec.h"

KSEQ_INIT(gzFile, gzread, gzseek)

enum dup_type {JUNK, HIGH, HAPLOTIG, PRIMARY, REPEAT, OVLP, UNKNOWN, UNDEF};
char *dup_type_s[] = {"JUNK", "HIGHCOV", "HAPLOTIG", "PRIMARY", "REPEAT", "OVLP", "UNKNOWN"};
typedef struct {
	uint32_t sn; //don't think there will be 2G contigs
	uint32_t s, e;
	int tp;
}dup_t;

typedef struct {
	char *name;
	char *tp;
	uint32_t s,e;
}dup_s;

typedef struct {
	char		*dup_fn;
	char		*fafn;
	char		*pref;
	uint32_t	ml;
	int			add_hap_pref;
	float		mlp;
	int			gs4dup;
	int			kh, skipm;
} opt_t;

typedef struct {size_t n, m; dup_t *a;} dup_v;

int print_dups(dup_t *dups, size_t n, sdict_t *dup_n)
{
	dup_t *dp = dups;
	size_t i;
	for ( i = 0; i < n; ++i ) 
		fprintf(stdout, "%s\t%u\t%u\n", dup_n->seq[dp[i].sn].name, dp[i].s, dp[i].e);
	return 0;
}
int print_dups2(dup_t *dups, size_t n, char *name)
{
	dup_t *dp = dups;
	size_t i;
	for ( i = 0; i < n; ++i ) 
		fprintf(stderr, "%s\t%u\t%u\t%s\n", name, dp[i].s, dp[i].e, ~dp[i].tp ? dup_type_s[dp[i].tp] : "0");
	return 0;
}
int which_tp(char *tp)
{
	int i;
	for ( i = JUNK; i < UNDEF; ++i) {
		if (!strcmp(dup_type_s[i], tp))
			return i;
	}
	return -1;
}
int parse_dup(char *s, int l, dup_s *k)
{
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) k->name= q;
		else if (t == 1) k->s = strtol(q, &r, 10);
		else if (t == 2) k->e = strtol(q, &r, 10);
		else if (t == 3) k->tp = q; 
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 2) return -1;
	return 0;
}


int dup_idx(dup_v dups, uint64_t *idx)
{
	size_t i, last;
	dup_t *dp = dups.a;
	size_t n = dups.n;

	for ( i = 1, last = 0; i <= n; ++i ) {
		if (i == n || dp[i].sn != dp[last].sn) {
			idx[dp[last].sn] = (uint64_t) (last << 32) | (i - last);
			last = i;
		}
	}	
	return 0;
}

int col_dups(char *fn, sdict_t *sn, dup_v *dups)
{
	kstream_t *ks;
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	dup_s d;
	kstring_t buf = {0, 0, 0};
	int dret;
	char *name = 0;
	uint32_t rid;	
	dup_t t;
	while (ks_getuntil(ks, KS_SEP_LINE, &buf, &dret) >= 0) {
		parse_dup(buf.s, buf.l, &d);
		if (!name || strcmp(name, d.name)) {
			if (name) { 
				free(name);
				dup_t t = (dup_t){rid, 0, 0, -1}; // add an end to the previous contig
				kv_push(dup_t, *dups, t);	
			}
			name = strdup(d.name);	
			rid = sd_put(sn, name, 0, 1);	
			t = (dup_t){rid, 0, 0, -1}; //add a start to the current contig
			kv_push(dup_t, *dups, t);	
		}
		t = (dup_t) {rid, d.s, d.e, which_tp(d.tp)};
		kv_push(dup_t, *dups, t);	
	}
	t = (dup_t) {rid, 0, 0, -1}; // add an end to the final contig
	kv_push(dup_t, *dups, t);	
	return 0;
}

int get_seqs_core(char *name, char *s, uint32_t l, dup_t  *dp, size_t n, uint32_t ml, FILE *hp, FILE* pp, int ahp, float mlp, int kh, uint32_t gs4dup, int skipm)
{
	//stop here, don't know how to get 
	size_t i, j;
	if (n <= 2) {
		fprintf(pp, ">%s\n%s\n", name, s);	
		return 0;			
	}
	dp[n-1].s = dp[n-1].e = l;
	//merge dups
	for (i = 1, j = 0; i <= n; ++i) {
		if (i == n) { ++j; break;}
		if (dp[i].s < dp[j].e + gs4dup) {
			dp[j].e = dp[i].e; 
			if (~dp[i].tp)	dp[j].tp = dp[i].tp;		
		} else 
			dp[++j] = dp[i]; 	
	}
	//skip middle of contigs
	if (skipm) {
		size_t z;
		for ( i = 0, z = 0; i < j; ++i) {
			if (dp[i].s == 0 || dp[i].e == l)	
				dp[z++] = dp[i];	
			else fprintf(stderr, "[W::%s] %s %u %u is skipped\n",__func__, name, dp[i].s, dp[i].e);
		}  
		j = z;
	}

	/*print_dups2(dp, j, name);*/
	/*char *seq = malloc(sizeof(char) * (l + 1 + (n-2) * 200));*/
	char *seq = malloc(sizeof(char) * (l + 1 + (n-2) * 200));
	char *hapseq = malloc(sizeof(char) * (l + 1 + (n-2) * 200));
	uint32_t poi = 0;
	for ( i = 0; i + 1 < j; ++i) {
		uint32_t st, ed;
		st = dp[i].e;
		ed = dp[i+1].s;
		if (ed - st) {
			if (poi)	{
				memset(seq + poi, 'N', 23); 	
				poi += 23;
			}  
			memcpy(seq + poi, s + st, ed - st);		
			poi += (ed - st);
			seq[poi] = 0;
		}
	}
	uint32_t happoi = 0;
	uint32_t tp = UNKNOWN;
	for (i = 0; i < j ; ++i) {
		uint32_t st, ed;
		st = dp[i].s;
		ed = dp[i].e;
		if (ed - st) {
			tp = dp[i].tp;
			if (happoi)  {
				memset(hapseq + happoi, 'N', 23); 	
				happoi += 23;
			} 
			memcpy(hapseq + happoi, s + st, ed - st);	
			happoi += (ed - st);
			hapseq[happoi] = 0;
		}
	}	
	if (poi > ml && (float) poi / l > mlp) {
		fprintf(pp, ">%s\n%s\n", name, seq);
		if (ahp) 
			fprintf(hp, ">hap_%s %s\n%s\n", name, dup_type_s[tp], hapseq);
		else 
			fprintf(hp, ">%s %s\n%s\n", name, dup_type_s[tp], hapseq);
	} else {
		/*fprintf(stderr, "Name %s KH: %d\tTP: %d, HIGH: %d\n", name, kh, tp, HIGH);*/
		if (kh && tp == HIGH) 
			fprintf(pp, ">%s %s\n%s\n", name, dup_type_s[tp], s);
		else {
			if (ahp) 	
				fprintf(hp, ">hap_%s %s\n%s\n", name, dup_type_s[tp], s);
			else
				fprintf(hp, ">%s %s\n%s\n", name, dup_type_s[tp], s);
		}	
	}
	free(seq);
	free(hapseq);
	return 0;
}


int get_seqs(opt_t *opts, dup_v dups, uint64_t *idx, sdict_t *sn)
{
	char *fafn = opts->fafn;
	uint32_t ml = opts->ml;
	char *outpref = opts->pref;
	int ahp = opts->add_hap_pref;
	float mlp = opts->mlp;
	int kh = opts->kh;
	uint32_t gs4dup = opts->gs4dup;
	int skipm = opts->skipm;	
	
	gzFile fp;
	kseq_t *seq;
	/*fp = gzdopen(fileno(stdin), "r");*/
	fp = fafn && strcmp(fafn, "-")? gzopen(fafn, "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	uint32_t sid;
	dup_t *dp = dups.a;
	char *hap_fn = malloc(sizeof(char) * (outpref?strlen(outpref) + 16 : 16)); //hap.fa
	char *pur_fn = malloc(sizeof(char) * (outpref? strlen(outpref) + 16: 16)); //purged.fa	
	if (outpref) {
		sprintf(hap_fn, "%s%s", outpref, ".hap.fa");	
		sprintf(pur_fn, "%s%s", outpref, ".purged.fa");
	} else {
		sprintf(hap_fn, "%s","hap.fa");	
		sprintf(pur_fn, "%s","purged.fa");
	} 
	
	FILE *hap_fp = fopen(hap_fn, "w");
	FILE *purged_fp = fopen(pur_fn, "w");
	while (kseq_read(seq) >= 0) {
		if (~(sid = sd_get(sn, seq->name.s))) {
			get_seqs_core(seq->name.s, seq->seq.s, seq->seq.l, &dp[(idx[sid] >> 32)], (uint32_t)idx[sid], ml, hap_fp, purged_fp, ahp, mlp, kh, gs4dup, skipm);
		} else 
			get_seqs_core(seq->name.s, seq->seq.s, seq->seq.l, 0, 0, ml, hap_fp, purged_fp, ahp, mlp, kh, gs4dup, skipm);
	} 
 	//add some basic statistics maybe 		
	kseq_destroy(seq);
	gzclose(fp);
	free(hap_fn); free(pur_fn);
	fclose(hap_fp);
	fclose(purged_fp);
	return 0;
}

//this bed file is sorted by name
int main(int argc, char *argv[])
{
	opt_t opts;
	int c;
	char *r;
	opts.ml = 10000;
	opts.pref = 0;
	opts.add_hap_pref = 1;
	opts.mlp = 0.05;
	opts.kh = 0, opts.skipm = 0;
	opts.gs4dup = 10000;
	char *program;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	while (~(c=getopt(argc, argv, "ap:cm:l:g:he"))) {
		switch (c) {
			case 'e': 
				opts.skipm = 1;
				break;
			case 'p': 
				opts.pref = optarg;
				break;
			case 'l': 
				opts.ml = strtol(optarg, &r, 10);
				break;
			case 'g': 
				opts.gs4dup = strtol(optarg, &r, 10);
				break;
			case 'c':
				opts.kh = 1;
				break;
			case 'm': 
				opts.mlp = strtof(optarg, &r);
				break;
			case 'a': 
				opts.add_hap_pref = 0;
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  [<options>] <DUPs.BED> <FASTA> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -e    BOOL     only remove sequences at the ends of a contig [FALSE]\n");	
				fprintf(stderr, "         -p    STR      prefix of output files [NULL]\n");	
				fprintf(stderr, "         -c    BOOL     keep high coverage contigs in the primary contig set [FALSE]\n");	
				fprintf(stderr, "         -a    BOOL     do not add prefix to haplotigs [FALSE]\n");	
				fprintf(stderr, "         -g    INT      maximum gap size between duplications [10k]\n");	
				fprintf(stderr, "         -l    INT      minimum primary contig length [10k]\n");	
				fprintf(stderr, "         -m    INT      minimum ratio of remaining primary contig length to the original contig length [0.05]\n");	
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}

	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] require a duplication list and fasta file!\n", __func__); goto help;
	}
	opts.dup_fn = argv[optind++];
	opts.fafn = argv[optind];
	sdict_t *sn = sd_init();
	dup_v dups = {0, 0, 0};	
	col_dups(opts.dup_fn, sn, &dups);	
	/*print_dups(&dups, sn);*/
	uint64_t *idx = malloc(sizeof(uint64_t) * sn->n_seq);
	dup_idx(dups, idx);
	/*print_dups(dups.a, dups.n, sn);*/
	/*get_seqs(opts.fafn, dups, idx, sn, opts.ml, opts.pref, opts.add_hap_pref, opts.mlp, opts.kh, opts.gs4dup);*/
	get_seqs(&opts, dups, idx, sn);
	kv_destroy(dups);
	sd_destroy(sn);
	free(idx);
	return 0;		
}


