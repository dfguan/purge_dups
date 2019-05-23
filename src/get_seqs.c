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

typedef struct {
	uint32_t sn; //don't think there will be 2G contigs
	uint32_t s, e;
	uint32_t bst_sn;
}dup_t;

typedef struct {
	char *name;
	uint32_t s,e;
	char *tp;
	char *bst_name;
}dup_s;


typedef struct {size_t n, m; dup_t *a;} dup_v;

int print_dups(dup_t *dups, size_t n, sdict_t *dup_n)
{
	dup_t *dp = dups;
	size_t i;
	for ( i = 0; i < n; ++i ) 
		fprintf(stdout, "%s\t%u\t%u\t%s\n", dup_n->seq[dp[i].sn].name, dp[i].s, dp[i].e, dp[i].bst_sn != 0XFFFFFFFF ? dup_n->seq[dp[i].bst_sn].name : "*");
	return 0;
}
int print_dups2(dup_t *dups, size_t n, char *name)
{
	dup_t *dp = dups;
	size_t i;
	for ( i = 0; i < n; ++i ) 
		fprintf(stdout, "%s\t%u\t%u\n", name, dp[i].s, dp[i].e);
	return 0;
}
int update_best_sn(dup_v *dups, sdict_t *sn)
{
	size_t n = dups->n;	
	dup_t *dp = dups->a;
	size_t i;
	for ( i = 0; i < n; ++i) 
		if (dp[i].bst_sn != 0XFFFFFFFF)  {
			while (sn->seq[dp[i].bst_sn].type == 1) dp[i].bst_sn = sn->seq[dp[i].bst_sn].best_hit;//should we check if best_sn == -1? this can't happen
		}
	return 0;
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
		else if (t == 4) k->bst_name = q;
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 4) return -1;
	return 0;
}


int dup_idx(dup_v dups, uint64_t *idx)
{
	size_t i, last;
	dup_t *dp = dups.a;
	size_t n = dups.n;
	for ( i = 1, last = 0; i <= n; ++i ) {
		if (i == n || dp[i].sn != dp[last].sn) {
			/*fprintf(stderr, "%d\t%d\n", dp[last].sn, i - last);*/
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
	uint32_t rid, brid;	
	while (ks_getuntil(ks, KS_SEP_LINE, &buf, &dret) >= 0) {
		parse_dup(buf.s, buf.l, &d);
		if (!name || strcmp(name, d.name)) {
			if (name) { 
				free(name); 
				dup_t t = (dup_t){rid, 0, 0, -1};
				kv_push(dup_t, *dups, t);	
			}
			name = strdup(d.name);	
			rid = sd_put(sn, name, 0, 1);	
			if (d.bst_name[0] != '*') {
				brid = sd_put(sn, d.bst_name, 0, 1);
				if (d.tp[0] != 'O') sn->seq[rid].type = 1;
				sn->seq[rid].best_hit = brid;		
			} else 
				brid = -1;
			/*brid = d.bst_name[0] != '*' ? sd_put(sn, name, 0, 1): -1;	*/
			/*if (!strcmp(d.tp, "HAPLOTIG")) sn->seq[rid].type = 1, sn->seq[rid].best_hit = brid; */
			dup_t t = (dup_t){rid, 0, 0, brid};
			kv_push(dup_t, *dups, t);	
		}
		brid = d.bst_name[0] != '*' ? sd_put(sn, d.bst_name, 0, 1): -1;	
		dup_t k = (dup_t) {rid, d.s, d.e, brid};
		kv_push(dup_t, *dups, k);	
	}
	dup_t k = (dup_t) {rid, 0, 0, -1};
	kv_push(dup_t, *dups, k);	
	return 0;
}

int get_seqs_core(char *name, char *s, uint32_t l, dup_t  *dp, size_t n, uint32_t ml, sdict_t *sn)
{
	size_t i;
	if (n <= 2) {
		char *dash_poi = strchr(name, '_');
		if (dash_poi) *dash_poi = 0;
		fprintf(stdout, ">%s\n%s\n", name, s);	
		if (dash_poi) *dash_poi = '_';
		return 0;			
	}
	dp[n-1].s = dp[n-1].e = l + 1;
	/*print_dups2(dp, n, name);*/
	/*char *seq = malloc(sizeof(char) * (l + 1 + (n-2) * 200));*/
	char *seq = malloc(sizeof(char) * (l + 1));
	char *hapseq = malloc(sizeof(char) * (l + 1));
	uint32_t poi = 0;
	for ( i = 0; i < n - 1; ++i) {
		uint32_t st, ed;
		st = dp[i].e;
		ed = dp[i+1].s;
		memcpy(seq + poi, s + st, ed - st - 1);		
		poi += (ed - st - 1);
		/*if (i != n - 2)	{*/
			/*memset(seq + poi, 'N', 200); 	*/
			/*poi += 200;*/
		/*}*/
		seq[poi] = 0;
	}
	if (poi) {
		char *dash_poi = strchr(name, '_');	
		if (dash_poi) *dash_poi = 0;
		fprintf(stdout, ">%s\n%s\n", name, seq);
		if (dash_poi) *dash_poi = '_';
	} 
	uint32_t happoi = 0, outpoi;
	for (i = 1; i < n - 1; ++i) {

		/*fprintf(stderr, "%s\t%d\t%s\n", name, i, dp[i].bst_sn != 0XFFFFFFFF ? sn->seq[dp[i].bst_sn].name : "null");*/
		uint32_t st, ed;
		st = dp[i].s;
		ed = dp[i].e;
		outpoi = happoi;
		memcpy(hapseq + happoi, s + st - 1, ed - st + 1);	
		happoi += (ed - st + 1);
		hapseq[happoi] = 0;
		if (dp[i].bst_sn != 0XFFFFFFFF) {
			char *dash_poi = strchr(sn->seq[dp[i].bst_sn].name, '_');	
			if (dash_poi) *dash_poi = 0;
			fprintf(stderr, ">%s_%04u\n%s\n", sn->seq[dp[i].bst_sn].name, ++sn->seq[dp[i].bst_sn].aux, hapseq + outpoi); 
			if (dash_poi) *dash_poi = '_';
		}
		else {
			char *dash_poi = strchr(name, '_');	
			if (dash_poi) *dash_poi = 0;
			fprintf(stderr, ">%s_0001\n%s\n", name, hapseq + outpoi); 
			if (dash_poi) *dash_poi = '_';
		}
	}	

	/*if (poi > ) {*/
		/*fprintf(stdout, ">%s\n%s\n", name, seq);*/
		/*fprintf(stderr, ">%s\n%s\n", name, hapseq);*/
	/*} else */
		/*fprintf(stderr, ">%s\n%s\n", name, s);*/
	free(seq);
	free(hapseq);
	return 0;
}


int get_seqs(char *fafn, dup_v dups, uint64_t *idx, sdict_t *sn, uint32_t ml)
{
	gzFile fp;
	kseq_t *seq;
	/*fp = gzdopen(fileno(stdin), "r");*/
	fp = fafn && strcmp(fafn, "-")? gzopen(fafn, "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	uint32_t sid;
	dup_t *dp = dups.a;
	while (kseq_read(seq) >= 0) {
		if (~(sid = sd_get(sn, seq->name.s))) {
			/*fprintf(stderr, "%u", (uint32_t) idx[sid]);*/
			get_seqs_core(seq->name.s, seq->seq.s, seq->seq.l, &dp[(idx[sid] >> 32)], (uint32_t)idx[sid], ml, sn);
		} else 
			get_seqs_core(seq->name.s, seq->seq.s, seq->seq.l, 0, 0, ml, sn);
	} 
 	//add some basic statistics maybe 		
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
typedef struct {
	char *dup_fn;
	char *fafn;
	uint32_t ml;
} opt_t;

//this bed file is sorted by name
int main(int argc, char *argv[])
{
	opt_t opts;
	int c;
	char *r;
	opts.ml = 10000;
	char *program;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	while (~(c=getopt(argc, argv, "l:h"))) {
		switch (c) {
			case 'l': 
				opts.ml = strtol(optarg, &r, 10);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  [<options>] <DUPs.BED> <FASTA> \n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -l    INT      minimum sequence length [10000]\n");	
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}

	if (optind + 2 > argc) {
		fprintf(stderr,"[E::%s] require duplication list and fasta file!\n", __func__); goto help;
	}
	opts.dup_fn = argv[optind++];
	opts.fafn = argv[optind];
	sdict_t *sn = sd_init();
	dup_v dups = {0, 0, 0};	
	col_dups(opts.dup_fn, sn, &dups);	
	update_best_sn(&dups, sn);
	/*print_dups(dups.a, dups.n, sn);*/
	uint64_t *idx = calloc(sn->n_seq, sizeof(uint64_t));
	dup_idx(dups, idx);
	/*int i = 0;*/
	/*for (i = 0; i < sn->n_seq; ++i) fprintf(stderr, "%s\t%d\n", sn->seq[i].name, (uint32_t)idx[i]);	*/
	get_seqs(opts.fafn, dups, idx, sn, opts.ml);
	free(idx);
	return 0;		
}


