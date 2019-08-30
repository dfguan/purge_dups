/*
 * =====================================================================================
 *
 *       Filename:  ngscstat.c
 *
 *    Description:  next generation seqencing coverage statistic
 *
 *        Version:  1.0
 *        Created:  17/03/2019 21:09:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>


#include "asset.h"
#include "bamlite.h"
#include "sdict.h"
#include "kvec.h"
#include "ksort.h"

typedef struct {
	int mq:15, rev:1, as:16;
	uint32_t s, e, tid;
}aln_inf_t;

typedef struct {
	/*int mq:15, rev:1, as:16;*/
	uint32_t tid:31, rev:1;
	uint32_t ntid:31, nrev:1;
	uint32_t s, e;
}aln2_inf_t;


#define cord_key(a) ((a).s)
KRADIX_SORT_INIT(cord, cors, cord_key, 4);
	
void col_cords(aln_inf_t  *fal, int min_mq, uint32_t max_is, cord_t *cc)
{
	if (fal->mq > min_mq) {
		uint32_t s = fal->s;
		uint32_t e = fal->e;
		if (e - s < max_is) {
			cors tmp = (cors) {s, e};
			cord_push(&cc[fal->tid], &tmp);	
		}
	}
}

uint32_t get_target_end(uint32_t *cigar, int n_cigar, uint32_t s)
{
	int i = 0;
	for ( i = 0; i < n_cigar; ++i) {
		uint8_t c  = bam_cigar_opchr(cigar[i]);
		if (c == 'M' || c == 'D') 
			s += cigar[i] >> BAM_CIGAR_SHIFT;
	}	
	return s;
}

int chl_col_ctgs(char *bam_fn, sdict_t *ctgs)
{
	bamFile fp;
	bam_header_t *h;
	bam1_t *b;
	fp = bam_open(bam_fn, "r"); //should check if bam is sorted
	if (fp == 0) {
		fprintf(stderr, "[E::%s] fail to open %s\n", __func__, bam_fn);
		return -1;
	}
	h = bam_header_read(fp);
	b = bam_init1();
	
	int i;
	for ( i = 0; i < h->n_targets; ++i) {
		char *name = h->target_name[i];
		uint32_t len = h->target_len[i];
		uint32_t lenl, lenr;
		lenl = lenr = len >> 1;
		sd_put(ctgs, name, len, 0);
	}
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);
	return 0;
}

/*bc_ary_t *proc_bam(char *srt_bam_fn, int min_as, int min_mq, uint32_t max_ins_len, sdict_t *ctgs, int opt)*/
int proc_bam(char *bam_fn, int min_mq, uint32_t max_is, sdict_t *ctgs, int opt, cord_t *cc)
{
	bamFile fp;
	bam_header_t *h;
	bam1_t *b;
	fp = bam_open(bam_fn, "r"); //should check if bam is sorted
	if (fp == 0) {
		fprintf(stderr, "[E::%s] fail to open %s\n", __func__, bam_fn);
		return -1;
	}
	
	b = bam_init1();

	char *cur_qn = NULL;
	long bam_cnt = 0;
	int is_set = 0;
	aln_inf_t aln;
	int aln_cnt = 0;
	uint8_t rev;
	while (1) {
		//segment were mapped 
		if (bam_read1(fp, b) >= 0 ) {
			if (!cur_qn || strcmp(cur_qn, bam1_qname(b)) != 0) {
				/*fprintf(stderr, "%d\t%d\t%d\n", aln_cnt, rev, aln.mq);*/
				if (aln_cnt == 2 && (rev == 1 || rev == 2)) col_cords(&aln, min_mq, max_is, cc);
				aln_cnt = 0;	
				rev = 0;
				is_set = 0;
				if (cur_qn) free(cur_qn); 
				cur_qn = strdup(bam1_qname(b));
			}
			if (b->core.flag & 0x4) continue; //not aligned
			++aln_cnt;
			rev = (rev << 1) | !!(b->core.flag & 0x10);
			if (b->core.isize > 0 && !is_set) {
				aln.s = opt ? get_target_end(bam1_cigar(b), b->core.n_cigar, b->core.pos) :  b->core.pos + 1;
				aln.mq = b->core.qual;
				aln.tid = b->core.tid;
				aln.e = aln.s + b->core.isize - 1;	
			} 
			if (opt && b->core.isize < 0) 
				aln.e = b->core.pos + 1;
			
			
			if ((++bam_cnt % 1000000) == 0) fprintf(stderr, "[M::%s] processing %ld bams\n", __func__, bam_cnt); 
		} else {
			/*fprintf(stderr, "%d\t%d\t%d\n", aln_cnt, rev, aln.mq);*/
			if (aln_cnt == 2 && (rev == 1 || rev == 2)) col_cords(&aln, min_mq, max_is, cc);
			break;	
		}
	}
	fprintf(stderr, "[M::%s] finish processing %ld bams\n", __func__, bam_cnt); 
	
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);
	return 0;
}


ctg_pos_t *col_pos(cord_t *cc, int n, sdict_t *ctgs)
{
	ctg_pos_t *d = ctg_pos_init();
	int k;
	for ( k = 0;k < n; ++k) {
		ctg_pos_push(d, k);
		cors *e = cc[k].coords;
		size_t n_e = cc[k].n;
		radix_sort_cord(e, e + n_e);
		size_t i = 0, j;
		uint32_t st,ed;
		
		while (i < n_e) {
			st = e[i].s, ed= e[i].e;
					/*fprintf(stderr, "h%s\t%u\t%u\n",ctgs->seq[k].name, st, ed);*/
			for ( j = i + 1; j <= n_e; ++j) {
					/*fprintf(stderr, "t%s\t%u\t%u\n",ctgs->seq[k].name, st, ed);*/
				if (j == n_e || e[j].s != e[i].s) {
					/*fprintf(stderr, "d%s\t%u\t%u\n",ctgs->seq[k].name, st, ed);*/
					pos_push(&d->ctg_pos[k], st << 1);
					pos_push(&d->ctg_pos[k], ed << 1 | 1);
					i = j;
					break;
				} else {
					ed = max(ed, e[j].e); //get the largest one, problem ? 
				}
			} 
		}	
	} 
	return d;	
}

/*ctg_pos_t *col_pos(bc_ary_t *bc_l, uint32_t min_bc, uint32_t max_bc, uint32_t min_inner_bcn, uint32_t min_mol_len, int n_targets)*/
	/*cord_t *cc = col_cords(bc_l, min_bc, max_bc, min_inner_bcn, max_span, min_mol_len, ctgs->n_seq, ctgs);	*/

int ngscstat(char *bam_fn[], int n_bam, int min_mq, uint32_t max_is, int max_cov, int opt, char *out_dir)
{
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting contigs \n", __func__);
#endif
	sdict_t *ctgs = sd_init();
	chl_col_ctgs(bam_fn[0], ctgs);
	if (!ctgs->n_seq) {
		fprintf(stderr, "[M::%s] No contigs found in bam file, quit\n", __func__);
		return 0;
	}
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing bam file\n", __func__);
#endif
	int i;
		
	cord_t *cc = calloc(ctgs->n_seq, sizeof(cord_t));

	for ( i = 0; i < n_bam; ++i) {
		if (proc_bam(bam_fn[i], min_mq, max_is, ctgs, opt, cc)) {
			return -1;	
		}	
	}	

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting positions\n", __func__);
#endif
	ctg_pos_t *d = col_pos(cc, ctgs->n_seq, ctgs);	
	cord_destroy(cc, ctgs->n_seq);	
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] calculating coverage\n", __func__);
#endif
	cov_ary_t *ca = cal_cov(d, ctgs);

	ctg_pos_destroy(d);
	if (!ca) {
		fprintf(stderr, "[W::%s] low quality alignments\n", __func__);
		return 0;	
	}
		
	/*csel_sup_reg(ca, min_cov_rat, min_cov, max_cov_rat, max_cov, ctgs);*/
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] selecting supported regions\n", __func__);
#endif
	char *type = "TX";
	char *desc = "10x data";

	print_coverage_stat(ca, max_cov, ctgs, type, out_dir);
	print_base_coverage(ca, ctgs, type, out_dir);
	
	cov_ary_destroy(ca, ctgs->n_seq); //a little bit messy
	sd_destroy(ctgs);

	fprintf(stderr, "Program finished successfully\n");
	return 0;

}



int main(int argc, char *argv[])
{
	int c;
	int min_mq = 30;
	uint32_t max_is=1000;
	int max_cov = 500;
	char *r;
	char *out_dir = ".";
	int option = 0; //the way to calculate molecule length //internal parameters not allowed to adjust by users
	while (~(c=getopt(argc, argv, "b:B:c:C:q:S:a:L:l:O:h"))) {
		switch (c) {
			case 'q':
				min_mq = atoi(optarg);
				break;
			case 'L':
				max_is = strtol(optarg, &r, 10);
				break;
			case 'O':
				out_dir = optarg;
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: ngscstat [options] <BAM_FILEs> ...\n");
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -q    INT      minimum alignment quality [30]\n");
				fprintf(stderr, "         -M    INT      maximum read depth [500]\n");
				/*fprintf(stderr, "         -S    INT      minimum aislignment score [0]\n");*/
				fprintf(stderr, "         -L    INT      maximum insert size [1000]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] require at least one bam file!\n", __func__); goto help;
	}
	
	char **bam_fn = argv+optind;
	int n_bam = argc - optind;
	fprintf(stderr, "Program starts\n");	
	ngscstat(bam_fn, n_bam,  min_mq, max_is, max_cov, option, out_dir);
	fprintf(stderr, "Program ends\n");	
	return 0;	
}




