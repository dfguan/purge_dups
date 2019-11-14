/*
 * =====================================================================================
 *
 *       Filename:  aa_pb.c
 *
 *    Description:  assess assembly with pacbio data
 *
 *        Version:  1.0
 *        Created:  15/09/2018 10:18:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "paf.h"
#include "asset.h"
/*ctg_pos_t *col_pos(char *paf_fn, int min_as, int min_mq, sdict_t *ctgs, uint32_t flank)*/
int col_pos(char *paf_fn, int min_mq, float min_mlr, int onlyp, sdict_t *ctgs, uint32_t flank, ctg_pos_t *d)
{
	paf_file_t *fp = paf_open(paf_fn);
	if (!fp) {
		fprintf(stderr, "[E::%s] %s doesn't exist!\n", __func__, paf_fn);
		return -1;
	}
	paf_rec_t r;
	while (paf_read(fp, &r) >= 0) {
		//primary tag always exist? it is an optional tga
		/*fprintf(stderr, "into %d\n", 1);	*/
		if ((onlyp && r.isprim) || !onlyp) {
			if (r.mq > min_mq || r.qe - r.qs > (uint32_t)(min_mlr * r.ql)) {
				/*fprintf(stderr, "into %d\n", 2);	*/
				int ind = sd_put(ctgs, r.tn, r.tl, 1);	
				/*fprintf(stderr, "insert %d\n", ind);	*/
				ctg_pos_push(d,ind);	
				/*if (r.rev) {*/
					/*uint32_t tmp = r.ts;*/
					/*r.ts = r.tl - r.te;*/
					/*r.te = r.tl - tmp;*/
				/*}*/
				if (r.te - r.ts > 2 * flank) {
					if (r.ts + flank < ctgs->seq[ind].len) r.ts += flank  + 1; //in other words, we only trust the "middle" part one base
					if (r.te > flank) r.te -= flank;
					pos_push(&d->ctg_pos[ind], r.ts << 1 );
					pos_push(&d->ctg_pos[ind], r.te << 1 | 1);
				}
				/*fprintf(stderr, "leave here\n");			*/
			}
		}
	}	
	paf_close(fp);	
	return 0;
}




/*int aa_pb(char *paf_fn, int min_cov, float min_cov_rat, int max_cov, float max_cov_rat, int min_mq, int min_as, uint32_t )*/
int aa_pb(char *paf_fn[], int n_paf, int min_cov, int max_cov, int min_mq, uint32_t flank, float min_mlr, int onlyp, char *out_dir)
{
	sdict_t *ctgs = sd_init();	
#ifdef DEBUG
	fprintf(stderr, "[M::%s] collecting positions from paf file\n", __func__);
#endif
	ctg_pos_t *d = ctg_pos_init();
	if (!d) {
		fprintf(stderr, "[E::%s] fail to initiate space for position arrays\n", __func__);
		return -1;
	}
	int i;
	for (i = 0; i < n_paf; ++i) {
		if(col_pos(paf_fn[i], min_mq, min_mlr, onlyp, ctgs, flank, d)) {
			return -1;	
		}	
	}	
	
#ifdef DEBUG
	fprintf(stderr, "[M::%s] calculating coverage for each base on genome\n", __func__);
#endif
	cov_ary_t *ca = cal_cov(d, ctgs);

	ctg_pos_destroy(d);
	if (!ca) {
		fprintf(stderr, "[W::%s] low quality alignments\n", __func__);
		return 0;	
	}
		
#ifdef DEBUG
	fprintf(stderr, "[M::%s] print coverage histogram for the contigs\n", __func__);
#endif
	/*sel_sup_reg(ca, min_cov_rat, min_cov, max_cov_rat, max_cov, ctgs);*/
	char *type = "PB";
	char *desc = "pacbio data";

	print_coverage_stat(ca, max_cov, ctgs, type, out_dir);
#ifdef DEBUG
	fprintf(stderr, "[M::%s] print coverage for each base of the contigs\n", __func__);
#endif
	print_base_coverage(ca, ctgs, type, out_dir);
#ifdef DEBUG
	fprintf(stderr, "[M::%s] print average coverage for each 1024 base of the contigs\n", __func__);
#endif
	print_coverage_wig(ca, ctgs, type, 1024, out_dir);
#ifdef DEBUG
	fprintf(stderr, "[M::%s] release coverage vector space\n", __func__);
#endif
	cov_ary_destroy(ca, ctgs->n_seq); //a little bit messy
#ifdef DEBUG
	fprintf(stderr, "[M::%s] release sequence vector space\n", __func__);
#endif
	sd_destroy(ctgs);
	fprintf(stderr, "Program finished successfully\n");
	return 0;
}


int main(int argc, char *argv[])
{
	int c;
	int max_cov = 500, min_cov = 10;
	int min_mq = -1;
	uint32_t flank = 0;	
	float min_mlr = 0;
	char *out_dir = ".";
	int onlyp = 1;
	while (~(c=getopt(argc, argv, "O:M:m:q:f:l:ph"))) {
		switch (c) {
			case 'M': 
				max_cov = atoi(optarg);
				break;
			case 'm':
				min_cov = atoi(optarg);
				break;
			case 'l':
				flank = atoi(optarg); 
				break;
			case 'f':
				min_mlr = atof(optarg);
				break;
			case 'q':
				min_mq = atoi(optarg); 
				break;
			case 'p':
				onlyp = 1; 
				break;
			case 'O':
				out_dir = optarg; 
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: aa_pb [options] <PAF_FILE> ...\n");
				fprintf(stderr, "Options:\n");	
				/*fprintf(stderr, "         -m    INT      minimum coverage [5]\n");	*/
				fprintf(stderr, "         -M    INT      maximum coverage [500]\n");
				fprintf(stderr, "         -f    FLOAT    minimum mapping length ratio [0]\n");	
				fprintf(stderr, "         -q    INT      minimum mapping quality [-1]\n");	
				fprintf(stderr, "         -l    INT      flanking space [0]\n");
				fprintf(stderr, "         -p    BOOL     only use primary alignments [TRUE]\n");
				fprintf(stderr, "         -O    STR      output directory [.]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] require at least one paf file!\n", __func__); goto help;
	}
	char **paf_fn = argv+optind;
	int n_paf = argc - optind;
	fprintf(stderr, "Program starts\n");	
	aa_pb(paf_fn, n_paf, min_cov, max_cov, min_mq, flank, min_mlr, onlyp, out_dir);
	return 0;	


}
