/*
 * =====================================================================================
 * *       Filename:  opt.c
 *
 *    Description:  opt.c
 *
 *        Version:  1.0
 *        Created:  02/05/2018 19:51:46
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
#include <string.h>
#include <getopt.h>

#include "opt.h"

int help()
{
	fprintf(stderr, "\nUsage: update [options] <PAF>\n");
	fprintf(stderr, "Options:\n");	
	/*fprintf(stderr, "         -r    FLOAT    minimum overlap ratio for an alignment [0.8]\n");	*/
	
	fprintf(stderr, "         -c    STR      base-level coverage file [NULL]\n");
	fprintf(stderr, "         -T    STR      cutoffs file [NULL]\n");
	fprintf(stderr, "         -f    INT      minimum fraction of haploid/diploid/bad/repetitive bases in a sequence [.8]\n")	;
	fprintf(stderr, "         -a    INT      minimum alignment score [70]\n");
	fprintf(stderr, "         -b    INT      minimum max match score [200]\n");
	fprintf(stderr, "         -2    BOOL     2 rounds chaining [FALSE]\n")	;
	fprintf(stderr, "         -m    INT      minimum matching bases for chaining [500]\n");
	fprintf(stderr, "         -M    INT      maximum gap size for chaining [20K]\n")	;
	fprintf(stderr, "         -G    INT      maximum gap size for 2nd round chaining [50K]\n");
	fprintf(stderr, "         -l    INT      minimum chaining score for a match [10K]\n")	;
	fprintf(stderr, "         -E    INT      maximum extension for contig ends [15K]\n")	;
	/*fprintf(stderr, "         -r    BOOL     read to reference alignment [FALSE]\n")	;*/
	fprintf(stderr, "         -h             help\n")	;
	return 0;
}

int parse_args(int argc, char *argv[], opt *o)
{
	o->max_gs = 20000;
	o->max_gs2rd = 50000;
	o->min_bl = 500; // mislabeled should be ml
	o->max_ext = 15000;
	o->min_dup_bl = 10000;
	o->s2s = 1;
	o->sr = 0;
	o->cov_fn = 0;
	o->ctg_gap = 15000;
	o->cut_fn = 0;
	o->min_frac = .80;
	o->min_bmf = 70;
	o->min_mmf = 200;
	int c;
	while ((c = getopt(argc, argv, "a:b:c:T:G:M:l:E:m:f:r2h")) != -1) {
		switch (c) {
			case 'a':
				o->min_bmf = atoi(optarg);
				break;
			case 'b':
				o->min_mmf = atoi(optarg);
				break;
			case 'c':
				o->cov_fn = optarg;
				break;
			case 'T':
				o->cut_fn = optarg;
				break;
			case 'G':
				o->max_gs2rd = atoi(optarg);
				break;
			case 'M':
				o->max_gs = atoi(optarg);
				break;
			case 'l':
				o->min_dup_bl = atoi(optarg);
				break;
			case 'E':
				o->max_ext = atoi(optarg);
				break;
			case 'm':
				o->min_bl = atoi(optarg);
				break;
			case 'f':
				o->min_frac = atof(optarg);
				break;
			case 'r':
				o->s2s = 0;
				break;
			case '2':
				o->sr = 1;
				break;
			case 'h':
				help();
				return 1;
			default:
				fprintf(stderr,"[E::%s] undefined option %c\n", __func__, c);
				help();
				return 1;
		}
	}	
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] paf file can't be omitted!\n", __func__);
		help();
		return 1;
	} 
	o->paf_fn = argv[optind];
	return 0;
}

