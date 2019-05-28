/*
 * =====================================================================================
 *
 *       Filename:  split_fa.c
 *
 *    Description:  split sequences by 'N's in fasta file 
 *
 *        Version:  1.0
 *        Created:  14/03/2019 08:27:05
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread, gzseek)

int split(char *name, char *seq, uint32_t l)
{
	uint32_t s =0, i = 0;
	while ( i <= l) {
		if (i == l || seq[i] == 'N' || seq[i] == 'n') {
			if (i != l) seq[i] = 0;
			if (i - s) 
				fprintf(stdout, ">%s:%u-%u\n%s\n", name, s + 1, i, seq + s);
			uint32_t j = i + 1;
			for ( j = i + 1; j < l && (seq[j] == 'N' || seq[j] == 'n'); ++j);	
			s = i = j;	
		}
		++i;	
	}	
	return 0;
}

int split_fa(char *fn, int split_by_n)
{
	gzFile fp;
	kseq_t *seq;
	fp = strcmp(fn, "-") == 0 ? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
	if (!fp) return 1;
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) 
		split_by_n ? split(seq->name.s, seq->seq.s, seq->seq.l) : fprintf(stdout, ">%s:%u-%u\n%s\n", seq->name.s, 1, seq->seq.l, seq->seq.s); // notice if larger than 4G error
 	//add some basic statistics maybe 		
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}



int main(int argc, char *argv[])
{

	char *program, *fafn; 
	int c;
	int split_by_n = 1;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	while (~(c=getopt(argc, argv, "nh"))) {
		switch (c) {
			case 'n': 
				split_by_n = 0;
				break;
			default:
help:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
				fprintf(stderr, "\nUsage: %s  [<options>] <STAT> ...\n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -n    BOOL    block split by N\n");	
				fprintf(stderr, "         -h            help\n");
				fprintf(stderr, "\n\nNotice: please set \"-n\" if you do not want break your scaffols into contigs.\n");
				return 1;	
		}		
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "[E::%s] require fa or fa.gz file", __func__); goto help;
	}
	fafn = argv[optind];
	return split_fa(fafn, split_by_n);
}
