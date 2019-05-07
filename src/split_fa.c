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

int split_fa(char *fn)
{
	gzFile fp;
	kseq_t *seq;
	fp = strcmp(fn, "-") == 0 ? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
	if (!fp) return 1;
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) 
		split(seq->name.s, seq->seq.s, seq->seq.l);
 	//add some basic statistics maybe 		
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}



int main(int argc, char *argv[])
{

	char *program, *fafn;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	if (argc < 2) {
		fprintf(stderr, "\n  Usage: %s <FA/FA.GZ>\n", program);
		return 1;
	}
	fafn = argv[1];
	return split_fa(fafn);
}
