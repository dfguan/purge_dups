/*
 * =====================================================================================
 *
 *       Filename:  cov.c
 *
 *    Description:  read coverage file 
 *
 *        Version:  1.0
 *        Created:  14/03/2019 17:21:48
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
#include <string.h>
#include "cov.h"
#include "kseq.h"
#include "asset.h"
#include "sdict.h"

KSTREAM_INIT(gzFile, gzread, gzseek, 0x10000);

typedef struct {
	char *seqn;
	uint32_t len;
}hdr_t;


int parse_hdr(char *s, int l, hdr_t *hdr)
{
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) hdr->seqn = q;
		else if (t == 1) hdr->len = strtol(q, &r, 10);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 1) return -1;
	return 0;
}

int parse_cov(char *s, int l, cov_t *cv)
{
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) cv->s = strtol(q, &r, 10);
		else if (t == 1) cv->e = strtol(q, &r, 10);
		else if (t == 2) cv->coverage = strtol(q, &r, 10);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 2) return -1;
	return 0;
}

//not good 
void *read_covs(char *fn, void *sn)
{
	sdict_t *osn = (sdict_t *)sn;
	kstream_t *ks;
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	kstring_t buf = {0, 0, 0};	
	cov_ary_t *ca = calloc(16, sizeof(cov_ary_t));
	int ca_n = 0, ca_m = 16, dret;
	hdr_t hd;	
	cov_t cv;
	cov_ary_t *cur = 0;
	while (ks_getuntil(ks, KS_SEP_LINE, &buf, &dret) >= 0) {
		if (buf.s[0] == '>')	{
			if (!parse_hdr(buf.s + 1, buf.l -1, &hd)) {
				int idx = sd_put(osn, hd.seqn, hd.len, 1);	
				if (idx >= ca_m) {
					ca_m = idx << 1;
					cov_ary_t *new_ca = calloc(ca_m, sizeof(cov_ary_t));
					memcpy(new_ca, ca, ca_n * sizeof(cov_ary_t));
					free(ca);
					ca = new_ca;	
				} 
				cur = &ca[idx];	
				cur->len = hd.len;
				if (idx >= ca_n) 
					ca_n = idx + 1;
			}		
		} else {
			parse_cov(buf.s, buf.l, &cv);
			cov_ary_push(cur, cv.s, cv.e, cv.coverage);
		}
	}
	gzclose(fp);
	free(buf.s);
	ks_destroy(ks);
	return ca;
}


