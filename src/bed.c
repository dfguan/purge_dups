/*
 * =====================================================================================
 *
 *       Filename:  bed.c
 *
 *    Description:  realization of bed functions
 *
 *        Version:  1.0
 *        Created:  27/09/2018 11:49:20
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
#include "bed.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, gzseek, 0x10000)

bed_file_t *bed_open(const char *fn)
{
	kstream_t *ks;
	gzFile fp;
	bed_file_t *pf;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	pf = (bed_file_t*)calloc(1, sizeof(bed_file_t));
	pf->fp = ks;
	pf->is_eof = 0;
	return pf;
}

int bed_close(bed_file_t *pf)
{
	kstream_t *ks;
	if (pf == 0) return 0;
	free(pf->buf.s);
	ks = (kstream_t*)pf->fp;
	gzclose(ks->f);
	ks_destroy(ks);
	free(pf);
	return 0;
}

int bed_init_track(char *s, int sl, bed_hdr_t *hdr)
{
	int i;
	char *q = s;

	for ( i = 0; i <= sl; ++i) {
		if (i < sl && (s[i] != '\t' && s[i] != ' ')) continue;
		s[i] = 0;
		if (!strncmp(q, "name", 4)) {
			q += 5;
			if (q[0] == '"') { //remove semicolon
				s[i-1] = 0;//user's responsiblity to check if semicolon is complete
				++q; 
			}	
			hdr->type = strdup(q);		
		} else if (!strncmp(q, "description", 11)) {
			q += 4;
			if (q[5] == '"') { //remove semicolon
				s[i-1] = 0;//user's responsiblity to check if semicolon is complete
				++q; 
			}	
			hdr->desc = strdup(q);		
		}
		q = s + i + 1;
	}	
	return 0;
}

int bed_hdr_read(bed_file_t *pf, bed_hdr_t *hdr)
{
	int ret, dret;	
	long int header_size = 0;	
	while((ret = ks_getuntil((kstream_t*)pf->fp, KS_SEP_LINE, &pf->buf, &dret)) >= 0) {
	if (!strncmp(pf->buf.s, "track", 5)) {
		header_size += pf->buf.l + 1;
		bed_init_track(pf->buf.s, pf->buf.l, hdr);
	} else if (strncmp(pf->buf.s, "browser", 7)) 
		break;
	else
		header_size += pf->buf.l + 1;
	}
	ks_seek((kstream_t *)pf->fp, header_size, SEEK_SET);//reset header 
	if (ret < 0) pf->is_eof = 1;
	return 0;
}

int bed_parse_lnks(int l, char *s, lnk_rec_t *pr) // s must be NULL terminated
{ // on return: <0 for failure; 0 for success; >0 for filtered
	char *q, *r;
	int i, t;
	
	/*pr->tech = NULL;*/
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) pr->ctgn = q;
		else if (t == 1) pr->is_l = *q == '+'?1:0; //1: left 0: right
		else if (t == 2) pr->ctgn2 = q;
		else if (t == 3) pr->is_l2 = *q == '+'?1:0;
		else if (t == 4) pr->wt = strtol(q, &r, 10);
		else if (t == 5) pr->llen = strtol(q, &r, 10);
		else if (t == 6) pr->rlen = strtol(q, &r, 10);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 4) return -1;
	return 0;
}
int bed_parse(int l, char *s, bed_rec_t *pr) // s must be NULL terminated
{ // on return: <0 for failure; 0 for success; >0 for filtered
	char *q, *r;
	int i, t;
	
	/*pr->tech = NULL;*/
	pr->le = pr->l_snp_n = pr->rs = pr->r_snp_n = 0;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) pr->ctgn = q;
		else if (t == 1) pr->len = strtol(q, &r, 10);
		else if (t == 2) pr->le = strtol(q, &r, 10);
		else if (t == 3) pr->rs = strtol(q, &r, 10);
		else if (t == 4) pr->l_snp_n = strtol(q, &r, 10);
		else if (t == 5) pr->r_snp_n = strtol(q, &r, 10);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 1) return -1;
	return 0;
}

int lnk_read(bed_file_t *pf, lnk_rec_t *r)
{
		int ret, dret;
read_more:
		ret = ks_getuntil((kstream_t*)pf->fp, KS_SEP_LINE, &pf->buf, &dret);
		if (ret < 0) return ret;

		ret = bed_parse_lnks(pf->buf.l, pf->buf.s, r);
		if (ret < 0)  goto read_more;
		return ret; 
}
int bed_read(bed_file_t *pf, bed_rec_t *r)
{
		int ret, dret;
read_more:
		ret = ks_getuntil((kstream_t*)pf->fp, KS_SEP_LINE, &pf->buf, &dret);
		if (ret < 0) return ret;

		ret = bed_parse(pf->buf.l, pf->buf.s, r);
		if (ret < 0)  goto read_more;
		return ret; 
}

