/*
 * =====================================================================================
 *
 *       Filename:  cdict.c
 *
 *    Description:  counter dictionary 
 *
 *        Version:  1.0
 *        Created:  21/10/2018 10:05:56
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

#include "khash.h"

#include "cdict.h"

KHASH_MAP_INIT_STR(str, uint32_t)
typedef khash_t(str) shash_t;

void cd_init(cdict_t *c)
{
	c->h = kh_init(str);
}

void cd_destroy(cdict_t *c)
{
	if (c) {
		uint32_t i;
		if (c->h) kh_destroy(str, (shash_t *)c->h);
		for ( i = 0; i < c->n_cnt; i += 2)  {
			if (c->cnts[i].is_l && c->cnts[i].name) 
				free(c->cnts[i].name);	
		}
		if (c->cnts) free(c->cnts);
	}
}

void cd_set_lim(cdict_t *c, uint32_t n, uint32_t min_wt)
{
	cdict_t *t;
	uint32_t i, j;
	for (i = 0; i < n; ++i)  {
		t = c + i;
		//maximum weight
		/*t->lim = 1;*/
		t->lim = 0;
		/*if (!t->n_cnt)  continue;*/
		/*uint32_t max = 5; */
			/*= t->cnts[0].cnt;*/
		/*fprintf(stderr, "%u\t", max);*/
		/*if (max <= 2) continue;*/
		/*max >>= 1;*/
		for(j = 0; j < t->n_cnt; ++j) {
			if (t->cnts[j].cnt >= min_wt)
				++t->lim;
			else
				break;				
		}	
	}
}

void cd_add2(cdict_t *c, const char *name, uint32_t is_l, uint32_t cnt, uint32_t snp_n)
{
	shash_t *h = (shash_t *)c->h;
	khint_t k;
	int absent;
	/*if (h) fprintf(stderr, "cd add");*/
	k = kh_put(str, h, name, &absent);
	if (absent) {
		if (c->n_cnt == c->m_cnt) {
			c->m_cnt = c->m_cnt ? c->m_cnt << 1 : 16;
			cd_cnt_t *ncnts = calloc(c->m_cnt, sizeof(cd_cnt_t));	
			if (c->cnts) memcpy(ncnts, c->cnts, sizeof(cd_cnt_t) * c->n_cnt);
			if (c->cnts) free(c->cnts);
			c->cnts = ncnts;
		}
		kh_key(h, k) = c->cnts[c->n_cnt].name = strdup(name); //
		kh_val(h, k) = c->n_cnt;  
		c->cnts[c->n_cnt].is_l = 0;
		
		c->cnts[c->n_cnt+1].name = c->cnts[c->n_cnt].name; // init two 
		c->cnts[c->n_cnt+1].is_l = 1;	
		
		c->cnts[c->n_cnt | is_l].cnt = cnt;
		c->cnts[c->n_cnt | is_l].snp_n = snp_n;
		c->n_cnt += 2;
	} else {
		uint32_t ind = kh_val(h, k);
		/*fprintf(stderr, "%u\t %s exist\n", ind, name);*/
		c->cnts[ind | is_l].cnt = cnt;
		c->cnts[ind | is_l].snp_n = snp_n;
	}
		
	/*fprintf(stderr, "%s\t%s\t%d\t%d\n", name, c->cnts[kh_val(h,k)].name, kh_val(h, k), absent);*/
	/*if (absent) {*/
		/*sd_seq_t *s;*/

		/*if (d->n_seq == d->m_seq) {*/
			/*d->m_seq = d->m_seq? d->m_seq<<1 : 16;*/
			/*d->seq = (sd_seq_t*)realloc(d->seq, d->m_seq * sizeof(sd_seq_t));*/
		/*}*/
		/*s = &d->seq[d->n_seq];*/
		/*kh_key(h, k) = s->name = strdup(name);*/
		/*kh_val(h, k) = d->n_seq++;*/
	/*} // TODO: test if len is the same;*/
	/*return kh_val(h, k);*/
}
/* 
void cd_add2(cdict_t *c, const char *name, uint32_t is_l, uint32_t cnt)
{
	shash_t *h = (shash_t *)c->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, name, &absent);
	if (absent) {
		if (c->n_cnt == c->m_cnt) {
			c->m_cnt = c->m_cnt ? c->m_cnt << 1 : 16;
			cd_cnt_t *ncnts = calloc(c->m_cnt, sizeof(cd_cnt_t));	
			if (c->cnts) memcpy(ncnts, c->cnts, sizeof(cd_cnt_t) * c->n_cnt);
			if (c->cnts) free(c->cnts);
			c->cnts = ncnts;
		}
		kh_key(h, k) = c->cnts[c->n_cnt].name = strdup(name); //
		kh_val(h, k) = c->n_cnt;  
		c->cnts[c->n_cnt].is_l = 0;
		
		c->cnts[c->n_cnt+1].name = c->cnts[c->n_cnt].name; // init two 
		c->cnts[c->n_cnt+1].is_l = 1;	
		
		c->cnts[c->n_cnt | is_l].cnt = cnt;
		c->cnts[c->n_cnt | is_l].snp_n = 1;
		c->n_cnt += 2;
	} else {
		uint32_t ind = kh_val(h, k);
		c->cnts[ind | is_l].cnt = cnt;
		c->cnts[ind | is_l].snp_n = 1;
	}
} 
*/
void cd_add(cdict_t *c, const char *name, uint32_t is_l, uint32_t snp_n)
{
	shash_t *h = (shash_t *)c->h;
	khint_t k;
	int absent;
	/*if (h) fprintf(stderr, "cd add");*/
	k = kh_put(str, h, name, &absent);
	if (absent) {
		/*fprintf(stderr, "%u\n", c->n_cnt);*/
		if (c->n_cnt == c->m_cnt) {
			c->m_cnt = c->m_cnt ? c->m_cnt << 1 : 16;
			cd_cnt_t *ncnts = calloc(c->m_cnt, sizeof(cd_cnt_t));	
			if (c->cnts) memcpy(ncnts, c->cnts, sizeof(cd_cnt_t) * c->n_cnt);
			if (c->cnts) free(c->cnts);
			c->cnts = ncnts;
		}
		kh_key(h, k) = c->cnts[c->n_cnt].name = strdup(name); //
		kh_val(h, k) = c->n_cnt; // necessary ? 
		c->cnts[c->n_cnt].is_l = 0;
		
		c->cnts[c->n_cnt+1].name = c->cnts[c->n_cnt].name; // init two 
		c->cnts[c->n_cnt+1].is_l = 1;	
		
		++c->cnts[c->n_cnt | is_l].cnt;
		c->cnts[c->n_cnt | is_l].snp_n = snp_n;
		c->n_cnt += 2;
	} else {
		uint32_t ind = kh_val(h, k);
		++c->cnts[ind | is_l].cnt;
		c->cnts[ind | is_l].snp_n = snp_n;
	}
		
	/*if (absent) {*/
		/*sd_seq_t *s;*/

		/*if (d->n_seq == d->m_seq) {*/
			/*d->m_seq = d->m_seq? d->m_seq<<1 : 16;*/
			/*d->seq = (sd_seq_t*)realloc(d->seq, d->m_seq * sizeof(sd_seq_t));*/
		/*}*/
		/*s = &d->seq[d->n_seq];*/
		/*kh_key(h, k) = s->name = strdup(name);*/
		/*kh_val(h, k) = d->n_seq++;*/
	/*} // TODO: test if len is the same;*/
	/*return kh_val(h, k);*/
}

void cd_norm(cdict_t *c)
{
	if (c) {
		uint32_t i;
		cd_cnt_t *t = c->cnts;
		for ( i = 0; i < c->n_cnt; ++i) 
			if (t[i].snp_n) t[i].cnt = t[i].cnt * 50000 / t[i].snp_n;
			else t[i].cnt = 0;	
		
	}	
}

int cmp_cnt(const void *a, const void *b) 
{
	cd_cnt_t *f = (cd_cnt_t *)a;
	cd_cnt_t *h = (cd_cnt_t *)b;
	if (f->cnt > h->cnt) return -1;
	else if (f->cnt == h->cnt) {
		if (f->snp_n > h->snp_n) return 1;
		else if (f->snp_n == h->snp_n) return 0;
		else return -1;
	}
	else return 1;
}

void cd_sort(cdict_t *c)
{
	if (c) {
		qsort(c->cnts, c->n_cnt, sizeof(cd_cnt_t), cmp_cnt);
	}	
}



