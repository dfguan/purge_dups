/*
 * =====================================================================================
 *
 *       Filename:  aa.c
 *
 *    Description:  common used functions 
 *
 *        Version:  1.0
 *        Created:  15/09/2018 19:12:48
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
#include <stdint.h>
#include "asset.h"

#include "ksort.h"
#define pos_key(a) (a)
KRADIX_SORT_INIT(pos, uint32_t, pos_key, 4)
void cov_ary_push(cov_ary_t *c, uint32_t s, uint32_t e, int cov)
{
	if (c->n >= c->m) {
		c->m = c->m ? c->m << 1 : 2;
		c->intv = realloc(c->intv, c->m * sizeof(cov_t));	
	}
	c->intv[c->n++] = (cov_t){s, e, cov};
	return ;
}

void cov_ary_destroy(cov_ary_t *ca, int n)
{
	if (ca && n) {
		int i;
		for ( i = 0; i < n; ++i) 
			if (ca[i].intv) free(ca[i].intv);
		free(ca);	
	}
}

void cord_push(cord_t *c, cors *cord)
{
	if (c->n >= c->m) {
		c->m = c->m ? c->m << 1 : 2;
		c->coords = realloc(c->coords, c->m * sizeof(cors));	
	} 
	c->coords[c->n++] = *cord;
}

void cord_push1(cord_t *c, cors *cord, char *ann)
{
	if (c->n >= c->m) {
		c->m = c->m ? c->m << 1 : 2;
		c->coords = realloc(c->coords, c->m * sizeof(cors));	
		c->anno = realloc(c->anno, c->m * sizeof(char *));	
	} 
	c->coords[c->n] = *cord;
	c->anno[c->n++] = ann;
}
void cord_destroy(cord_t *c, int n)
{
	if (c) {
		int i;
		for ( i = 0; i < n; ++i ) {
			if (c[i].coords) free(c[i].coords);
			if (c[i].anno) {
				int j ;
				for (j=0; j < c[i].n; ++j) free(c[i].anno[j]);
				free(c[i].anno);
			}	
		}
		free(c);	
	}
}

ctg_pos_t *ctg_pos_init()
{
	return (ctg_pos_t*)calloc(1, sizeof(ctg_pos_t));
}

void ctg_pos_push(ctg_pos_t *_d,int s)
{
	if (s >= _d->n) { //problem occurs when not add by order // this is a new one
		if (_d->n >= _d->m){
			_d->m = _d->m ?  _d->m << 1 : 16;
			/*_d->m = s ? 2 : s << 1; //in case of jumping */
			_d->ctg_pos = realloc(_d->ctg_pos, sizeof(pos_t) * _d->m);
			/*if (!_d->ctg_pos)  fprintf(stderr, "can");*/
		} //should calloc memcpy by yourself
		_d->ctg_pos[s] = (pos_t){0,0,0}; 
		++_d->n;
	} 
	return;
}



void pos_push(pos_t *ps, uint32_t _p)
{
	if (ps->n >= ps->m) {
		ps->m = ps->m ? ps->m << 1 : 2;
		ps->p = realloc(ps->p, sizeof(uint32_t) * ps->m);
	}
	ps->p[ps->n++] = _p;
}
void pos_destroy(pos_t *ps)
{
	if (ps && ps->p) free(ps->p);
}

void ctg_pos_destroy(ctg_pos_t *_d)
{
	if (_d && _d->ctg_pos) {
		int i;
		for ( i = 0; i < _d->n; ++i) {
			pos_destroy(&_d->ctg_pos[i]);
		}
		free(_d->ctg_pos);	
		free(_d);
	}
}

/*char *strjoin(...)*/
/*{*/
	/*va_start()*/

/*}*/

/*void dump_coverage(cov_ary_t *ca, sdict_t *ctgs, char *tp)*/
/*{*/
	

/*}*/

void proc_covs(uint32_t s, uint32_t e, int *covs, int w, int freq, uint32_t *w_idx, int coverage, FILE *fp, char *name)
{
	uint32_t _w_idx = *w_idx;
	int col_len = w / freq;
	for (; s < e; ++s) {
		if (s / w > _w_idx) {
			int j;
			fprintf(fp, "%s\t%u",name, _w_idx);
			for ( j = 0; j < freq; ++j) 
				fprintf(fp, "\t%u", covs[j]);
			fprintf(fp, "\n");
			memset(covs, 0, sizeof(uint32_t) * freq);
			++_w_idx;	
		}
		uint32_t covs_idx = (s - _w_idx * w) / col_len;
		if (covs[covs_idx] < coverage) covs[covs_idx] = coverage;
	}		
	*w_idx = _w_idx;
}


void print_maxcol_cov(cov_ary_t *ca, sdict_t* ctgs, char *tp, int w, int freq, char *out_dir)
{
	char *col = malloc(strlen(tp) + strlen(out_dir) + 10);
	strcpy(col, out_dir);
	strcat(col, "/");
	strcat(col, tp);
	strcat(col, ".cov.col");
	FILE *fp = fopen(col, "w");
	if (!fp) return;
	uint32_t i;
	int col_len = w / freq;
	int* covs = (int *)calloc(freq, sizeof(int));
	for ( i = 0; i < ctgs->n_seq; ++i) {
		if (ca[i].n) {
			uint32_t j, w_idx = 0;
			for ( j = 0; j < ca[i].n; ++j) {
				/*if (j == 0) s = 0,e = ca[i].intv[j].s;*/
				/*else s = ca[i].intv[j-1].e, e = ca[i].intv[j].s;*/
					
				/*proc_covs(s, e, covs, w, freq, &w_idx, 0, fp, ctgs->seq[i].name);				   */
				proc_covs(ca[i].intv[j].s, ca[i].intv[j].e, covs, w, freq, &w_idx, ca[i].intv[j].coverage, fp, ctgs->seq[i].name);	
				if (j == ca[i].n - 1) {
					uint32_t s, e;
					s = ca[i].intv[j].e, e = ctgs->seq[i].len;
					proc_covs(s, e, covs, w, freq, &w_idx, 0, fp, ctgs->seq[i].name);
					e = (e - w_idx * w + col_len - 1) / col_len;
					if (e) fprintf(fp, "%s\t%u", ctgs->seq[i].name, w_idx);
					for (s = 0; s < e; ++s)
					   fprintf(fp, "\t%u", covs[s]);
					fprintf(fp, "\n");	
				}
			}		
		}	
	}
}

void print_coverage(cov_ary_t *ca, sdict_t* ctgs, char *tp)
{

	char *bdgname = malloc(strlen(tp) + 9);
	strcpy(bdgname, tp);
	strcat(bdgname, ".cov.bedg");
	FILE *fp = fopen(bdgname, "w");
	if (!fp) return;
	fprintf(fp, "track type=\"bedGraph\" name=\"%s\"\n", tp);	

	/*FILE *fp_cvg = fopen()*/
	int i, j;
	for ( i = 0; i < ctgs->n_seq; ++i) 
		if (ca[i].n) 
			for (j = 0; j < ca[i].n; ++j)  
				fprintf(fp,"%s\t%u\t%u\t%d\n", ctgs->seq[i].name, ca[i].intv[j].s, ca[i].intv[j].e, ca[i].intv[j].coverage);
	fclose(fp);	
	free(bdgname);	

}


void print_base_coverage(cov_ary_t *ca, sdict_t* ctgs, char *tp, char *out_dir)
{
	char *wigname = malloc(strlen(tp) + strlen(out_dir) + 12);
	strcpy(wigname, out_dir);
	strcat(wigname, "/");
	strcat(wigname, tp);
	strcat(wigname, ".base.cov");
	FILE *fp = fopen(wigname, "w");
	if (!fp) return;
	int i, j;
	uint32_t s, e;
	for ( i = 0; i < ctgs->n_seq; ++i) {
		if (ca[i].n) {
			fprintf(fp, ">%s\t%u\n", ctgs->seq[i].name, ctgs->seq[i].len);
			for (j = 0; j < ca[i].n; ++j) {
				/*if (j == 0) s = 0, e = ca[i].intv[j].s;*/
				/*else s = ca[i].intv[j-1].e, e = ca[i].intv[j].s;*/
				/*fprintf(fp, "%u\t%u\t%d\n", s, e-1, 0);	*/
				int coverage = ca[i].intv[j].coverage;
				s = ca[i].intv[j].s;
				e = ca[i].intv[j].e;
				fprintf(fp, "%u\t%u\t%d\n", s,e, coverage);	
				if (j == ca[i].n - 1) {
					s = ca[i].intv[j].e, e = ctgs->seq[i].len;
					if (s!=e) fprintf(fp, "%u\t%u\t%d\n",s+1,e, 0);	
				}
			}
		} else {
			fprintf(fp, ">%s\t%u\n", ctgs->seq[i].name, ctgs->seq[i].len);
			fprintf(fp, "%u\t%u\t%d\n", 1, ctgs->seq[i].len, 0);
		}
	}
	fclose(fp);	
	free(wigname);	

}

void print_coverage_stat(cov_ary_t *ca, int max_cov, sdict_t* ctgs, char *tp, char *out_dir)
{
	char *wigname = malloc(strlen(tp) + strlen(out_dir) + 10);
	strcpy(wigname, out_dir);
	strcat(wigname, "/");
	strcat(wigname, tp);
	strcat(wigname, ".stat");
	FILE *fp = fopen(wigname, "w");
	if (!fp) return;
	int i, j;
	uint32_t *freq = (uint32_t *)calloc(max_cov + 1, sizeof(uint32_t));
	if (!freq) {
		fprintf(stderr, "[E::%s] fail to allocate memory space for frequency vector, quit!\n", __func__);
		return;
	}
	for ( i = 0; i < ctgs->n_seq; ++i) {
		if (ca[i].n) {
			for (j = 0; j < ca[i].n; ++j) {
				int coverage = ca[i].intv[j].coverage;
				uint32_t s = ca[i].intv[j].s;
				uint32_t e = ca[i].intv[j].e;
				if (coverage >= max_cov) coverage = max_cov;
				freq[coverage] += (e-s+1);
			}
		}
	}
	for ( i = 0; i <= max_cov; ++i) 
		fprintf(fp, "%d\t%u\n", i, freq[i]);
	fclose(fp);	
	/*fprintf(stderr, "releasing frequency vector\n");*/
	free(freq);
	/*fprintf(stderr, "releasing filename vector\n");*/
	free(wigname);	
	/*fprintf(stderr, "leave %s\n", __func__);*/
}

void print_coverage_wig(cov_ary_t *ca, sdict_t* ctgs, char *tp, uint32_t ws, char *out_dir)
{
	char *wigname = malloc(strlen(tp) + strlen(out_dir) + 10);
	strcpy(wigname, out_dir);
	strcat(wigname, "/");
	strcat(wigname, tp);
	strcat(wigname, ".cov.wig");
	FILE *fp = fopen(wigname, "w");
	if (!fp) return;
	fprintf(fp, "track type=\"wiggle_0\" name=\"%s\"\n", tp);	

	int i, j;
	uint64_t *total_coverage = NULL;
	uint32_t n, m;
	n = m = 0;
	for ( i = 0; i < ctgs->n_seq; ++i) {
		if (ca[i].n) {
			n = 0;
			uint32_t tot_cov_cnt = (ctgs->seq[i].len + ws - 1)/ws;	
			/*fprintf(stderr, "%s\t%d\t%d\n", ctgs->seq[i].name, ca[i].intv[ca[i].n-1].e, m);*/
			if (tot_cov_cnt > m) {
				if (total_coverage) free(total_coverage);
				total_coverage = calloc(tot_cov_cnt, sizeof(uint64_t));
				m = tot_cov_cnt;
			} else memset(total_coverage, 0, sizeof(uint64_t) * tot_cov_cnt);
			for (j = 0; j < ca[i].n; ++j) {
				uint32_t s_idx, e_idx;
				uint32_t s, e;
				int coverage = ca[i].intv[j].coverage;
				s = ca[i].intv[j].s - 1;
				e = ca[i].intv[j].e - 1;
				s_idx = s / ws;			
				e_idx = e / ws;
				if (s_idx == e_idx) 
					total_coverage[s_idx] += (uint64_t)(e - s + 1) * coverage;
				else {
					total_coverage[s_idx] += (uint64_t)((s_idx + 1) * ws - s) * coverage;
					for (++s_idx; s_idx < e_idx; ++s_idx) total_coverage[s_idx] += (uint64_t) ws * coverage;
					total_coverage[s_idx] += (uint64_t)(e - s_idx * ws + 1)*coverage;	
				}	
			}
			if (tot_cov_cnt > 2) {
				fprintf(fp, "fixedStep chrom=%s start=1 step=%d span=%d\n", ctgs->seq[i].name, ws, ws);
				uint32_t z;
				for (z = 0; z < tot_cov_cnt - 2; ++z) fprintf(fp, "%u\n", total_coverage[z]/ws);
				uint32_t lastnbases = ctgs->seq[i].len - (tot_cov_cnt - 2) * ws;
				fprintf(fp, "variableStep chrom=%s span=%d\n", ctgs->seq[i].name, lastnbases);
				fprintf(fp, "%u %u\n", (tot_cov_cnt - 2) * ws + 1, (total_coverage[z] + total_coverage[z+1]) /lastnbases);
			} else {
				uint32_t lastnbases = ctgs->seq[i].len;
				fprintf(fp, "variableStep chrom=%s span=%d\n", ctgs->seq[i].name, lastnbases);
				fprintf(fp, "%u %u\n", 1, total_coverage[0] /lastnbases);
			} 
		}
	}
	fclose(fp);	
	free(wigname);	
}

void sel_sup_reg(cov_ary_t *ca, int min_cov, int max_cov, sdict_t* ctgs, char *tp, char *desc)
{
	//print a header 
	fprintf(stdout, "track name=\"%s\" description=\"%s\"\n", tp,desc);	

	/*FILE *fp_cvg = fopen()*/
	int i, j;
	for ( i = 0; i < ctgs->n_seq; ++i) {
		if (ca[i].n) {
			/*int set_cov = ca[i].tot_cov/ctgs[i].seq->len;	*/
			int set_min_cov = min_cov;
			int set_max_cov = max_cov;
			/*int set_min_cov = min(min_cov_rat *set_cov, min_cov);*/
			/*int set_max_cov = max(max_cov_rat * set_cov, max_cov);*/
			uint32_t pe, ps; 
			uint8_t is_set = 0;
			for (j = 0; j < ca[i].n; ++j)  {
				/*fprintf(stderr,"%s\t%u\t%u\t%d\n", ctgs->seq[i].name, ca[i].intv[j].s, ca[i].intv[j].e, ca[i].intv[j].coverage);*/
				if (!(ca[i].intv[j].coverage < set_min_cov || ca[i].intv[j].coverage > set_max_cov)) {
					if (!is_set) {//need to look forward to see if it can go furthur
						ps = ca[i].intv[j].s;
						pe = ca[i].intv[j].e;
						is_set = 1;
					} else if (ca[i].intv[j].s - pe < CONT_THRES) {
						pe = ca[i].intv[j].e;	
					} else {
						fprintf(stdout, "%s\t%u\t%u\n", ctgs->seq[i].name, ps, pe);
						ps = ca[i].intv[j].s;
						pe = ca[i].intv[j].e;	
					}
		
				}
			}
			if (is_set) fprintf(stdout, "%s\t%u\t%u\n", ctgs->seq[i].name, ps, pe);
		}
	}
}

cov_ary_t *cal_cov(ctg_pos_t *d, sdict_t* ctgs)
{
	cov_ary_t *ca = calloc(ctgs->n_seq, sizeof(cov_ary_t));
	int i,j;
	for (i = 0; i < d->n; ++i) {
		ca[i].len = ctgs->seq[i].len;
		long tot_cov =0, cov = 0;
		//sort first
		/*fprintf(stderr, "enter sort\n");*/
		pos_t *ps = &d->ctg_pos[i];
		radix_sort_pos(ps->p, ps->p + ps->n);
		/*fprintf(stderr, "leave sort\n");*/
		/*fprintf(stdout, "%u\n", ps->n);*/
		int s = 1, e = 0;
		/*fprintf(stderr, "enter cal\n");*/
		for (j = 0; j < ps->n; ++j) {
			if (ps->p[j] & 1) {
				e = ps->p[j] >> 1;
				if (e >= s) {
					cov_ary_push(&ca[i], s, e, cov);
					tot_cov = (s -e + 1) * cov;
				}		
				s = (ps->p[j] >> 1) + 1;	
				--cov;	
			} else {
				e = (ps->p[j] >> 1) - 1;
				if (e >= s) {
					cov_ary_push(&ca[i], s, e, cov);
					tot_cov = (s -e + 1) * cov;
				}		
				s = ps->p[j] >> 1;	
				++cov;	
			}
		}
		if (e < ca[i].len) {
			e = ca[i].len;
			cov_ary_push(&ca[i], s, e, cov);
		}
		/*fprintf(stderr, "leave cal\n");*/
		ca[i].tot_cov = tot_cov;	
	}
	return ca;
}
