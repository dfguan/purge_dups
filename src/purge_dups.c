/*
 * =====================================================================================
 *
 *       Filename:  eg.c
 *
 *    Description:  enrich GFA by using miniasm assembly
 *
 *        Version:  1.0
 *        Created:  02/05/2018 19:40:57
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
/*#include <omp.h>*/

#include "kvec.h"
#include "ksort.h"
#include "paf.h"
#include "cov.h"
#include "sdict.h"
#include "asset.h"
#include "opt.h"

#define eg_idx_qn(a, b) (((a).qns>>32) < ((b).qns>>32))
#define eg_idx_tn(a, b) (((a).tns >> 32) < (((b).tns >> 32)))
/*#define INT_MATCH 1*/
/*#define OVLP_MATCH 2*/
/*#define CON_MATCH 3*/
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))




typedef struct {
	uint64_t qns, tns;
	uint32_t qe, te;
	uint32_t qcov:30, qtg:2, tcov:30, ttg:2; //0:S 1:M 2:E 3:ALL
	uint32_t ql, tl;
	uint32_t ml:30, rev:1, del:1;
	uint32_t bl:30, tail:1, con:1; //mapped at the end of reference, contained by another query
	uint8_t  qtp:4, ttp:4; // 0:low,1:hap,2:dip,3:high //take 4 bytes due arrangements don't know how to solve, maybe not a problem.
}eg_hit_t;

#define hit_key(a) ((a).qns)
KRADIX_SORT_INIT(hit, eg_hit_t, hit_key, 8)

	enum dup_type {JUNK, HAPLOTIG, PRIMARY, REPEAT, OVLP, UNKNOWN};
char *dup_type_s[] = {"JUNK", "HAPLOTIG", "PRIMARY", "REPEAT", "OVLP", "UNKNOWN"};
typedef struct {
	uint32_t sn:28, tp:3, del:1; //don't think there will be 2G contigs
	uint32_t s, e;
	uint32_t bst_sn;
}dup_t;

typedef struct {size_t n, m; dup_t *a;} dup_v;
typedef struct {size_t n, m; eg_hit_t *a;} eg_hit_v;
typedef struct {size_t n, m; cov_ary_t *a;} cov_ary_v;
typedef struct {size_t n, m; eg_hit_t *rht; uint64_t *idx;} eg_hits_t;

/**
 * @func	print_hit
 * @brief	print single hit
 *
 */


/*int print_hit(eg_hit_t *r, gfa_t *gf)*/
/*{*/
	/*size_t i = 0; eg_hit_t * rht = r; fprintf(stderr, "%s\t%u\t%u\t%u\t%llu\t%u\t%u\t%u\t%d\t%d\n", gf->seg[rht[i].qns >> 32].name, rht[i].ql, (uint32_t)rht[i].qns, rht[i].qe, rht[i].tns>>32 , rht[i].tl, (uint32_t) rht[i].tns , rht[i].te, rht[i].tail, rht[i].con);*/
	/*return 0;*/
/*}*/
int print_dups(dup_v *dups, sdict_t *dup_n)
{
	dup_t *dp = dups->a;
	size_t n = dups->n;
	size_t i;
	for ( i = 0; i < n; ++i ) 
		fprintf(stdout, "%s\t%u\t%u\t%s\t%s\n", dup_n->seq[dp[i].sn].name, dp[i].s, dp[i].e, dup_type_s[dp[i].tp], dp[i].bst_sn != 0xFFFFFFFF ? dup_n->seq[dp[i].bst_sn].name:"*");
	return 0;
}

//print a single hit
int print_hit(eg_hit_t *rht, sdict_t *tn)
{
	size_t i = 0;
	/*fprintf(stderr, "Order: QID\tQLEN\tQS\tQE\tTID\tTLEN\tTS\tTE\n");*/
	fprintf(stderr, "%s\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%d\t%d\n", tn->seq[rht[i].qns >> 32].name, rht[i].ql, (uint32_t)rht[i].qns, rht[i].qe, rht[i].rev?'-':'+', tn->seq[rht[i].tns>>32].name, rht[i].tl, (uint32_t) rht[i].tns , rht[i].te, rht[i].tail, rht[i].con);
	return 0;
}
/**
 * @func   print_hits
 * @brief  print
 */

int print_hits(eg_hit_t *rht, size_t s, size_t e, sdict_t *sn, char *tag)
{
	size_t i;
	/*fprintf(stderr, "Order: QID\tQLEN\tQS\tQE\tTID\tTLEN\tTS\tTE\n");*/
	fprintf(stderr, "%u hits\n", e - s);
	for (i = s; i < e; ++i) { 
		/*if (!rht[i].del) fprintf(stdout, "%u: %s\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%d\t%d\n", i, tn->seq[rht[i].qns >> 32].name, rht[i].ql, (uint32_t)rht[i].qns, rht[i].qe, rht[i].rev?'-':'+', tn->seq[rht[i].tns>>32].name, rht[i].tl, (uint32_t) rht[i].tns , rht[i].te, rht[i].tail, rht[i].con);*/
		if (!rht[i].del) fprintf(stderr, "%s: %s\t%u\t%u\t%u\t%u\t%u\t%c\t%s\t%u\t%u\t%u\t%u\t%u\t%d\t%d\n",tag, sn->seq[rht[i].qns >> 32].name, rht[i].qtg, rht[i].ql, (uint32_t)rht[i].qns, rht[i].qe, rht[i].qcov, rht[i].rev?'-':'+', sn->seq[rht[i].tns>>32].name, rht[i].ttg, rht[i].tl, (uint32_t) rht[i].tns , rht[i].te, rht[i].tcov, rht[i].bl, rht[i].ml);
	}
	return 0;
}


/**
 * @func   f_qn
 * @brief  return qn 
 *
 *
 */


uint32_t f_qn (eg_hit_t *r)
{
	return r->qns >> 32;
}
/**
 * @func   f_tn
 * @brief  return tn 
 */
uint32_t f_tn (eg_hit_t *r)
{
	return r->tns >> 32;
}	

/**
 * @func  index
 * @brief index alignments according to fun
 */

uint64_t *hit_index(eg_hit_t *rht, size_t n_rht, size_t n_idx, uint32_t (*f)(eg_hit_t *))
{
	uint64_t *idx = (uint64_t *)calloc(n_idx, sizeof(uint64_t)); //don't use malloc here cause some of the n_idx doesn't have value or pass **idx here, 
	if (!idx) {
		fprintf(stderr, "[E::%s] failed to allocate memory space, required %lu\n", __func__, n_idx * sizeof(uint64_t));
		exit(1);
	}
	
	/*fprintf(stderr, "%p\n",idx);*/
	size_t i, last;
	n_idx = 0;
	for (i=1, last = 0; i <= n_rht; ++i) {
		if (i == n_rht || f(rht+i) != f(rht + last)) {
			/*fprintf(stderr, "%u\t%u\n", i, last);*/
			idx[f(rht+last)] = (uint64_t)last << 32 | (i - last); //don't forget the left bracket here !!!
			/*fprintf(stderr, "%u\n", f(rht+i-1));	*/
			last = i;
			/*++*n_idx;*/
		}
	}
	/*for(i = 0; i <*n_idx; ++i)	fprintf(stderr, "this:%u\t%u\t%u\n", i, idx[i]>>32, (uint32_t)idx[i]);*/
	return idx;
}

/**
 * @func    cmp_qtn
 * @brief   compare query and target's names then query start position
 *
 */


int cmp_qtn (const void *r, const void *s) 
{
	eg_hit_t *p = (eg_hit_t *)r;
	eg_hit_t *q = (eg_hit_t *)s;
	/*uint32_t pl = p->ql;*/
	/*uint32_t ql = q->ql;*/
	/*if (pl < ql) return -1; */
	/*else if (pl > ql) return 1;*/
	/*else {*/
		uint32_t pn = p->qns >> 32;
		uint32_t qn = q->qns >> 32;
		if (pn == qn) {
			uint32_t ptn = p->tns >> 32;
			uint32_t qtn = q->tns >> 32;	
			if (ptn > qtn) return -1;
			else if (ptn == qtn) {
				if (p->qns > q->qns) return 1;
				else if (p->qns == q->qns) return 0;
				else return -1;
			} 
			else return 1;	
		} else if (pn > qn) return -1;
		else return 1;
	/*}*/
}
/**
 * @func    cmp_qtgse
 * @brief   compare query's 
 *
 */

void insertionSort(uint32_t a[], sdict_t *sn) 
{ 
    int i,j; 
	uint32_t key;
	int n = sn->n_seq;
    for (i = 1; i < n; i++) { 
        key = a[i]; 
        j = i - 1; 
 		for ( j = i - 1; j >= 0 && sn->seq[a[j]].len < sn->seq[key].len; --j) a[j+1] = a[j]; 
        a[j + 1] = key; 
    } 
	return ;
} 


int cmp_qtgse (const void *r, const void *s) 
{
	eg_hit_t *p = (eg_hit_t *)r;
	eg_hit_t *q = (eg_hit_t *)s;
	uint32_t pn = p->qns >> 32;
   	uint32_t qn = q->qns >> 32;	
	if (pn == qn) {
		if (p->qtg > q->qtg) return 1;
		else if (p->qtg == q->qtg) {
			uint64_t pse = (p->qns << 32) | p->qe;
			uint64_t qse = (q->qns << 32) | q->qe;		
			if (pse > qse) return 1; 
			else if (pse == qse) return 0;
			else return -1;	
		} else return -1;
	} else if (pn > qn) return 1;
	else return -1;
}
/**
 * @func    cmp_q
 * @brief   compare query's name
 *
 */

int cmp_qn(const void *r, const void *s) 
{
	eg_hit_t *p = (eg_hit_t *)r;
	eg_hit_t *q = (eg_hit_t *)s;
	uint32_t pn = p->qns >> 32;
   	uint32_t qn = q->qns >> 32;	
	if (pn == qn) return 0;
	else if (pn > qn) return 1;
	else return -1;
}

int cmp_q (const void *r, const void *s) 
{
	eg_hit_t *p = (eg_hit_t *)r;
	eg_hit_t *q = (eg_hit_t *)s;
	if (p->qns == q->qns) {
		if (p->qe > q->qe) return -1;
		else if (p->qe == q->qe) return 0;
		else return 1;
	} else if (p->qns > q->qns) return 1;
	else return -1;
}

/**
 * @func    cmp_t
 * @brief   compare target's 
 *
 */
//tn   
int cmp_t(const void *r, const void *s) 
{
	eg_hit_t *p = (eg_hit_t *)r;
	eg_hit_t *q = (eg_hit_t *)s;
	if (p->tns == q->tns) {
		if (p->te > q->te) return -1;
		else if (p->te == q->te) return 0;
		else return 1;
	} else if (p->tns > q->tns) return 1;
	else return -1;
	/*fprintf(stderr, "%d\t",z);*/
}

/**
 * @func  mt_sort
 * @brief mulitple threads sort according to tn,ts,te using openmp
 */

int mt_sort(eg_hit_t *rht, uint64_t *idx, size_t n_ind, int (*cmp) (const void *r, const void *s))
{
	size_t i;
	/*#pragma parallel for number_threads(10) //10 threads*/
	for (i=0; i < n_ind; ++i) qsort(&rht[idx[i]>>32], (uint32_t)idx[i], sizeof(eg_hit_t), cmp);		
	return 0;	
}


/**
 * @func  set_cov 
 * @brief set covered alignments 
 *
 */

/*int set_cov(rg_hit_t *rht)*/
/*{*/
	/*return 0;*/
/*}*/

int merge_seg_core2(eg_hit_t *rht, size_t s, size_t e, int gs, sdict_t *sn)
{
    /*
     qry -----------------
	 ref
	 |   \  \
	 |   \ \
	 |   
	 |
	 |
	 */ 
	size_t i, j, k;
 	//score: ml backtrace: bl 
	for ( i = s; i < e; ++i) {
		rht[i].tail = 1;
		rht[i].bl = -1;
		uint32_t qml = rht[i].ml;
		rht[i].del = 0;
		/*fprintf(stdout, "(%u, %u)\n", i, __LINE__ );*/
		/*print_hit(rht+i,tn);*/
		for ( j = i - 1; j + 1 > s; --j) {
			if (rht[i].rev == rht[j].rev) {
				// j u < i v
				if (rht[i].rev) {
					/*fprintf(stdout, "(%u, %u, %u)\n", j, i, __LINE__ );*/
					/*print_hit(rht+j,tn);*/
					/*print_hit(rht+i, tn);*/
					if (rht[j].qns < rht[i].qns && rht[j].qe < rht[i].qe && rht[j].te > rht[i].te && rht[j].tns > rht[i].tns) {
					/*if (rht[j].qe <= (uint32_t)rht[i].qns + 2000 && (uint32_t)rht[j].tns + 2000 >= rht[i].te ) {*/
						int qgs = (uint32_t)rht[i].qns - rht[j].qe;	
						/*uint32_t qgs = (uint32_t)rht[i].qns > rht[j].qe ? (uint32_t)rht[i].qns - rht[j].qe : 0;	*/
					/*fprintf(stdout, "pass1 (%u, %u, %u)\n", j, i, __LINE__ );*/
						int tgs = (uint32_t)rht[j].tns - rht[i].te;//any of them is almost nested
						if (qgs < 0) {
							int qlen = rht[j].qe - (uint32_t)rht[j].qns;								
							int qlen1 = rht[i].qe - (uint32_t)rht[i].qns;
							if ((float)(-qgs)/qlen > .95 || (float)(-qgs)/qlen1 > .95) continue; 	
						} 
					  	if (tgs < 0) {
							int tlen = rht[j].te - (uint32_t)rht[j].tns;								
							int tlen1 = rht[i].te - (uint32_t)rht[i].tns;
							if ((float)(-qgs)/tlen > .95 || (float)(-tgs)/tlen1 > .95) continue; 	
						}	
						if (tgs > gs || qgs > gs) continue;
					/*fprintf(stdout, "pass (%u, %u, %u)\n", j, i, __LINE__ );*/
						//otherwise
					} else 
						continue;
				} else {
					/*fprintf(stderr, "%u (%u, %u, %u)\n", is_sr, j, i, __LINE__ );*/
					/*print_hit(rht+j,sn);*/
					/*print_hit(rht+i, sn);*/
					if (rht[j].qns < rht[i].qns && rht[j].qe < rht[i].qe && rht[j].te < rht[i].te && rht[j].tns < rht[i].tns) {
					/*if (rht[j].qe <= (uint32_t)rht[i].qns + 2000 && (uint32_t)rht[i].tns + 2000 >= rht[j].te) {*/
						int qgs =  (uint32_t)rht[i].qns - rht[j].qe;
						int tgs =  (uint32_t) rht[i].tns - rht[j].te;
						if (qgs < 0) {
							int qlen = rht[j].qe - (uint32_t)rht[j].qns;								
							int qlen1 = rht[i].qe - (uint32_t)rht[i].qns;
							if ((float)(-qgs)/qlen > .95 || (float)(-qgs)/qlen1 > .95) continue; 	
						} 
					  	if (tgs < 0) {
							int tlen = rht[j].te - (uint32_t)rht[j].tns;								
							int tlen1 = rht[i].te - (uint32_t)rht[i].tns;
							if ((float)(-qgs)/tlen > .95 || (float)(-tgs)/tlen1 > .95) continue; 	
						}	
						if (qgs > gs || tgs > gs) continue;
					} else 
						continue;	
					/*fprintf(stderr, "pass (%u, %u, %u)\n", j, is_sr, __LINE__ );*/
				}		
				rht[j].tail = 0;
				rht[j].del = 1;
				if (rht[i].ml < rht[j].ml + qml) {
					rht[i].ml = rht[j].ml + qml;
					rht[i].bl = j;
				}	
			}		
		}
	}
	//update coordinate
	uint32_t max = 0;	
	uint32_t max_idx = -1;
	uint32_t sent = 0X3FFFFFFF;	
	for (i = s; i < e; ++i) {
		if (rht[i].tail) {
			//backtrace to find start 
			if (rht[i].ml > max) {
				max = rht[i].ml;
				max_idx = i;
			} 
			for (k = i; rht[k].bl != sent; k = rht[k].bl); 
			rht[i].qns = rht[k].qns;
			if (rht[i].rev) rht[i].te = rht[k].te;
			else rht[i].tns = rht[k].tns;	
		}
	}	
	return max_idx;
}


int merge_seg_core(eg_hit_t *rht, size_t s, size_t e, uint32_t max_gs, sdict_t *sn)
{
	size_t i, j;
	uint32_t rtn;
	/*print_hits(rht, s,e, sn, "merge_core");*/
	for ( i = j = s; i <= e; ++i) {
		if (i == e || (rht[j].tns >> 32) != (rht[i].tns >> 32)) {
			/*fprintf(stdout, "index %d\t%d\n", j, i);*/

			/*print_hits(rht, j,i, tn);*/
			rtn = merge_seg_core2(rht, j, i, max_gs, sn);
			/*if (~rtn) print_hit(rht+rtn, tn);*/
			j = i;	
		} 
	}
	return 0;
}

int get_bm_mm_core2(eg_hit_t *rht, size_t s, size_t e, sdict_t *sn, int *bmf, int *mmf)
{
	uint32_t bml = 0;
	uint32_t mml = 0;
	uint32_t ql = rht[s].ql;
	
	size_t i;
	uint32_t qs = 1, qe = 0;
	for ( i = s; i <= e; ++i) {
		if (i == e) {
			bml += (qe + 1 - qs);
			break;	
		}
		uint32_t cur_qs = (uint32_t)rht[i].qns;
		mml += (rht[i].qe + 1 - cur_qs);
		if ( cur_qs > qe) {
			bml += (qe + 1 - qs);
			qs = cur_qs;
			qe = rht[i].qe;
		} else if (rht[i].qe > qe){
			qe = rht[i].qe;	
		}
	}
	*bmf = bml * 100 / ql;
	*mmf = mml * 100 / ql;
	fprintf(stderr, "bmfore: %s\t%s\t%f\t%f\n", sn->seq[rht[s].qns >> 32].name, sn->seq[rht[s].tns >>32].name, (float) bml * 100/ql, (float) mml * 100/ql);
	return 0;
}


int get_bm_mm_core(eg_hit_t *rht, size_t s, size_t e, sdict_t *sn, int *sid, int *bm, int *mm)
{
	size_t i, j;
	uint32_t rtn;
	/*print_hits(rht, s,e, sn, "merge_core");*/
	int tmp_sid[2] = {-1, -1}, tmp_bm[2] = {-1, -1}, tmp_mm[2] = {-1, -1};
	int _bm, _mm;
	for ( i = j = s; i <= e; ++i) {
		if (i == e || (rht[j].tns >> 32) != (rht[i].tns >> 32)) {
			/*fprintf(stderr, "GETBM: %s\t%s\tis %s deleted\n", sn->seq[rht[j].qns >> 32].name, sn->seq[rht[j].tns>>32].name, sn->seq[rht[j].tns >> 32].del && sn->seq[rht[j].tns >> 32].del2 ? " ":" not");*/
			if (!sn->seq[rht[j].tns >> 32].del && !sn->seq[rht[j].tns >> 32].del2) {
				rtn = get_bm_mm_core2(rht, j, i, sn, &_bm, &_mm);
				if (_bm > tmp_bm[0])  tmp_bm[1] = tmp_bm[0], tmp_mm[1] = tmp_mm[0], tmp_sid[1] = tmp_sid[0], tmp_bm[0] = _bm, tmp_mm[0] = _mm,tmp_sid[0]= (rht[j].tns >> 32);
				else if (_bm > tmp_bm[1]) tmp_bm[1] = _bm, tmp_mm[1] = _mm, tmp_sid[1]= (rht[j].tns >> 32);	
			} 
			j = i;	
		} 
	}
	*bm = tmp_bm[0];
	/*fprintf(stderr,"bm: %d\t%d\t%d\n", tmp_mm[0], tmp_mm[1], tmp_sid[0]);*/
	*mm = tmp_mm[0] > tmp_mm[1] ? tmp_mm[0] : tmp_mm[1];
	*sid = tmp_sid[0];
	fprintf(stderr,"bm: %s\t%s\t%d\t%d\n", sn->seq[rht[s].qns >> 32].name, tmp_sid[0] == -1 ? "null" : sn->seq[tmp_sid[0]].name, *bm, *mm);
	return 0;
}
int flt_by_bm_mm(eg_hit_t *rht, uint64_t *idx, size_t n_idx, sdict_t *sn, dup_v *dups, int min_bm, int min_mm)
{
	size_t j;
	/*#pragma parallel for number_threads(4)*/
	int bm, mm, sid;
	uint32_t *srt_idx  = malloc(sizeof(uint32_t) * n_idx);
	for ( j = 0; j < n_idx; ++j ) srt_idx[j] = j;
	insertionSort(srt_idx, sn);
	while (1) {
		for (j = 0; j < n_idx; ++j) {
			//use type to indicate the candidates 
			if (!sn->seq[j].del && !sn->seq[j].del2 && sn->seq[j].type == 7 && sn->seq[j].aux != PRIMARY)  {
				get_bm_mm_core(rht , idx[j] >> 32, (idx[j] >> 32) + (uint32_t)idx[j], sn, &sid, &bm, &mm); 
				if (bm >= min_bm) {
					if (mm >= min_mm) {
						/*fprintf(stderr, "bm: %s\t%u\t%u\t%s\n", sn->seq[j].name, bm, mm, "REPEAT");*/
						/*dup_t k = {(uint32_t)j, REPEAT, 0, 1, sn->seq[j].len};*/
						/*kv_push(dup_t, *dups, k);*/
						sn->seq[j].type = REPEAT;
						sn->seq[j].best_hit = sid;
					} else {
						/*fprintf(stderr, "bm: %s\t%u\t%u\t%s\n", sn->seq[j].name, bm, mm, "HAPLOTIGS");*/
						/*dup_t k = {(uint32_t)j, HAPLOTIG, 0, 1, sn->seq[j].len};*/
						/*kv_push(dup_t, *dups, k);*/
						sn->seq[j].type = HAPLOTIG;
						sn->seq[j].best_hit = sid;
					}
				} else {
					sn->seq[j].type = PRIMARY;
					sn->seq[j].best_hit = sid;
				}	
			}
		}
		int check = 0;
		for (j = 0; j < n_idx; ++j) {
			uint32_t seq_idx = srt_idx[j];
			/*fprintf(stderr, "seq_idx: %d\n", seq_idx);*/
			if (!sn->seq[seq_idx].del && sn->seq[seq_idx].aux != PRIMARY && !sn->seq[seq_idx].del2) {
				uint32_t cur_tp = sn->seq[seq_idx].type;
				uint32_t best_hit_tp = sn->seq[seq_idx].best_hit  == 0xFFFFFFFF ? 7 : sn->seq[sn->seq[seq_idx].best_hit].type;
				//if both types are haplotigs or repeat 
				if (cur_tp == HAPLOTIG || cur_tp == REPEAT)  {
					if (best_hit_tp == HAPLOTIG || best_hit_tp == REPEAT) {
						fprintf(stderr, "RESET: %s\t%s\n", sn->seq[seq_idx].name, sn->seq[sn->seq[seq_idx].best_hit].name);
						sn->seq[sn->seq[seq_idx].best_hit].type = 7;
						sn->seq[sn->seq[seq_idx].best_hit].del2 = 0;
					}
				   sn->seq[seq_idx].del2 = 1;	
				   check = 1;
				}  
			}
		}
		if (!check) break;
	}
	fprintf(stderr, "[M::%s] check overpuring\n", __func__);
	//check overpurging 
	for (j = 0; j < n_idx; ++j) {
		uint32_t seq_idx = srt_idx[j];
		if (!sn->seq[seq_idx].del && sn->seq[seq_idx].aux != PRIMARY) {
			if (sn->seq[seq_idx].del2 && sn->seq[seq_idx].best_hit != 0xFFFFFFFF && sn->seq[sn->seq[seq_idx].best_hit].del2)  {

				get_bm_mm_core(rht , idx[seq_idx] >> 32, (idx[seq_idx] >> 32) + (uint32_t)idx[seq_idx], sn, &sid, &bm, &mm); 
				if (bm >= min_bm) {
					if (mm >= min_mm) {
						/*fprintf(stderr, "bm: %s\t%u\t%u\t%s\n", sn->seq[j].name, bm, mm, "REPEAT");*/
						/*dup_t k = {(uint32_t)j, REPEAT, 0, 1, sn->seq[j].len};*/
						/*kv_push(dup_t, *dups, k);*/
						sn->seq[seq_idx].type = REPEAT;
						sn->seq[seq_idx].best_hit = sid;
					} else {
						/*fprintf(stderr, "bm: %s\t%u\t%u\t%s\n", sn->seq[j].name, bm, mm, "HAPLOTIGS");*/
						/*dup_t k = {(uint32_t)j, HAPLOTIG, 0, 1, sn->seq[j].len};*/
						/*kv_push(dup_t, *dups, k);*/
						sn->seq[seq_idx].type = HAPLOTIG;
						sn->seq[seq_idx].best_hit = sid;
					}
				} else {
					sn->seq[seq_idx].type = PRIMARY;
					sn->seq[seq_idx].best_hit = sid;
					sn->seq[seq_idx].del2 = 0;
				}	
			}
		}
	}
	
	free(srt_idx);
	//write dups 
	for (j = 0; j < n_idx; ++j) {
		if (sn->seq[j].del2)	{
			fprintf(stderr, "TYPE: %s\t%s\n", sn->seq[j].name, dup_type_s[sn->seq[j].type]);
			dup_t k = (dup_t){(uint32_t)j, sn->seq[j].type, 0, 1, sn->seq[j].len, sn->seq[j].best_hit};	
			kv_push(dup_t, *dups, k);
			sn->seq[j].del = 1;
		}
	}	
	return 0;
}
/**
 * @func    merge_segs
 * @brief   merge segment if they are aligned closed to each other 
 * @alg     construct alignment graph, and find paths.  
 */
int merge_segs(eg_hit_t *rht, uint64_t *idx, size_t n_idx, uint32_t max_gs, sdict_t *sn)
{
	size_t j;
	/*#pragma parallel for number_threads(4)*/
	for (j = 0; j < n_idx; ++j) {
		merge_seg_core(rht , idx[j] >> 32, (idx[j] >> 32) + (uint32_t)idx[j], max_gs, sn); 
	}
	return 0;
}
/**
 * @func   flt_hits
 * @brief  remove nested alignments 
 * @comm   not as efficient as I want but have no idea to improve it
 */
#define MAX_LEFT 5000 //maybe use fraction
int flt_hits(eg_hit_t *rht, size_t n_rht)
{
	size_t i = 0, j = 0;
	for (i = 0; i < n_rht; ++i) {
		if (!rht[i].del) 
			for (j = i + 1; j < n_rht; ++j) {
				if (rht[j].del) continue;
				uint32_t qs = (uint32_t)rht[j].qns;  
				if ((rht[i].qns >> 32) != (rht[j].qns >> 32) || qs > rht[i].qe) break;
				if (rht[j].qe <= rht[i].qe && (rht[j].tns >> 32) == (rht[i].tns >> 32) && (uint32_t)rht[j].tns >= (uint32_t)rht[i].tns && rht[j].te <= rht[i].te)  
					rht[j].del = 1;			
			} 
	}
	return 0;
}

int flt_hits2(eg_hit_t *rht, size_t n_rht)
{
	size_t i = 0, j = 0;
	for (i = 0; i < n_rht; ++i) {
		if (!rht[i].del) 
			for (j = i + 1; j < n_rht; ++j) {
				uint32_t ts = (uint32_t)rht[j].tns;  
				if (ts > rht[i].te) break;
				if (ts >= (uint32_t)rht[i].tns && rht[j].te <= rht[i].te)  
					rht[j].del = 1;			
			} 
	}
	return 0;
}
/** 
 * @func	is_cont
 * @brief	if aln q is contained in aln p 
 * @notice	not helpful don't use
 */


int is_cont(eg_hit_t *p, eg_hit_t *q) 
{
	uint32_t max_ts = max((uint32_t) p->tns, (uint32_t)q->tns);
	uint32_t min_te = min(p->te, q->te);
	uint32_t p_qs, p_qe, q_qs, q_qe;
	int32_t p_ext5, p_ext3, q_ext5, q_ext3;
	p_qs = (uint32_t) p->qns + max_ts - (uint32_t) p->tns;
	p_qe = p->qe + min_te - p->te;
	q_qs = (uint32_t) q->qns + max_ts - (uint32_t) q->tns;
	q_qe = q->qe + min_te - q->te;
	if (p->rev) 
		p_ext5 = p->ql - p->qe, p_ext3 = p_qs;
	else
		p_ext3 = p->ql - p->qe, p_ext5 = p_qs;
	if (q->rev) 
		q_ext5 = q->ql - q->qe, q_ext3 = q_qs;
	else
		q_ext3 = q->ql - q->qe, q_ext5 = q_qs;
	return !(p_ext5-q_ext5) || !(p_ext3-q_ext3) || ((p_ext5 - q_ext5) ^ (p_ext3 - q_ext3)) >= 0; 
}


/**
 * @func   set_cont
 * @brief  mark contained queries
 */

int set_cont(eg_hit_t *rht, uint64_t *idx, size_t n_idx)
{
	size_t i = 0, last;
		
	for (i = 0; i < n_idx; ++i) {
		for (last = idx[i]>>32, i = last + 1; i < (idx[i]>>32) + (uint32_t)idx[i]; ++i) {
			/*if (is_cont(rht + i, rht + i-1))*/
			if (rht[i].te < rht[i-1].te)
				rht[i].con = 1;
			else {
				rht[i].con = 0;	
				last = i;
			}
		}
	}
	return 0;	
}

/**
 * @func   cleanup_hits 
 * @brief  clean up hits after 
 */

size_t cleanup_hits(eg_hit_t *rht, size_t n_rht, sdict_t *sn)
{
	size_t i, j;
	for (i = 0, j = 0; i < n_rht; ++i) 
		if (!rht[i].del && !sn->seq[rht[i].qns>>32].del && !sn->seq[rht[i].tns >> 32].del) 
			rht[j++] = rht[i];
	return j;
}

int double_hits(eg_hits_t *rhts)
{
		eg_hit_t *p;
		eg_hit_v h = {rhts->n, rhts->m, rhts->rht};
		size_t i = 0;
		size_t n_rht = rhts->n;
		eg_hit_t *rht = rhts->rht;
		for ( i = 0; i < n_rht; ++i) {
			kv_pushp(eg_hit_t, h, &p);
			p->qns = rht[i].tns;	
			p->ql = rht[i].tl; p->qe = rht[i].te;
			p->qtg = rht[i].ttg;
			p->tns = rht[i].qns;
			p->tl = rht[i].ql; p->te = rht[i].qe;
			p->ttg = rht[i].qtg;
			p->rev = rht[i].rev; p->ml = rht[i].ml; p->bl = rht[i].bl = i; //this is identifier now	
			p->con = rht[i].con; p->del = rht[i].del;	p->tcov = rht[i].qcov; p->qcov = rht[i].tcov;
			p->qtp = rht[i].ttp, p->ttp = rht[i].qtp;
		}	
		rhts->n *= 2;
		return 0;

}

/**
 * @func     update_cords
 * @brief    update coordinates for hits 
 * @notice   keep ml bl unchanged
 */
int update_cords(eg_hit_t *rht, size_t n_rht, uint32_t ext)
{
	size_t i;
	for (i = 0; i < n_rht; ++i) {
		uint32_t qs = (uint32_t) rht[i].qns;
		uint32_t ts = (uint32_t) rht[i].tns;
		
		if (rht[i].ql - rht[i].qe < ext && qs < ext) {
			rht[i].qtg = 3;
			rht[i].qns = (rht[i].qns >> 32 << 32 ) | 1;
			rht[i].qe = rht[i].ql;
		} else if (rht[i].ql - rht[i].qe < ext) {
			rht[i].qtg = 2;
			rht[i].qe = rht[i].ql;	
		} else if (qs < ext) {
			rht[i].qtg = 0;
			rht[i].qns = (rht[i].qns >> 32 << 32 ) | 1;
		} else 
			rht[i].qtg = 1;	
			
		if (rht[i].tl - rht[i].te < ext && ts < ext) {
			rht[i].ttg = 3;
			rht[i].tns = (rht[i].tns >> 32 << 32 ) | 1;
			rht[i].te = rht[i].tl;
		} else if (rht[i].tl - rht[i].te < ext) {
			rht[i].ttg = 2;
			rht[i].te = rht[i].tl;	
		} else if (ts < ext) {
			rht[i].ttg = 0;
			rht[i].tns = (rht[i].tns >> 32 << 32 ) | 1;
		} else 
			rht[i].ttg = 1;	
	}
	return 0;
}

/**
 * @func    hit_read
 * @brief   read hits from paf file *
 */
int hit_read(char *paf_fn, eg_hits_t *rhts, sdict_t *sn, uint32_t min_match, int is_s2s)
{
	
	paf_file_t *fp;
	paf_rec_t r;
	eg_hit_v h = {0, 0, 0};
	fp = paf_open(paf_fn);
	if (!fp) {
		fprintf(stderr, "[E::%s] can not open PAF file %s\n",__func__, paf_fn);
		exit(1);	
	}	
	while (paf_read(fp, &r) >= 0) {
		/*if (r.ml >= min_match && (r.ql < r.tl || (r.ql == r.tl && strcmp(r.qn, r.tn) >= 0))) {*/
		if ((r.ql < r.tl || (r.ql == r.tl && strcmp(r.qn, r.tn) > 0))) {
			eg_hit_t *p;
			kv_pushp(eg_hit_t, h, &p);
			p->qns = (uint64_t)sd_put(sn, r.qn, r.ql, 1)<<32 | (r.qs + 1);// one base	
			p->ql = r.ql; p->qe = r.qe; //one base
			p->tns = (uint64_t)sd_put(sn, r.tn, r.tl, 0) << 32 | (r.ts + 1); //one base
			p->tl = r.tl; p->te = r.te; //one base
			p->rev = r.rev; p->ml = r.ml; p->bl = r.bl; 	
			p->con = 0; p->del = 0;	p->tcov = p->qcov = 0;
		}
	}
	paf_close(fp);
	radix_sort_hit(h.a, h.a + h.n);
	rhts->rht = h.a;
	rhts->n = h.n;
	rhts->m = h.m;
	return 0; 
}	



eg_hits_t *eg_init()
{
	return (eg_hits_t *)calloc(1, sizeof(eg_hits_t));
}


int eg_destroy(eg_hits_t *r)
{
	if (r) {
		if (r->idx) free(r->idx);
		if (r->rht) free(r->rht);
		free(r);
	}
	return 0;

}

/* 
int purge_contigs_core(eg_hit_t *rht, size_t s, size_t e, sdict_t *sn, dup_v *dups, uint32_t ctg_gap)
{
	size_t i;
	uint32_t st = (uint32_t)rht[s].qns, ed = rht[s].qe;
	for ( i = s + 1; i < e; ++i ) {
		if ((uint32_t)rht[i].qns > ed  + ctg_gap) //replace with params;
	   		break;
		else if (rht[i].qe > ed)
			ed = rht[i].qe;
	}
	if (sn->seq[rht[s].qns >> 32].len  + st == ed + 1) {
		sn->seq[rht[s].qns >> 32].del = 1;
		dup_t k = (dup_t) {(uint32_t)(rht[s].qns >>32), PRIMARY, 0, st, ed, };	
		kv_push(dup_t, *dups, k);
	}
	return 0;	
}
 */
/*  
int purge_contigs(eg_hits_t *rhts, sdict_t *sn, dup_v *dups, uint32_t ctg_gap)
{
	eg_hit_t *rht = rhts->rht;
	uint64_t *idx = rhts->idx;
	size_t n_idx = sn->n_seq;
	size_t j;
	for (j = 0; j < n_idx; ++j) {
		purge_contigs_core(rht , idx[j] >> 32, (idx[j] >> 32) + (uint32_t)idx[j], sn, dups, ctg_gap); 
	}
	return 0;
}
*/
int flt_hits_by_ml(eg_hit_t *rht, size_t n, uint32_t min_len)
{
	size_t i;
	for ( i = 0; i < n; ++i ) 
		if (rht[i].ml < min_len) 
			rht[i].del = 1;
	return 0;
}

int get_cuts(char *fn, uint32_t *cutoffs)
{
	if (!fn) return 1;
	FILE *fp = fopen(fn, "r");
	if (!fp) return 1;
	int count = fscanf(fp, "%d\t%d\t%d\t%d\t%d\t%d\n", cutoffs, cutoffs + 1, cutoffs + 2, cutoffs + 3, cutoffs + 4, cutoffs + 5);	
	fclose(fp);
	if (count == 5) cutoffs[5] = 160;
	if (count < 5) return 1;
	else
		return 0;
}

uint32_t bin_srch(cov_ary_t *ca, uint32_t ca_n, uint32_t p)
{
	if (ca_n == 0) return -1;
	uint32_t l = 0, h = ca_n -1;
	while (l < h) {
		uint32_t m = (l + h) >> 1;
		/*fprintf(stderr, "%u %u %u %p %u %u\n", l, h, m, &ca->intv[m], ca->intv[m].s, ca->intv[m].e);*/
		if (ca->intv[m].s > p)
			h = m - 1;
		else if (ca->intv[m].e < p) 
			l = m + 1;
		else 
			l = h = m;	
	}		
	return l;
}
typedef struct {
	char *ctgn;
	uint32_t s, e, nl;
}name_t;

int parse_name(char *s, int l, name_t *nt)
{
	char *q, *r;
	int i, t;
	for (i = t = 0, q = s; i <= l; ++i) {
		if (i < l && s[i] != ':' && s[i] != '-') continue;
		s[i] = 0;
		if (t == 0) nt->ctgn = q, nt->nl = i;
		else if (t == 1) nt->s = strtol(q, &r, 10), s[i] = '-';
		else if (t == 2) nt->e = strtol(q, &r, 10);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 2) return -1;
	return 0;

}
int classify_seq(cov_ary_t *ca, sdict_t *sn, sdict_t *osn, uint32_t* cutoffs, float min_frac)
{
	size_t i;
	name_t nt;
	uint32_t low_cut = cutoffs[0];
	uint32_t dip_cut = cutoffs[3];
	uint32_t high_cut = cutoffs[5];

	for ( i = 0; i < sn->n_seq; ++i ) {
		sd_seq_t * cur_seq = &sn->seq[i];
		uint32_t hap_n, dip_n, high_n, low_n; 
		hap_n = dip_n = high_n = low_n = 0;

		char *name = cur_seq->name;
		parse_name(name, strlen(name), &nt);	
		
		uint32_t s1s = nt.s;
		uint32_t s1e = nt.e;
		uint32_t sid = sd_get(osn, nt.ctgn);
		name[nt.nl] = ':';
		
		uint32_t bs, be, s;
		cur_seq->aux = JUNK;
		if (~sid && ca[sid].n) {
			if (s1s == 1 && s1e == cur_seq->len) 
				bs = 0, be = ca[sid].n - 1;
			else {
				bs = bin_srch(&ca[sid], ca[sid].n, s1s);	
				be = bin_srch(&ca[sid], ca[sid].n, s1e);	
			}
			
			s = s1s;
			uint32_t _cov;
			for ( ; bs < be; ++bs ) {
				uint32_t e = ca[sid].intv[bs].e;
				_cov = ca[sid].intv[bs].coverage;
				if (_cov <= low_cut) low_n += ( e - s + 1);
			   	else if (_cov <= dip_cut) hap_n += (e - s + 1);
				else if (_cov > high_cut) high_n += (e - s + 1);
				else dip_n += (e - s  + 1);	
				s = ca[sid].intv[bs + 1].s;			
			}
			_cov = ca[sid].intv[be].coverage;
			if (_cov <= low_cut) low_n += ( s1e - s + 1);
			else if (_cov <= dip_cut) hap_n += (s1e - s + 1);
			else if (_cov > high_cut) high_n += (s1e - s + 1);
			else dip_n += (s1e - s  + 1);	
			if ((float) low_n / (s1e - s1s + 1) > min_frac) cur_seq->aux = JUNK;
			else if ((float)dip_n / (s1e - s1s + 1) > min_frac) cur_seq->aux = PRIMARY;
			else if ((float)high_n / (s1e - s1s + 1) > min_frac)	cur_seq->aux = REPEAT;
			else cur_seq->aux = HAPLOTIG;
		}
		fprintf(stderr, "%s\t%u\t%u\t%u\t%u\n", cur_seq->name, low_n, hap_n, dip_n, high_n);
	}
	return 0;
}

int cal_cov_4reg(eg_hit_t *rht, size_t n_rht, cov_ary_t *ca, sdict_t *sn, sdict_t *osn, uint32_t* cutoffs, float min_frac)
{
	size_t i;
	name_t nt;
	uint32_t low_cut = cutoffs[0];
	uint32_t dip_cut = cutoffs[3];
	uint32_t high_cut = cutoffs[5];
	for ( i = 0; i < n_rht; ++i ) {
		print_hit(rht+i, sn);
		uint32_t hap_n, dip_n, high_n, low_n; 
		hap_n = dip_n = high_n = low_n = 0;
		char *name = sn->seq[rht[i].qns >> 32].name;
		parse_name(name, strlen(name), &nt);	
		uint32_t s1s = nt.s + (uint32_t)rht[i].qns - 1;
		uint32_t s1e = nt.s + rht[i].qe - 1;
		uint32_t sid = sd_get(osn, nt.ctgn);
		name[nt.nl] = ':';
		uint32_t bs, be, s;
		uint64_t coverage;
		rht[i].qtp = JUNK;
		if (~sid) {
			bs = bin_srch(&ca[sid], ca[sid].n, s1s);	
			be = bin_srch(&ca[sid], ca[sid].n, s1e);	

			s = s1s;
			coverage = 0;
			uint32_t _cov;
			for ( ; bs < be; ++bs ) {
				uint32_t e = ca[sid].intv[bs].e;
				_cov = ca[sid].intv[bs].coverage;
				coverage += _cov * (e -s + 1);
				if (_cov <= low_cut) low_n += ( e - s + 1);
			   	else if (_cov <= dip_cut) hap_n += (e - s + 1);
				else if (_cov > high_cut) high_n += (e - s + 1);
				else dip_n += (e - s  + 1);	
				s = ca[sid].intv[bs + 1].s;			
			}
			_cov = ca[sid].intv[be].coverage;
			coverage += (s1e - s + 1) * _cov;
			rht[i].qcov = coverage / (s1e - s1s + 1);
			if (_cov <= low_cut) low_n += ( s1e - s + 1);
			else if (_cov <= dip_cut) hap_n += (s1e - s + 1);
			else if (_cov > high_cut) high_n += (s1e - s + 1);
			else dip_n += (s1e - s  + 1);	
			if ((float) low_n / (s1e - s1s + 1) > min_frac) rht[i].qtp = JUNK;
			else if ((float)hap_n / (s1e - s1s + 1) > min_frac) rht[i].qtp = HAPLOTIG;
			else if ((float)dip_n / (s1e - s1s + 1) > min_frac)	rht[i].qtp = PRIMARY;
			else rht[i].qtp = REPEAT;
		}

		hap_n = dip_n = high_n = low_n = 0;
		name = sn->seq[rht[i].tns >> 32].name;
		parse_name(name, strlen(name), &nt);	
		s1s = nt.s + (uint32_t)rht[i].tns - 1;
		s1e = nt.s + rht[i].te - 1;
		sid = sd_get(osn, nt.ctgn);
		name[nt.nl] = ':';
		rht[i].ttp = JUNK;
		if (~sid) {
			bs = bin_srch(&ca[sid], ca[sid].n, s1s);	
			be = bin_srch(&ca[sid], ca[sid].n, s1e);	
			s = s1s;
			coverage = 0;
			
			uint32_t _cov;
			for ( ; bs < be; ++bs ) {
				uint32_t e = ca[sid].intv[bs].e;
				_cov = ca[sid].intv[bs].coverage;
				coverage += _cov * (e -s + 1);
				if (_cov <= low_cut) low_n += ( e - s + 1);
			   	else if (_cov <= dip_cut) hap_n += (e - s + 1);
				else if (_cov > high_cut) high_n += (e - s + 1);
				else dip_n += (e - s  + 1);	
				s = ca[sid].intv[bs + 1].s;			
			}
			_cov = ca[sid].intv[be].coverage;
			coverage += (s1e - s + 1) * _cov;
			rht[i].tcov = coverage / (s1e - s1s + 1);
			if (_cov <= low_cut) low_n += ( s1e - s + 1);
			else if (_cov <= dip_cut) hap_n += (s1e - s + 1);
			else if (_cov > high_cut) high_n += (s1e - s + 1);
			else dip_n += (s1e - s  + 1);	
			if ((float) low_n / (s1e - s1s + 1) > min_frac) rht[i].ttp = JUNK;
			else if ((float)hap_n / (s1e - s1s + 1) > min_frac) rht[i].ttp = HAPLOTIG;
			else if ((float)dip_n / (s1e - s1s + 1) > min_frac)	rht[i].ttp = PRIMARY;
			else rht[i].ttp = REPEAT;
		}
	}	
	return 0;
}

//

int flt_hits4(eg_hit_t *rht, size_t n, sdict_t *sn)
{
	size_t i, j;
	uint8_t *del = calloc(n>>1, sizeof(uint8_t));
	for ( i = 0, j = 0; i <=n ; ++i) {
		if (i == n || (rht[i].qns >> 32) != (rht[j].qns >> 32)) {
			//find out the largest 
			size_t z;
			/*fprintf(stderr, "%s %u %u\n", sn->seq[rht[j].qns>>32].name, maxsloc, mineloc);*/
			for ( z = j; z < i; ++z) {
				size_t t;
				for ( t = z + 1; t < i; ++t) {
					//check if z is removed 
					if (rht[t].qe > rht[z].qe) break; //we only look for nested queries 
					if (del[rht[t].bl]) continue;	
					//check if t is nested in z ? 
					if (rht[z].qe >= rht[t].qe && (float)(rht[t].qe - (uint32_t)rht[t].qns)/(rht[z].qe - (uint32_t)rht[z].qns) <.85)  //don't remove repetitive elements
						del[rht[t].bl] = 1;
				}
			}
			j = i;
		}
	}
	for ( i = 0; i < n; ++i) 
		if (del[rht[i].bl]) {
			rht[i].del = 1;
		}
	free(del);	
	return 0;
}

//compare matched length remove smaller one
/*  
int flt_hits4(eg_hit_t *rht, size_t n, sdict_t *sn) 
{
	size_t i, j;
	uint8_t *del = calloc(n>>1, sizeof(uint8_t));
	uint32_t maxsloc = 1;
	uint32_t mineloc = rht[0].ql;
	for ( i = 0, j = 0; i <=n ; ++i) {
		if (i == n || (rht[i].qns >> 32) != (rht[j].qns >> 32)) {
			//find out the largest 
			size_t z;
			for ( z = j; z < i; ++z) {
				if (del[rht[z].bl]) continue;
				if (rht[z].qtg == 0) {
					if ((float)(rht[z].qe - (uint32_t)rht[z].qns + 1)/maxsloc < .85) 
						del[rht[z].bl] = 1;	
				} else 	if (rht[z].qtg == 1) {
					if (rht[z].qe <= maxsloc) del[rht[z].bl] = 1; 
					else if ((uint32_t)rht[z].qns >= mineloc) del[rht[z].bl] = 1;	
				} else {
					if ((float)(rht[z].qe - (uint32_t)rht[z].qns)/(rht[z].ql - mineloc) < .85) {
						del[rht[z].bl] = 1;	
					}
				}	
			}
			if (i == n) break;
			else {
				maxsloc = 1;
				mineloc = rht[i].ql;
				j = i--;
			}
		} else {
			if (rht[i].qtg == 0) {
				if (rht[i].qe > maxsloc) 
					maxsloc = rht[i].qe;
			} else if (rht[i].qtg == 2) {
				if ((uint32_t)rht[i].qns < mineloc) 
					mineloc = (uint32_t)rht[i].qns;	
			}
		}
	}
	for ( i = 0; i < n; ++i) 
		if (del[rht[i].bl]) {
			rht[i].del = 1;
		}
	free(del);	
	return 0;
}

*/
int flt_hits3(eg_hit_t *rht, size_t n, uint32_t *cutoffs, int usecuts, sdict_t *sn, dup_v *dups) 
{
	size_t i;
	uint32_t  dip_cut = cutoffs[3];
	/*uint32_t high_cut = cutoffs[5];	*/
	/*char *name;*/
	//purge contigs;
	for ( i = 0; i < n; ++i ) {
			if (!usecuts && rht[i].qtg == 3) {

				sn->seq[rht[i].qns >> 32].del = 1;
				/*name = sn->seq[rht[i].qns >> 32].name;*/
				/*parse_name(name, strlen(name), &nt);*/
				fprintf(stderr, "removed:\t%s\t%u\t%u\t%s\t%u\t%u\n", sn->seq[rht[i].qns>>32].name, rht[i].qcov, cutoffs[2], sn->seq[rht[i].tns>>32].name, rht[i].tcov, cutoffs[2]);
				dup_t k = (dup_t) {(uint32_t)(rht[i].qns >>32), HAPLOTIG, 0, 1, rht[i].qe, (uint32_t)(rht[i].tns >> 32)};	
				/*dup_t k = (dup_t) {sd_get(osn, nt.ctgn), nt.s, nt.e};	*/
				kv_push(dup_t, *dups, k);
				rht[i].del = 1;
			} else if (!usecuts && rht[i].ttg == 3) {
				sn->seq[rht[i].tns >> 32].del = 1;	
				fprintf(stderr, "removed:\t%s\t%u\t%u\t%s\t%u\t%u\n", sn->seq[rht[i].tns>>32].name, rht[i].tcov, cutoffs[2], sn->seq[rht[i].qns>>32].name, rht[i].qcov, cutoffs[2]);
				/*name = sn->seq[rht[i].qns >> 32].name;*/
				/*parse_name(name, strlen(name), &nt);*/
				dup_t k = (dup_t) {(uint32_t)(rht[i].tns >>32), HAPLOTIG, 0,  1, rht[i].te,(uint32_t)(rht[i].tns >> 32)};	
				kv_push(dup_t, *dups, k);
				rht[i].del = 1;
			} else if (usecuts && (rht[i].ttg == 3 || rht[i].qtg == 3))  //ususally not possible
				rht[i].del = 1;
	}
	for ( i = 0; i < n ; ++i ) {
		if (rht[i].qtg == 3 || rht[i].ttg == 3) continue; 
		if (usecuts && (rht[i].qcov > dip_cut || rht[i].tcov > dip_cut)) 
			rht[i].del = 1;
		else if (sn->seq[rht[i].qns >> 32].del || sn->seq[rht[i].tns>>32].del) 
			rht[i].del = 1;
	}	
	return 0;
}

int purge_dups2(eg_hit_t *rht, size_t n,sdict_t *sn, dup_v *dups) //second round to purge continous query 
{
	size_t i, j;
	uint32_t max_idx = 0;
	for ( i = 0; i < n; ++i) if (rht[i].bl > max_idx) max_idx = rht[i].bl;
	uint8_t *mask = calloc(max_idx + 1, sizeof(char)); //not real occurence if occures mulitple times set to 0; 
	for ( i =1, j = 0; i <= n; ++i ) {
		if (i == n || (rht[i].qns >> 32) != (rht[j].qns >> 32) || rht[i].qtg != rht[j].qtg) {
			if (i - j != 1) {
				if (rht[j].qtg == 1) {
					while ( j < i - 1) {
						float ovlp = rht[j].qe > (uint32_t)rht[j+1].qns? rht[j].qe - (uint32_t)rht[j+1].qns : 0;
						if (ovlp /(rht[j].qe - (uint32_t)rht[j].qns) >= .85 || ovlp/(rht[j+1].qe - (uint32_t)rht[j+1].qns) >= .85) 
							mask[rht[j].bl] = mask[rht[j+1].bl] = 1; 	
						++j;
					}
				} else 
					for (;j < i; ++j) mask[rht[j].bl] = 1; 
			} 
			j = i;
		} 
	}	
	for ( i = 0; i < n; ++i) {
		if ((rht[i].qns >> 32) < (rht[i].tns >> 32)) {
			/*fprintf(stderr, "sn: %s %d %d %s %d %d\n",sn->seq[rht[i].qns >> 32].name, occ_cnt[((rht[i].qns>>32)<<2) | rht[i].qtg ], rht[i].qtg, sn->seq[rht[i].tns>>32].name, rht[i].ttg, occ_cnt[((rht[i].tns >> 32) << 2) | rht[i].ttg]);*/
			if (!mask[rht[i].bl]) {
				if (sn->seq[rht[i].qns >> 32].len < sn->seq[rht[i].tns >>32].len) {
					dup_t k = (dup_t) {(uint32_t)(rht[i].qns >>32), OVLP, 0, (uint32_t)rht[i].qns, rht[i].qe, (uint32_t)(rht[i].tns >> 32)};	
					kv_push(dup_t, *dups, k);
				} else {
					dup_t k = (dup_t) {(uint32_t)(rht[i].tns >>32), OVLP, 0, (uint32_t)rht[i].tns, rht[i].te,(uint32_t)(rht[i].qns >> 32)};	
					kv_push(dup_t, *dups, k);
				}
			}	
		}
	}
	free(mask);
	return 0;
}
int purge_dups(eg_hit_t *rht, size_t n,sdict_t *sn, dup_v *dups) //second round to purge continous query 
{
	size_t i, j;

	uint8_t *occ_cnt = calloc(sn->n_seq * 4, sizeof(char)); //not real occurence if occures mulitple times set to 0; 
	for ( i =1, j = 0; i <= n; ++i ) {
		if (i == n || (rht[i].qns >> 32) != (rht[j].qns >> 32) || rht[i].qtg != rht[j].qtg) {
			if (i - j == 1) 
				occ_cnt[((rht[j].qns >> 32) << 2) | rht[j].qtg] = 1;
			else if (rht[j].qtg == 1) {
				while ( j < i - 1) {
					float ovlp = rht[j].qe > (uint32_t)rht[j+1].qns? rht[j].qe - (uint32_t)rht[j+1].qns : 0;
					if (ovlp /(rht[j].qe - (uint32_t)rht[j].qns) >= .85 || ovlp/(rht[j+1].qe - (uint32_t)rht[j+1].qns) >= .85)
						break;		
					++j;
				}
				if (j + 1 == i)	
					occ_cnt[(rht[j].qns >> 32) << 2 | rht[j].qtg] = 1;
				else 
					occ_cnt[(rht[j].qns >> 32) << 2 | rht[j].qtg] = 0;
			} else 
					occ_cnt[(rht[j].qns >> 32) << 2 | rht[j].qtg] = 0;
			j = i;
		} 
	}	
	
	for ( i = 0; i < n; ++i) {
		if ((rht[i].qns >> 32) < (rht[i].tns >> 32)) {
			/*fprintf(stderr, "sn: %s %d %d %s %d %d\n",sn->seq[rht[i].qns >> 32].name, occ_cnt[((rht[i].qns>>32)<<2) | rht[i].qtg ], rht[i].qtg, sn->seq[rht[i].tns>>32].name, rht[i].ttg, occ_cnt[((rht[i].tns >> 32) << 2) | rht[i].ttg]);*/
			if (occ_cnt[((rht[i].qns >> 32) << 2) | rht[i].qtg] == 1 && occ_cnt[((rht[i].tns >> 32) << 2) | rht[i].ttg] == 1) {
				if (sn->seq[rht[i].qns >> 32].len < sn->seq[rht[i].tns >>32].len) {
					dup_t k = (dup_t) {(uint32_t)(rht[i].qns >>32), HAPLOTIG, 0, (uint32_t)rht[i].qns, rht[i].qe, (uint32_t)(rht[i].tns >> 32)};	
					kv_push(dup_t, *dups, k);
				} else {
					dup_t k = (dup_t) {(uint32_t)(rht[i].tns >>32), HAPLOTIG, 0, (uint32_t)rht[i].tns, rht[i].te, (uint32_t) (rht[i].qns >> 32)};	
					kv_push(dup_t, *dups, k);
				}
			}	
		}
	}
	free(occ_cnt);
	return 0;
}


int cmp_dupt(const void *a, const void *b)
{
	dup_t *p = (dup_t *)a;
	dup_t *q = (dup_t *)b;
	if (p->sn > q->sn) return 1;
	else if (p->sn == q->sn) {
		if (p->s > q->s) 
			return 1;
		else if (p->s == q->s) return 0;
	   	else return -1;	
	} else return -1;
}


int update_dup_cords(dup_v *dups, sdict_t *sn, sdict_t *dup_n) 
{
	char *name; 
	name_t nt;
	dup_t *dp = dups->a;
	size_t n = dups->n;
	size_t i;
	for ( i = 0; i < n; ++i) {
		name = sn->seq[dp[i].sn].name;
		parse_name(name, strlen(name), &nt);
		uint32_t idx = sd_put(dup_n, nt.ctgn, 0, 1);
		name[nt.nl] = ':';
		dp[i].sn = idx;
		dp[i].s = nt.s + dp[i].s - 1;
		dp[i].e = nt.s + dp[i].e - 1;
		dp[i].del = 0;
		if (dp[i].bst_sn != 0XFFFFFFFF) {
			name = sn->seq[dp[i].bst_sn].name;
			parse_name(name, strlen(name), &nt);
			idx = sd_put(dup_n, nt.ctgn, 0, 1);
			dp[i].bst_sn = idx;
		}
	}
	return 0;	
}	

int merge_dups(dup_v *dups)
{
	dup_t *dp = dups->a;
	size_t n = dups->n;
	size_t i,j;
	if (!n) return 0;
	uint32_t s = dp[0].s, e= dp[0].e, tp = dp[0].tp, bst_sn = dp[0].bst_sn;
	for ( i = 1, j= 0; i <= n; ++i) {
		if (i == n || dp[i].sn != dp[j].sn || dp[i].s > e) {
			dp[j].s = s, dp[j].e = e, dp[j].del = 0, dp[j].tp = tp, dp[j].bst_sn = bst_sn;
			if (i != n) s=dp[i].s, e= dp[i].e, tp = dp[i].tp, bst_sn = dp[i].bst_sn;
		   	j = i;	
		} else {
			if (dp[i].e > e) 
				e = dp[i].e;
			if (dp[i].tp != tp) //I guess this wont happend, but what if report to me? best_sn != j.best_sn ? kind ignore this right now would be a problem? 
				tp = UNKNOWN;
			dp[i].del = 1;
		}
	}
	
	for ( i = 0, j = 0; i < n; ++i) {
		if (!dp[i].del) dp[j++] = dp[i];
	}
	dups->n = j;
	return 0;
}

int flt_by_tp(sdict_t *sn, dup_v *dups)
{
	size_t i, n = sn->n_seq;
	for ( i = 0; i < n; ++i) {
		if (sn->seq[i].aux == JUNK ) {
			fprintf(stderr, "tp:%s\n", sn->seq[i].name);
			dup_t k = (dup_t){(uint32_t)i, JUNK, 0, 1, sn->seq[i].len, 0xFFFFFFFF};
			kv_push(dup_t, *dups, k);
			sn->seq[i].del = 1;
		} else if (sn->seq[i].aux == REPEAT) {
			fprintf(stderr, "tp:%s\n", sn->seq[i].name);
			dup_t k = (dup_t){(uint32_t)i, REPEAT, 0, 1, sn->seq[i].len, 0xFFFFFFFF};
			kv_push(dup_t, *dups, k);
			sn->seq[i].del = 1;
		}
	}
	return 0;
}

int main(int argc, char *argv[])
{
	opt opts;
	if (parse_args(argc, argv, &opts))    return 1;		
	
#ifdef DEBUG
	fprintf(stderr, "[M::%s] finish parsing params\n", __func__);
#endif	
	/*fprintf(stderr,"[M::%s] parsing gfa...\n", __func__);*/
	//read gfa from gfa file
	/*gfa_t *gf = gfa_read(opts.gfa_fn);*/
	/*gfa_print(gf, stdout, 1);*/
	//read hits from  paf files
	sdict_t *sn = sd_init();
	eg_hits_t *rhts = eg_init();
	hit_read(opts.paf_fn, rhts, sn, opts.min_bl, opts.s2s);		
#ifdef DEBUG
	fprintf(stderr, "[M::%s] finish reading hits\n", __func__);
#endif	
	//apply a similar strategy as purge_haplotigs to remove haplotigs
	uint32_t cutoffs[6];
	dup_v dups = {0, 0, 0};
	size_t n_ind = sn->n_seq;
	cov_ary_t *ca = 0;
	sdict_t *osn = 0;
	if (!get_cuts(opts.cut_fn, cutoffs)) {
#ifdef DEBUG
	fprintf(stderr, "[M::%s] finish reading cutoffs\n", __func__);
#endif	
		osn = sd_init();
		ca = read_covs(opts.cov_fn, osn);
#ifdef DEBUG
	fprintf(stderr, "[M::%s] finish reading coverages\n", __func__);
#endif	

		classify_seq(ca, sn, osn, cutoffs, opts.min_frac);		
#ifdef DEBUG
	fprintf(stderr, "[M::%s] finish classifying sequences\n", __func__);
#endif	
		//remove dups and junks 
		flt_by_tp(sn, &dups);	
#ifdef DEBUG
	fprintf(stderr, "[M::%s] finish filtering by tags\n", __func__);
#endif	
		rhts->n = cleanup_hits(rhts->rht, rhts->n, sn);
#ifdef DEBUG
	fprintf(stderr, "[M::%s] finish cleaning\n", __func__);
#endif
	rhts->idx = hit_index(rhts->rht, rhts->n, n_ind, f_qn);
	mt_sort(rhts->rht, rhts->idx, n_ind, cmp_qtn);
#ifdef DEBUG
	fprintf(stderr, "[M::%s] finish sorting\n", __func__);
#endif	
#ifdef DEBUG
	fprintf(stderr, "[M::%s] finish indexing sequences\n", __func__);
#endif	
	flt_by_bm_mm(rhts->rht, rhts->idx, n_ind, sn, &dups, opts.min_bmf, opts.min_mmf);	
#ifdef DEBUG
	fprintf(stderr, "[M::%s] finish reassigning sequences by bm and mm\n", __func__);
#endif	
		rhts->n = cleanup_hits(rhts->rht, rhts->n, sn);
	}  
	//merge left hits
	flt_hits_by_ml(rhts->rht, rhts->n, opts.min_bl);
	rhts->n = cleanup_hits(rhts->rht, rhts->n, sn);
	
	size_t rnd = opts.sr ? 2 : 1;
	size_t i;
	int max_gs = opts.max_gs;
	for ( i = 0; i < rnd; ++i ) {
		if (rhts->idx) free(rhts->idx);
		rhts->idx = hit_index(rhts->rht, rhts->n, n_ind, f_qn);
		mt_sort(rhts->rht, rhts->idx, n_ind, cmp_qtn);
		merge_segs(rhts->rht, rhts->idx, n_ind, max_gs, sn);
		rhts->n = cleanup_hits(rhts->rht, rhts->n, sn);
		qsort(rhts->rht, rhts->n, sizeof(eg_hit_t), cmp_q);	
		flt_hits(rhts->rht, rhts->n);
		rhts->n = cleanup_hits(rhts->rht, rhts->n, sn);
		max_gs = opts.max_gs2rd;
	}
	
	flt_hits_by_ml(rhts->rht, rhts->n, opts.min_dup_bl); //replace by parameter
	rhts->n = cleanup_hits(rhts->rht, rhts->n, sn);
	print_hits(rhts->rht, 0, rhts->n, sn, "merged");
	update_cords(rhts->rht, rhts->n, opts.max_ext); // same as above
	
	if (ca) {
		cal_cov_4reg(rhts->rht, rhts->n, ca, sn, osn, cutoffs, opts.min_frac);		
		flt_hits3(rhts->rht, rhts->n,cutoffs, 1, sn, &dups);
		sd_destroy(osn);
	} else 
		flt_hits3(rhts->rht, rhts->n,cutoffs, 0, sn, &dups);

	
	rhts->n = cleanup_hits(rhts->rht, rhts->n, sn);
	double_hits(rhts);
	//this part purge haplotigs again
	qsort(rhts->rht, rhts->n, sizeof(eg_hit_t), cmp_q);	
	print_hits(rhts->rht, 0, rhts->n, sn, "bef");
	flt_hits4(rhts->rht, rhts->n, sn);
	rhts->n = cleanup_hits(rhts->rht, rhts->n, sn);
	qsort(rhts->rht, rhts->n, sizeof(eg_hit_t), cmp_qtgse);	
	print_hits(rhts->rht, 0, rhts->n, sn, "aft");
	
	//this part create undirected graph
	purge_dups2(rhts->rht, rhts->n, sn, &dups); //second round to purge continous query 
	//time to update dups 
	sdict_t *dup_n = sd_init();
	/*fprintf(stderr, "dups %d %p\n", __LINE__, dups.a);*/
	update_dup_cords(&dups, sn, dup_n);	
	qsort(dups.a, dups.n, sizeof(dup_t), cmp_dupt);	
	merge_dups(&dups);
	print_dups(&dups, dup_n);
	sd_destroy(dup_n);
	kv_destroy(dups);
	if (rhts) eg_destroy(rhts);	
	sd_destroy(sn);	

	/*fprintf(stderr,"[M::%s] parsing paf...\n", __func__);*/
	return 0;	
}




