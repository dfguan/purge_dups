/*
 * =====================================================================================
 *
 *       Filename:  build_graph.c
 *
 *    Description:  build graph with links information
 *
 *        Version:  1.0
 *        Created:  19/11/2018 19:39:30
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
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "bed.h"
#include "graph.h"
#include "sdict.h"
#include "cdict.h"

sdict_t *col_ctgs(char *fn)
{
	bed_file_t *bf = bed_open(fn);
	if (!bf) return 0;
	sdict_t *ctgs = sd_init();		
	bed_rec_t r;
	while (bed_read(bf, &r) >= 0) 
		sd_put(ctgs, r.ctgn, 0, r.len, 1);
	bed_close(bf);
	return ctgs;
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

int get_links(char *links_fn, cdict_t *cds, sdict_t *ctgs)
{	
	bed_file_t *bf = bed_open(links_fn);
	if (!bf) return 0;
	lnk_rec_t r;
	uint32_t line_n = 0;
	name_t nt, nt2;
	while (lnk_read(bf, &r) >= 0) {
		uint32_t ind1, ind2;
		parse_name(r.ctgn, strlen(r.ctgn), &nt);
		parse_name(r.ctgn2, strlen(r.ctgn2), &nt2);
		if (!strcmp(nt.ctgn, nt2.ctgn)) continue; 
		ind1 = sd_get(ctgs, nt.ctgn);
		ind2 = sd_get(ctgs, nt2.ctgn);
		if (nt.e != ctgs->seq[ind1].len || nt2.e != ctgs->seq[ind2].len) continue;	
		/*if (r.is_l) */
			/*ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, 0);		*/
		/*else*/
			/*ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, 0, r.llen);		*/
		/*if (r.is_l2)*/
			/*sd_put2(ctgs, r.ctgn2, 0, 0, 0, r.rlen, 0);		*/
		/*else*/
			/*sd_put2(ctgs, r.ctgn2, 0, 0, 0, 0, r.rlen);		*/
		/*uint32_t ind2 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, r.rlen);		*/
		line_n += 1;
		cd_add2(&cds[ind1<<1|r.is_l], nt2.ctgn, r.is_l2, r.wt, r.llen);	//this has been normalized	
		cd_add2(&cds[ind2<<1|r.is_l2], nt.ctgn, r.is_l, r.wt, r.rlen);	//this has been normalized	
	} 
	bed_close(bf);
	return 0;	
}
graph_t *build_graph(cdict_t *cds, sdict_t *ctgs)
{
	graph_t *g = graph_init();
	
	uint32_t n = ctgs->n_seq << 1;
	uint32_t i, j;
	//create nodes
	for ( i = 0; i < ctgs->n_seq; ++i) 
		add_node(g, ctgs->seq[i].name, 0, ctgs->seq[i].len);
	//create edges	
	for ( i = 0; i < n; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint8_t is_l = i & 1;	
		cdict_t *c = cds + i;	
		if (!c->n_cnt) continue;	
		for (j = 0; j < c->lim; ++j) {
			char *name2 = c->cnts[j].name;
			if (strcmp(name1, name2) == 0) continue;
			//hsortand shaking
			/*fprintf(stderr, "try hand shaking\n");*/
			uint32_t ind = sd_get(ctgs, name2) << 1 | c->cnts[j].is_l;
			uint32_t k;
			uint8_t hand_shaking = 0;
			uint32_t ocnt = 0;
			for ( k = 0; k < cds[ind].lim; ++k) {
					if (strcmp(name1, cds[ind].cnts[k].name) == 0 && cds[ind].cnts[k].is_l == is_l) {
						hand_shaking = 1;
						ocnt = cds[ind].cnts[k].cnt;
						break;
					}
			}	
			/*if (hand_shaking) fprintf(stderr, "I hand shaking\n");*/
			ocnt = 1;
			if (hand_shaking) add_dedge(g, name1, is_l, name2, c->cnts[j].is_l, ocnt * c->cnts[j].cnt);	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.
		}		
	}	
	return g;
}



int buildg(void *sn, void *_og, char *edge_fn, int min_wt, char *prefix)
{
	sdict_t *ctgs = (sdict_t *) sn;
	graph_t *og = (graph_t *)_og; 
#ifdef DEBUG
	fprintf(stderr, "[M::%s] collecting contigs from faidx file\n", __func__);
#endif
	if (!ctgs) return 1;
	uint32_t n_ctg = ctgs->n_seq;	
	/*fprintf(stderr, "%u\n", ctgs->n_seq);*/
	cdict_t* cds = calloc(n_ctg<<1, sizeof(cdict_t)); 
	uint32_t n_cds = n_ctg<<1;
	uint32_t i;
	for ( i = 0; i < n_cds; ++i) cd_init(cds+i); 
#ifdef DEBUG
	fprintf(stderr, "[M::%s] collecting links\n", __func__);
#endif
	get_links(edge_fn, cds, ctgs);
	/*anothernorm(cds, ctgs);*/
	/*return 0;*/
	for ( i = 0; i < n_cds; ++i) cd_sort(cds+i); 
	cd_set_lim(cds, n_cds, min_wt); 
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
#ifdef DEBUG
	fprintf(stderr, "[M::%s] building graph\n", __func__);
#endif
	graph_t *g = build_graph(cds, ctgs);
#ifdef DEBUG
	fprintf(stderr, "[M::%s] processing graph\n", __func__);
#endif
	process_graph(g);
			
#ifdef DEBUG
	fprintf(stderr, "[M::%s] merging graph\n", __func__);
#endif
	merge_graph(og, g, 1);
#ifdef DEBUG
	fprintf(stderr, "[M::%s] output graph\n", __func__);
#endif

	for (i = 0; i < n_cds; ++i) {
		/*fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);*/
		cd_destroy(cds +i);	

	}
	/*fprintf(stderr, "leave\n");*/
	if (cds) free(cds);
#ifdef DEBUG
	fprintf(stderr, "[M::%s] releasing memory\n", __func__);
#endif
	return 0;

}

