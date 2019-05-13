/*
 * =====================================================================================
 *
 *       Filename:  graph.c
 *
 *    Description:  realization of graph functions  
 *
 *        Version:  1.0
 *        Created:  21/10/2018 12:55:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <zlib.h>

#include "graph.h"
#include "khash.h"
#include "ksort.h"
#include "kdq.h"
#include "kseq.h"
#include "kvec.h"



KDQ_INIT(uint32_t)


KSEQ_INIT(gzFile, gzread, gzseek)	
/*KSTREAM_INIT(gzFile, gzread, gzseek, 0x10000)*/



#define edge_key(a) ((a).v)
KRADIX_SORT_INIT(edge, edge_t, edge_key, 4)

KHASH_MAP_INIT_STR(str, uint32_t)

typedef khash_t(str) shash_t;

uint8_t rc_table[128]={
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,84,4,71,4,4,4,67,4,4,4,4,4,4,67,4,
	4,4,4,4,65,4,4,4,4,4,4,4,4,4,4,4,4,
	84,4,71,4,4,4,67,4,4,4,4,4,4,67,4,4,
	4,4,4,65,4,4,4,4,4,4,4,4,4,4,4
};

graph_t *graph_init(void) 
{	
	graph_t *g = (graph_t *)calloc(1, sizeof(graph_t));	
	g->h = kh_init(str);
	g->as.h = kh_init(str);
	return g;
}

void free_node(vertex_t *n)
{
	if (!n) return;
	if (n->name) free(n->name);
	if (n->seq) free(n->seq);
}


void graph_destroy(graph_t *g)
{
	/*fprintf(stderr, "enter\n");*/
	if (g) {
		//free nodes
		uint32_t n_nd = g->vtx.n;
		uint32_t i;
		for ( i = 0 ; i < n_nd; ++i) free_node(&g->vtx.vertices[i]);
		free(g->vtx.vertices);
		if (g->h) kh_destroy(str, (shash_t *)g->h);		
		//free edges
		if (g->eg.edges) free(g->eg.edges);
		if (g->eg.edge_idx) free(g->eg.edge_idx);
		//free paths
		if (g->pt.paths) {
			for (i = 0; i < g->pt.n; ++i) free(g->pt.paths[i].ns);
			free(g->pt.paths);
		} 
		//free asms
		if (g->as.h) kh_destroy(str, (shash_t *)g->as.h);
		if (g->as.asms) {
			for ( i = 0; i < g->as.n; ++i) free(g->as.asms[i].pn);
			free(g->as.asms);	
		}
		free(g);
	}

	/*fprintf(stderr, "leave\n");*/
}


int out_vetices(graph_t *g, FILE *fout)
{

	vertex_t *vs = g->vtx.vertices;
	uint32_t n_vs = g->vtx.n;
	uint32_t i;
	for ( i = 0; i < n_vs; ++i) 
		fprintf(fout, "S\t%s\t%s\t%u\n",  vs[i].name, vs[i].seq ? "ATG" : "*", vs[i].len);
	
	return 0;
}

int out_edges(graph_t *g, int all, FILE *fout)
{
	vertex_t *vs = g->vtx.vertices;
	edge_t *edg = g->eg.edges;
	path_t *pt = g->pt.paths;
	uint32_t n_edges = g->eg.n;		
	if (all) n_edges += g->eg.n_del;	
	uint32_t i;
	for ( i = 0; i < n_edges; ++i) {
		uint32_t v = edg[i].v, w = edg[i].w;
		if (v>>2 > w >> 2) fprintf(fout, "L\t%s\t%c\t%s\t%c\t%s\twt:%u\n", v>>1 & 1 ? pt[v>>2].name : vs[v>>2].name, v&1?'+':'-', w>>1 & 1 ? pt[w>>2].name : vs[w>>2].name, w&1?'+':'-', "*", edg[i].wt); // + head of sequence - tail of sequqnce	
	}
	return 0;
}

int out_paths(graph_t *g, FILE *fout)
{
	vertex_t *vs = g->vtx.vertices;
	path_t *p = g->pt.paths;
	uint32_t n_p = g->pt.n;
	uint32_t i;	
	for ( i = 0; i < n_p; ++i) {
		uint32_t j;
		if (!p[i].name) 
			fprintf(fout, "P\t%c%06u\t", p[i].is_circ?'c':'u',i); //check before output use the same name?
		else 
			fprintf(fout, "P\t%s\t", p[i].name);
		uint32_t v;	
		/*fprintf(stderr, "%d\n", p[i].n);*/
		for ( j = 0; j < p[i].n; ++j)  {
			v = p[i].ns[j];
			fprintf(fout, "%s%c%c", v>>1 & 1 ? p[v>>2].name : vs[v>>2].name, v&1?'+':'-', j + 1 == p[i].n ? '\n' : ','); // + head of sequence - tail of sequqnce	
		}
	}
	return 0;
}

int out_asms(graph_t *g, FILE *fout)
{
	asm_t *as = g->as.asms;
	path_t *ps = g->pt.paths;
	uint32_t n_a = g->as.n;
	uint32_t i, j;
	for ( i = 0; i < n_a; ++i) {
		fprintf(fout, "A\t%s\t%u\t", as[i].name, as[i].n);
		uint32_t t;
		for ( j = 0; j < as[i].n; ++j ) {
			t = as[i].pn[j];
			fprintf(fout, "%s%c",ps[t>>1].name, j + 1 == as[i].n ? '\n':',');
		}
	}
	fprintf(fout, "C\t%s\n", as[g->as.casm].name);
	return 0;
}

int out_graph(graph_t *g)
{
	fprintf(stdout, "H\tVN:Z:1.0\n");	
	out_vetices(g, stdout);
	out_edges(g, 0, stdout);
	out_paths(g, stdout);
	return 0;
}

uint32_t get_name2id(graph_t *g, char *nm)
{
	shash_t *h = (shash_t *)g->h;
	khint_t k = kh_get(str, h, nm);
	return k == kh_end(h) ? -1 : kh_val(h, k);
}

int update_seqs(graph_t *g, char *name, char *s, uint32_t l)
{
	vertices_t * ns = &g->vtx;
	uint32_t idx = get_name2id(g,  name); 
	vertex_t *n = &ns->vertices[idx>>1];
	n->seq = s;	
	n->len = l;	
	return 0;
}

//make mistakes when break the order S->P->L 
uint32_t add_node(graph_t *g, char* name, char *seq, uint32_t len)
{
	shash_t *h = (shash_t *)g->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, name, &absent);
	if (absent) {
		vertices_t * ns = &g->vtx;
		if (ns->n == ns->m) {
			ns->m = ns->m ? ns->m << 1 : 16;
			ns->vertices = realloc(ns->vertices, sizeof(vertex_t) * ns->m);	
		}
		vertex_t *n = &ns->vertices[ns->n]; 
		if (seq) n->seq = strdup(seq);
		else n->seq = 0;
		n->len = len;
		kh_key(h, k) =n->name = strdup(name);	
		kh_val(h, k) = (ns->n<<1 | 0); //0 to identify vertices
		++ns->n;	
	} else {
		uint32_t ind = kh_val(h, k);
		if (ind & 1) {
			//this is a path
			/*fprintf(stderr, "[W::%s] %s is a path\n", __func__, name);*/
			return ind;		
		} 
		ind = ind >> 1;
		if (seq && !g->vtx.vertices[ind].seq) g->vtx.vertices[ind].seq = strdup(seq);
		if (len && !g->vtx.vertices[ind].len) g->vtx.vertices[ind].len = len; 	
	}	
	return kh_val(h, k);
}
//only used for new graph don't mix path with vertex
int idx_edge(graph_t *g)
{
	/*fprintf(stderr, "indexing edges\n");*/
	uint32_t n_vtx = g->vtx.n;
	if (!g->eg.edge_idx) g->eg.edge_idx = calloc(n_vtx<<1, sizeof(uint64_t));
	
	edges_t *es = &g->eg;
	uint32_t n_edges = es->n;
	if (!es->is_srt) {
		radix_sort_edge(es->edges, es->edges+n_edges);
		es->is_srt = 1;
	} 
	uint64_t *idx = g->eg.edge_idx;
	uint32_t i, j;
	/*fprintf(stderr, "enter %u edges\n",n_edges);*/
	if (n_edges) 	
		for ( i = 1, j = 0; i <= n_edges; ++i) {
				/*fprintf(stderr, "%u\n", i);*/
			if (i == n_edges || es->edges[i].v != es->edges[j].v) {
				/*fprintf(stderr, "%u\n", es->edges[j].v);*/
				idx[vtx_idx(es->edges[j].v)] = (uint64_t) j << 32 | (i - j);	
				j = i;	
			}
		}
	/*fprintf(stderr, "leave\n");*/
}

int add_edge1(graph_t *g, edge_t *e)
{
	edges_t *es = &g->eg;
	if (es->n == es->m) {
		es->m = es->m ? es->m << 1 : 16;
		es->edges = realloc(es->edges, sizeof(edge_t) * es->m);
	}
	es->edges[es->n++] = *e;
}

int add_udedge(graph_t *g, char *sname, uint32_t sl, char *ename, uint32_t er, uint32_t wt)
{
	uint32_t sind = add_node(g, sname, 0, 0);
	uint32_t eind = add_node(g, ename, 0, 0);
	
	edge_t e = (edge_t) {sind << 1 | sl, eind << 1 | er, wt, 0, 0};
	edge_t re = (edge_t) {eind << 1 | er, sind << 1 | sl, wt, 0, 0}; //undirect graph
	
	add_edge1(g, &e);
	add_edge1(g, &re);
	return 0;
}


int add_dedge(graph_t *g, char *sname, uint32_t sl, char *ename, uint32_t er, uint32_t wt)
{
	uint32_t sind = add_node(g, sname, 0, 0);
	uint32_t eind = add_node(g, ename, 0, 0);
	edge_t e = (edge_t) {sind << 1 | sl, eind << 1 | er, wt, 0, 0};
	/*edge_t re = (edge_t) {eind << 1 | er, sind << 1 | sl, wt, 0, 0}; //undirect graph*/
	/*fprintf(stderr, "%s %s\n", sname, ename);	*/
	/*fprintf(stderr, "EDGE %u\t%u\n",sind<<1|sl, eind << 1 | er);	*/
	add_edge1(g, &e);
	/*add_edge1(g, &re);*/
	return 0;
}


uint32_t add_path(graph_t *g, char *name,  uint32_t *nodes, uint32_t n, uint32_t is_circ) 
{
	
	paths_t *ps = &g->pt;
	shash_t *h = (shash_t *)g->h;
	char *pname = name;
	uint32_t pn = ps->n;
	if (!pname) {
		//this is a new path it's name hasn't been initiated yet
		pname = malloc(sizeof(char) * 8); // path name length is 7;
		do {
			sprintf(pname, "%c%06u%c", is_circ ? 'c':'u', pn, 0);	
			++pn;	
		} while (kh_get(str, h, pname) != kh_end(h)); 
	} else {
		pname = malloc(sizeof(char) * (strlen(name) + 1));
		strcpy(pname, name);
	}
	
	khint_t k;
	int absent;
	k = kh_put(str, h, pname, &absent);
	if (absent) {
		if (ps->n == ps->m) {
			ps->m = ps->m ? ps->m << 1 : 16;
			ps->paths = realloc(ps->paths, sizeof(path_t) * ps->m);
		}	
		path_t *p = &ps->paths[ps->n];
		p->ns = malloc(n*sizeof(uint32_t));
		memcpy(p->ns, nodes, sizeof(uint32_t) * n);
		p->n = n;
		kh_key(h, k) =p->name = strdup(pname);	
		kh_val(h, k) = (ps->n<<1 | 1);
		++ps->n;	
	} else 
		fprintf(stderr, "[W::%s] path name %s has been added\n", __func__, pname);
	if (pname) free(pname);
	return kh_val(h, k);
}

int add_asm(graph_t *g, char *name,  uint32_t *nodes, uint32_t n) 
{
	asms_t *as = &g->as;
	shash_t *h = (shash_t *)as->h;
	
	char *aname = name;
	uint32_t an = as->n;
	if (!aname) {
		//this is a new path it's name hasn't been initiated yet
		aname = malloc(sizeof(char) * 7);// path name length is 7
		do {
			sprintf(aname, "a%05u%c", an, 0);	
			++an;	
		} while (kh_get(str, h, aname) != kh_end(h)); 
	} else {
		aname = malloc(sizeof(char) * (strlen(name) + 1));
		strcpy(aname, name);
	}

	khint_t k;
	int absent;
	k = kh_put(str, h, aname, &absent);
	if (absent) {
		if (as->n == as->m) {
			as->m = as->m ? as->m << 1 : 16;
			as->asms = realloc(as->asms, sizeof(asm_t) * as->m);
		}	
		asm_t *a = &as->asms[as->n];
		a->pn = malloc(n*sizeof(uint32_t));
		memcpy(a->pn, nodes, sizeof(uint32_t) * n);
		a->n = n;
		kh_key(h, k) =a->name = strdup(aname);	
		kh_val(h, k) = as->n++;
	} else 
		fprintf(stderr, "[W::%s] assembly has been added\n", __func__);
	if (aname) free(aname);
	return kh_val(h, k);
}

int del_r(graph_t *g, edge_t *a) // remove reverse edge
{
	uint32_t v_idx = vtx_idx(a->w);
	edge_t *e = edges(g, v_idx);
	uint32_t ne = edge_n(g, v_idx);
	edge_t *i;
	for ( i = e; i < e + ne; ++i) {
		if (i->w == a->v) {
			i->is_del = 1;
			break;
		}
	}
	return 0;
}

//only keeps maximum weighted edges
int clean_edges(graph_t *g)
{	
	uint32_t n_vtx = g->vtx.n; 
	uint32_t i;
	/*out_edges(g);*/
	uint32_t n_del = 0;
	for (i = 0; i < n_vtx<<1; ++i) {
		edge_t *e = edges(g, i);
		uint32_t en = edge_n(g, i);	
		edge_t *a;
		//find maximum weight
		//if there are mulitple maximum weights break them all we can't decide which way to choose
		uint32_t mwt = 0;
		uint32_t n_mwt = 0;	
		for ( a = e; a < e + en; ++a) {
			if (a->is_del) continue;
			if (a->wt > mwt) {
				mwt = a->wt;
				n_mwt = 1;	
			} else if (a->wt == mwt) {
				++n_mwt;
			} 	
		} 

		for (a = e; a < e + en; ++a) {
			if (a->is_del) continue;
			if (n_mwt > 1 || a->wt != mwt) {
				++n_del;
				a->is_del = 1;
				del_r(g, a);			
			} 		
		}
	}	
	g->eg.n_del = n_del;
	/*fprintf(stderr, "%u nodes %u del\n", n_vtx, n_del);	*/
	return 0;
}

int join_ends(graph_t *g)
{
	uint32_t n_vtx = g->vtx.n;
	
	uint32_t i;
	for ( i = 0; i < n_vtx << 1; i+=2) {
		edge_t e = (edge_t) {(i<<1), (i<<1)+1, 1000, 0, 0};	
		add_edge1(g, &e);
		edge_t re = (edge_t) {(i<<1)+1, i<<1, 1000, 0, 0};
		add_edge1(g, &re);
	}
	fprintf(stderr, "[M::%s] %u edges\n", __func__, g->eg.n);
	g->eg.is_srt = 0; //set unsorted
	return 0;
}

int update_graph(graph_t *g) // update index
{
	uint32_t n_edges = g->eg.n;		
	edge_t *es = g->eg.edges;
	edge_t *end = es + n_edges; 
	uint32_t n_del = 0;
	
	for (; es < end; ++es) {
		if (es->is_del) {
			edge_t tmp = *es;
			*es = *(end - 1);
			*(end - 1) = tmp;
			--end;
			--es;	
			++n_del;
		}	
	}
	fprintf(stderr, "[M::%s] %u vertices\n", __func__, g->vtx.n);	
	fprintf(stderr, "[M::%s] %u edges, %u del\n", __func__, n_edges, n_del);	
	g->eg.n -= n_del;
	idx_edge(g);
	return 0;
}


int vis_r(graph_t *g, edge_t *a)
{
	edge_t *e = edges(g, vtx_idx(a->w));
	uint32_t ne = edge_n(g, vtx_idx(a->w));
	edge_t *i;
	for ( i = e; i < e + ne; ++i) {
		if (i->w == a->v) {
			i->is_vis = 1;
			break;
		}
	}
	return 0;
}

//no paths graph
int srch_path(graph_t *g)
{
	fprintf(stderr, "[M::%s] traversing graph\n", __func__);	
	uint32_t n_vtx = g->vtx.n;
	/*fprintf(stderr, "%d\n", n_vtx);		*/
	uint8_t *mark = calloc(n_vtx << 1, sizeof(uint8_t));
	
	uint32_t i;
	
	kdq_t(uint32_t) *q;
	q = kdq_init(uint32_t);
	uint8_t is_circ;
	for ( i = 0; i < n_vtx<<1; ++i) {
		if (mark[i]) continue;
		uint32_t p, s;
		p = s = i; kdq_clean(q); is_circ = 0;
		while (1) {
			if (mark[p]) break;
			else {
				mark[p] = 1;
				/*fprintf(stderr, "%d\n", p);*/
				kdq_push(uint32_t, q, p);				
			}
		 //traverse forwardly
			
			edge_t *es = edges(g, p);
			uint32_t es_n = edge_n(g, p);
			uint32_t j;
			for (j = 0; j < es_n; ++j) if (!(es + j)->is_vis) break;
			/*fprintf(stderr, " leave \n");*/
			if (j == es_n) break;
			else {
				es = es + j;
				es->is_vis = 1;
				vis_r(g, es);
				p = vtx_idx(es->w);
			}
		}
			/*fprintf(stderr, "ENTER %d %d\n", __LINE__, kdq_size(q));*/
		if (p != s || kdq_size(q) == 1) {
			p = s;
			while (1) { // traverse backwardly 
				edge_t *es = edges(g, p);
				uint32_t es_n = edge_n(g, p);
				uint32_t j;
				for (j = 0; j < es_n; ++j) if (!(es + j)->is_vis) break;
				/*while (j < es_n) if (!(es + j)->is_vis) break;*/
				if (j == es_n) break;
				else {
					es = es + j;
					es->is_vis = 1;
					vis_r(g, es);
					p = vtx_idx(es->w);
				}
				if (mark[p]) break;
				else {
					mark[p] = 1;
					/*fprintf(stderr, "%d\n", p);*/
					kdq_unshift(uint32_t, q, p);
				}	
			}	
		} else is_circ = 1;
		//add path
			/*fprintf(stderr, "ENTER %d %d\n", __LINE__, kdq_size(q));*/
		if (kdq_size(q) >= 2) {
			paths_t *pths = &g->pt;
			if (pths->n == pths->m) {
				pths->m = pths->m ? pths->m << 1 : 16;
				pths->paths = realloc(pths->paths, sizeof(path_t) * pths->m);
			}
			path_t *pth = &pths->paths[pths->n++];
			pth->name = 0;
			pth->ns = malloc(sizeof(uint32_t) * (kdq_size(q) >> 1));	
			/*fprintf(stderr, "ENTER %d %d\n", __LINE__, kdq_size(q));*/
			uint32_t k, m;
			for ( k = 0, m = 0; k < kdq_size(q); k += 2, ++m) pth->ns[m] = kdq_at(q, k);
			/*for ( k = 0; k < m; ++k) fprintf(stderr, "PID: %d\t%d\n", k, pth->ns[k]);*/
			pth->n = m;
			pth->is_circ = is_circ;
		}
		/*fprintf(stderr, "LEAVE %d\n", __LINE__);*/
	}
	return 0;
}

int merge_graph(graph_t *g, graph_t *c, int all)
{
	vertex_t *vs = c->vtx.vertices;
	uint32_t n_vs = c->vtx.n;
	
	uint32_t i;
	uint32_t *crspid = (uint32_t *)malloc(sizeof(uint32_t) * n_vs);
	for ( i = 0; i < n_vs; ++i) 
		crspid[i] = (add_node(g, vs[i].name, vs[i].seq, vs[i].len)); //could return a path id 
	edge_t *edg = c->eg.edges;
	uint32_t n_edges = c->eg.n;		
	if (all) n_edges += c->eg.n_del;	
	for ( i = 0; i < n_edges; ++i) {
		uint32_t v = edg[i].v, w = edg[i].w;
		add_dedge(g, vs[v>>2].name, v&1, vs[w>>2].name, w&1, edg[i].wt);//
	}
	/*fprintf(stderr, "path\n");	*/
	path_t *p = c->pt.paths;
	uint32_t n_p = c->pt.n;
	kvec_t(uint32_t) pids;
	kv_init(pids);

	for ( i = 0; i < n_p; ++i) {
		uint32_t j;
		/*fprintf(stderr, "old: %d\t%d\n", j, p[i].n);*/
		/*for (j = 0; j < p[i].n; ++j) fprintf(stderr, "old: %d\t%d\t%d\n", j, p[i].ns[j], crspid[p[i].ns[j]>>1]);*/
		for (j = 0; j < p[i].n; ++j) p[i].ns[j] = (crspid[p[i].ns[j]>>1] << 1) | (p[i].ns[j] & 1);
		/*for (j = 0; j < p[i].n; ++j) fprintf(stderr, "new: %d\t%d\n", j, p[i].ns[j]);*/
		/*fprintf(stderr, "new: %d\t%d\n", j, p[i].n);*/
		uint32_t pid = add_path(g, 0, p[i].ns, p[i].n, p[i].is_circ);
		kv_push(uint32_t, pids, pid);
	}
	//add asm
	/*fprintf(stderr, "asm\n");	*/
	uint32_t aid = add_asm(g, 0,  pids.a, pids.n);	
	g->as.casm = aid;
	graph_destroy(c);
	kv_destroy(pids);
	if (crspid) free(crspid);
	return 0;
}




int process_graph(graph_t *g)
{
	idx_edge(g);
	/*out_edges(g);*/
	clean_edges(g);
	join_ends(g);
	update_graph(g);
	srch_path(g);
	return 0;
}

// GFA IO
int add_s(graph_t *g, char *s)
{	
	/*fprintf(stderr, "enters\n");*/
	char *name = 0;
	char *seq = 0;
	char *p, *q;
	int i;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) name = q;
			else if (i == 1) {
				seq = q[0] == '*' ? 0 : strdup(q);
				break;
			}	
			++i, q = p + 1;	
			if (c == 0) break;	
		}
	}	
	uint32_t len = seq == 0 ? 0 : strlen(seq);
	add_node(g, name, seq, len);
	/*fprintf(stderr, "leaves\n");*/
	return 0;

}

int add_e(graph_t *g, char *s)
{
	char *n1, *n2;
	char d1, d2;
	char *p, *q;
	int i;
	uint32_t wt;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) n1 = q;
			else if (i == 1) d1 = q[0];
		   	else if (i == 2) n2 = q;
			else if (i == 3) d2 = q[0];
			else if (i == 5) wt = strtoul(q+3, NULL, 10);	
			++i, q = p + 1;	
			if (c == 0) break;	
		}
	}	

	add_dedge(g, n1, d1 == '+', n2, d2 == '+', wt);
	return 0;
}

int add_p(graph_t *g, char *s)
{
	char *name;	
	char *p, *q;
	int i;
	kvec_t(uint32_t) ns;
	kv_init(ns);
	char *nodes_str;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) name = q;
			else if (i == 1) nodes_str = q;
			++i, q = p + 1;	
			if (c == 0) break;	
		}
	}	
	for (p = q = nodes_str;; ++p) {
		int e = *p; 
		if (*p == 0 || *p == ',') {
			*p = 0;
				/*int c = q[ql - 1];*/
			int c =*(p-1);
			*(p-1) = 0;
			/*q[ql-1] = 0;*/
			uint32_t n_id = add_node(g, q, 0, 0);
			/*fprintf(stderr, "node: %s  %u\n",q, n_id);*/
			n_id = n_id << 1 | (c == '+');
			kv_push(uint32_t, ns, n_id);
			q = p + 1;
			
		}
	  	if (!e)	break;
	}
	/*fprintf(stderr, "enterp %d\n", ns.n);*/
	
	if (ns.n) add_path(g, name, ns.a, ns.n, name[0]=='c');
	kv_destroy(ns);
	/*fprintf(stderr, "leavep\n");*/
	return 0;
}

int add_a(graph_t *g, char *s)
{
	char *name;	
	char *p, *q;
	int i;
	kvec_t(uint32_t) ns;
	kv_init(ns);
	char *nodes_str;
	int node_n; 
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) name = q;
			else if (i == 1) node_n = atoi(q);
			else if (i == 2) nodes_str = q;
			++i, q = p + 1;	
			if (c == 0) break;	
		}
	}	
	for (p = q = nodes_str;; ++p) {
		int e = *p; 
		if (*p == 0 || *p == ',') {
			*p = 0;
				/*int c = q[ql - 1];*/
			/*int c =*(p-1);*/
			/**(p-1) = 0;*/
			/*q[ql-1] = 0;*/
			uint32_t n_id = get_name2id(g, q);
			/*fprintf(stderr, "node: %s  %u\n",q, n_id);*/
			/*n_id = n_id << 1 | (c == '-');*/
			kv_push(uint32_t, ns, n_id);
			q = p + 1;
			
		}
	  	if (!e)	break;
	}
	/*fprintf(stderr, "enterp %d\n", ns.n);*/
	if (ns.n) add_asm(g, name, ns.a, ns.n);
	kv_destroy(ns);
	/*fprintf(stderr, "leavep\n");*/
	return 0;
}

graph_t  *load_gfa(char *fn)
{
	graph_t *g = graph_init();
	
	kstream_t *ks;
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	kstring_t buf = {0, 0, 0};
	int dret;
	while (ks_getuntil(ks, KS_SEP_LINE, &buf, &dret) >= 0) {
		if (buf.s[0] == 'S') add_s(g, buf.s);
		else if (buf.s[0] == 'L') add_e(g, buf.s);	
		else if (buf.s[0] == 'P') add_p(g, buf.s);
	}
	return g;
}

int read_seq(graph_t *g, char *seqfn)
{
	gzFile fp = seqfn && strcmp(seqfn, "-") ? gzopen(seqfn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 1;
	kseq_t *seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) 
		add_node(g, seq->name.s, seq->seq.s, seq->seq.l);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}


int cp_seq(char *s, char *t, uint32_t len, int is_not_rc)
{
	if (is_not_rc) 
		memcpy(s, t, len * sizeof(char));
	else {
		uint32_t i;
		for ( i = 0; i < len; ++i) s[i] = rc_table[t[len -i - 1]];	
	}
	return 0;
}

uint32_t *parse_path(graph_t *g, uint32_t pid, uint32_t *n)
{
	path_t *pts = g->pt.paths;
	kvec_t(uint32_t) eles, ctgids;
	kv_init(eles);
	kv_init(ctgids);
	kv_push(uint32_t, eles, pid<<1);
	/*fprintf(stderr, "Enter %d %d\n", __LINE__, pts->n);*/
	while (eles.n > 0) {
		pid = kv_pop(eles);		
	/*fprintf(stderr, "Enter %d %s %s %d %d\n", __LINE__, pts[pid>>2].name, (pid>>1) & 1 ? "PATH" : "NODE", pts[pid>>2].n, pid);*/
	/*int j;	*/
	/*for (j = 0; j < pts[pid>>2].n; ++j) */
		/*fprintf(stderr, "Enter %d %s %d\n", __LINE__, (pts[pid>>2].ns[j]>>1) & 1 ? "PATH" : "NODE", pts[pid>>2].ns[j]);*/
		
		if ((pid >> 1) & 1) {
			//is a path
			path_t *pt = &pts[pid>>2];
			int i;
			if (pid & 1)  // == reverse complementary 
				for (i = 0; i < pt->n; ++i) 
					kv_push(uint32_t, eles, pt->ns[i]^1);	
			 else  //forward
				for ( i = pt->n - 1; i >= 0; --i) 
					kv_push(uint32_t, eles, pt->ns[i]);	
		}else  // is a seq
			kv_push(uint32_t, ctgids, pid);	
	}	
	/*fprintf(stderr, "LEAVE %d\n", __LINE__);*/
	kv_destroy(eles);
	*n = kv_size(ctgids);
	return ctgids.a;
}

int get_path(graph_t *g, uint32_t min_l, char *fn)
{	
	path_t *ps = g->pt.paths;
	vertex_t *vs = g->vtx.vertices;
	asm_t *as = &g->as.asms[g->as.casm];
	uint32_t i, j;
	FILE *fout = fn ? fopen(fn, "w") : stdout;
	for ( i = 0; i < as->n; ++i) {
		uint32_t m;
		uint32_t *p = parse_path(g, as->pn[i], &m);
	   	char *ref_nm = ps[as->pn[i]>>1].name;
		uint32_t ref_len = 0; 	
		for ( j = 0; j < m; ++j) {
			uint32_t seq_len = vs[p[j] >> 2].len;
			if (seq_len) {
				ref_len += seq_len;//200 'N' s	
				if (j != m - 1) ref_len += 200;
			} 		
		}
		char *ref_seq = NULL;
		if (ref_len && ref_len >= min_l) {
			ref_seq = malloc(sizeof(char) * (ref_len+1));
			char *s = ref_seq;
		
			for (j = 0; j < m; ++j) {
				cp_seq(s, vs[p[j]>> 2].seq, vs[p[j]>>2].len, p[j] & 1) , s += vs[p[j]>>2].len;
				if (j != m - 1) memset(s,'N', 200), s += 200;
			}
			*s = 0;
		} 	
		fprintf(fout, ">%s %u\n",ref_nm, ref_len);
		if (ref_seq) {
			fprintf(fout, "%s\n",ref_seq);
			free(ref_seq);		
		}
		free(p);
	}
	if (fn) fclose(fout); 
	return 0;
}

int set_c(graph_t *g, char *s)
{
	int i;
	char *p, *q, *name;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) name = q;
			++i, q = p + 1;	
			if (c == 0) break;	
		}
	}	
	if (i == 0) {
		fprintf(stderr, "[E::%s] current assembly is not set\n", __func__);
		return 1;
	}
	shash_t *h = (shash_t *)g->as.h;
	khint_t k = kh_get(str, h, name);
	if (k == kh_end(h)) {
		fprintf(stderr, "[E::%s] wrong assembly id \n", __func__);
		return 1;			
	}
   	g->as.casm = kh_val(h, k);
	return 0;
}

graph_t *load_sat(char *fn) 
{
	graph_t *g = graph_init();
	kstream_t *ks;
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	kstring_t buf = {0, 0, 0};
	int dret;
	while (ks_getuntil(ks, KS_SEP_LINE, &buf, &dret) >= 0) {
		if (buf.s[0] == 'S') add_s(g, buf.s);
		else if (buf.s[0] == 'L') add_e(g, buf.s);	
		else if (buf.s[0] == 'P') add_p(g, buf.s);
		else if (buf.s[0] == 'A') add_a(g, buf.s);
		else if (buf.s[0] == 'C') set_c(g, buf.s);
	}
	return g;	
}

int dump_sat(graph_t *g, char *fn)
{
	FILE *fout = fn ? fopen(fn, "w") : stdout;
	fprintf(fout, "H\tVN:Z:1.0\n");	
	out_vetices(g, fout);
	/*fprintf(stderr, "output paths\n");*/
	out_paths(g, fout);
	/*fprintf(stderr, "output asms\n");*/
	out_edges(g, 0, fout);
	out_asms(g, fout);
	if (fn) fclose(fout);
	return 0;
}
