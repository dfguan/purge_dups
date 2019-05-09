/*
 * =====================================================================================
 *
 *       Filename:  build_graph.h
 *
 *    Description:  header for build graph 
 *
 *        Version:  1.0
 *        Created:  19/11/2018 19:53:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef BUILD_GRAPH_H
#define BUILD_GRAPH_H

#ifdef __cplusplus
extern "C" {
#endif
int buildg(void *ctgs, void *g, char *edge_fn, int min_wt, char *prefix);
#ifdef __cplusplus
}
#endif

#endif
