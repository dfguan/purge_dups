/*
 * =====================================================================================
 *
 *       Filename:  cov.h
 *
 *    Description:  header 
 *
 *        Version:  1.0
 *        Created:  14/03/2019 17:22:03
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef COV_H
#define COV_H

#include <stdint.h>
#include <sys/types.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif
#ifdef __cplusplus
extern "C" {
#endif
void* read_covs(char *fn, void *sn);
void* read_cuts(char *fn);
#ifdef __cplusplus
}
#endif

#endif
