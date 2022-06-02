/*
 * =====================================================================================
 *
 *       Filename:  calcuts.c
 *
 *    Description:  calculate coverage for haploid and diploid
 *
 *        Version:  1.0
 *        Created:  13/03/2019 20:43:10
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
#include <string.h>
#include <math.h>
#include <zlib.h>
#include <getopt.h>

#include "kvec.h"

#define swap(T, x, y) do {T tmp = x; x = y; y = tmp; } while(0)
#define MAX_DEPTH  500
#define LOWEST_CUT 5

typedef struct {
	float min_frac;
	int   min_mc;
	char *stat_fn;
	int isdip;
	int low_cov;
	int upper_cov;
	int dip_cov;
	int fhord;
} opt_t;

typedef struct {
	int idx_s:31, ispeak:1, idx_e: 31, del:1;
	uint32_t cnt;	
} locopt_t;
float norm_cdf(int x, float p, int n) {
    float mean = n * p;
    float sd = sqrt(n * p * (1 - p));
    return 0.5 * (1 + erf((x - mean)/(sd * sqrt(2))));
}
//read depth counts from a bedlike file produced by pbcstat (DG).   
uint32_t *read_counts(const char *fn)
{
	uint32_t max_depth = MAX_DEPTH;
	uint64_t *depth2cnt = (uint64_t *) calloc(max_depth + 1, sizeof(uint64_t));
	if (!depth2cnt) return 0;
	FILE *fp = strcmp(fn, "-") == 0 ? stdin : fopen(fn, "r");
	int idx;
	uint64_t cnt, mcnt;
    mcnt = 0;
	while (EOF != fscanf(fp, "%d\t%lu\n", &idx, &cnt)) {
		depth2cnt[idx] = cnt;
        if (cnt > mcnt)
            mcnt = cnt;
    }
    if (mcnt > UINT32_MAX) {
        int scale = 0;
        while (mcnt > UINT32_MAX) {
            ++scale;
            mcnt >>= 1;
        }
        for (idx = 0; idx <= max_depth; ++idx)
            depth2cnt[idx] >>= scale;
    }
	if (strcmp(fn, "-") != 0) fclose(fp);	
	uint32_t *cnts = (uint32_t *) calloc(max_depth + 1, sizeof(uint32_t));
    for (idx = 0; idx <= max_depth; ++idx)
        cnts[idx] = (uint32_t) depth2cnt[idx];
    free(depth2cnt);
    return cnts;
}

int get_mean(uint32_t *depth2cnt)
{
	int i;
	uint64_t sum = 0, totalc = 0;
	for (i = 1; i < MAX_DEPTH; ++i) {
		sum += (uint64_t) depth2cnt[i]* i;
		totalc += depth2cnt[i];	
	}
	return sum/totalc;
}

int get_max(uint32_t *depth2cnt)
{
	int max_idx;
	uint32_t max_cnt = 0;
	int i;
	for ( i = LOWEST_CUT + 1; i < MAX_DEPTH; ++i) 
		if (depth2cnt[i] > max_cnt) max_idx = i, max_cnt = depth2cnt[i];
	return max_idx;
}

int cmplocopt(const void *p, const void *q) 
{
	/*locopt_t *a = (locopt_t *)p;*/
	/*locopt_t *b = (locopt_t *)q;*/
	return (int)(((locopt_t *)q)->cnt - ((locopt_t *)p)->cnt);
	/*if (a->cnt > b->cnt) return -1;*/
	/*else if (a->cnt == b->cnt) return 0;*/
	/*else return 1;*/
}

int print_loopt(locopt_t *lopt, int n)
{
	int i;
	for ( i = 0; i < n; ++i) 
	fprintf(stderr, "%s\t%d\t%d\t%u\t%s\n",lopt[i].ispeak ? "PEAK" : "VALLEY", lopt[i].idx_s, lopt[i].idx_e, lopt[i].cnt, lopt[i].del? "del":"not del");
	return 0;
}


int outliers(uint32_t *depth2cnt)
{
	int i;
	uint32_t n = 0;
	for ( i = 1; i <= MAX_DEPTH; ++i) n += depth2cnt[i];
	uint32_t m = n >> 1;	
	uint32_t p[] = {m, n, m};
	uint32_t z[] = {m>>1, m, ((m>>1) + 1 + n) >> 1};		
	float quarts[] = {0, 0, 0};
	int j = 1; uint32_t num = 0;
	for ( i = 0; i < 3; ++i) {
		for ( ; j <= MAX_DEPTH; ++j) 
			if (z[i] > num) num += depth2cnt[j]; 	
			else break;
		if (p[i] & 1) 
			quarts[i] = j;	
		else if (z[i] + 1 <= num) 
			quarts[i] = j;
		else 
			quarts[i] = (float) j + .5;	
	}
	float iqr = quarts[2] - quarts[0];
	fprintf(stderr, "%f\t%f\t%f\n", quarts[0], quarts[1], quarts[2]);
	fprintf(stderr, "%f\t%f\n", quarts[0] - 1.5 * iqr, quarts[2] + 1.5 * iqr);	
	return 0;
}


/*
 * @fun        calcuts
 * @params     depth-count array
 * @alg        calculate peaks and possible valleys given a Gaussion Distribution through first and second derivatives 
 *
 */ 

#define pi 3.14
#define square(x) ((x)*(x))
#include <math.h>

uint32_t gaussion_val(double x, double m, double mu, double delta)
{
	double f = 1 / sqrt(2 * pi * delta) * exp(-square(x-mu)/(2*delta));
	return (uint32_t)(f * m);
}

int gaussion_est(uint32_t *depth2cnt, int s, int e, double *mu, double *delta)
{
	double sum = 0;
	double totalc = 0;
	int i;
	for ( i = s; i < e; ++i) {
		double t = depth2cnt[i]; 
		totalc += t;
		sum += t * i; // count * n
	}
		
	double sqsum = 0;
	double avg = sum / totalc;	
	for ( i = s; i < e; ++i) {
		double t = depth2cnt[i]; 
		sqsum += t * square((double)(i) - avg); // count * (n)^2
	}
	*mu = avg;
	*delta = sqsum /totalc;
	return 0;
}
int get_gaussion(uint32_t *depth2cnt, uint32_t max_idx, uint32_t max_cnt) 
{
	int i;	
	for ( i = 0; i <= MAX_DEPTH && (float)depth2cnt[i]/max_cnt < .2; ++i); 
		int s = i;
		for ( ; i <= MAX_DEPTH && (float)depth2cnt[i]/max_cnt >= .2; ++i); 
		int e = i;	
		double mu, delta;
		double tot_cnt = 0;
		/*s = 3;*/
		/*e = MAX_DEPTH + 2;*/
		gaussion_est(depth2cnt, s, e, &mu, &delta);
		fprintf(stderr, "%u \t %lf\t %lf\t %lf\t %lf\n", i, mu, delta, mu - 3 * sqrt(delta), mu + 3 * sqrt(delta));	
	return 0;	
}

int calcuts(uint32_t *depth2cnt, int *cutoffs, int min_mc, float min_frac, int fhord) 
{
	int max_idx = get_max(depth2cnt);
	uint32_t max_cnt = depth2cnt[max_idx];
	/*outliers(depth2cnt);*/
	/*get_gaussion(depth2cnt, max_idx, max_cnt);	*/
	/*uint32_t all_len = 0;	*/
	/*float *acc_cnt = malloc(sizeof(float) * (MAX_DEPTH));*/
	int i, j;
	/*for ( i = 1; i <= MAX_DEPTH; ++i) {*/
		/*all_len += depth2cnt[i];	*/
	/*}*/
	/*uint32_t tmp_len = 0;*/
	/*for (i = 0; i <= MAX_DEPTH; ++i) {*/
		/*tmp_len = 0;*/
		/*if ( i < max_idx)  for ( j = 1; j <=i;++j) tmp_len += depth2cnt[j];*/
		 /*else for ( j = i; j <= MAX_DEPTH; ++j) tmp_len += depth2cnt[j];*/
		/*acc_cnt[i] = (float) tmp_len / all_len;	*/
	/*}		*/
	/*for ( i = 1; i <= MAX_DEPTH; ++i) {*/
		/*fprintf(stderr, "%d\t%f\t%f\n", i, acc_cnt[i], (float)depth2cnt[i]/max_cnt);*/
	/*}*/
	
	
	//first derivatives
	int *derts1 = malloc(sizeof(int) * (MAX_DEPTH+1)); 
	int *derts1_fit = malloc(sizeof(int) * (MAX_DEPTH+1)); 
	derts1[0] = 1;
	for ( i = 0; i <= MAX_DEPTH; ++i) 
		derts1[i+1] = (int)(depth2cnt[i+1] - depth2cnt[i]); 	
	//fitting
	derts1_fit[0] = derts1[0];
	derts1_fit[1] = derts1[1];
	//local vote to change to new value upper 5 next 5
	derts1_fit[MAX_DEPTH] = derts1[MAX_DEPTH];	
	for ( i = 2; i < MAX_DEPTH - 1; ++i) {
		int pos = 0, neg = 0;
		int z;
		for ( z = i - 1; z > 0 && z > i - 6; --z) derts1[z] >= 0 ? ++pos : ++neg;	
		for ( z = i + 1; z < MAX_DEPTH && z < i + 6; ++z) derts1[z] >= 0 ? ++pos : ++neg;	
		if ((pos > neg && derts1[i] < 0) || (neg > pos && derts1[i] > 0))
			derts1_fit[i] = (derts1[i-1] + derts1[i+1]) >> 1; 	
		else
			derts1_fit[i] = derts1[i];
	} 
	//second derivatives is not being used 
	int *derts2 = malloc(sizeof(int) *MAX_DEPTH);
	for (i = 0; i < MAX_DEPTH; ++i) 
		derts2[i] = derts1[i+1] - derts1[i];
	
	kvec_t(locopt_t)  locopts;
	kv_init(locopts);
	/*for ( i = 0; i <= MAX_DEPTH; ++i) */
		/*fprintf(stderr, "F1: %d\t%d\n",i, derts1[i]);*/
	/*for ( i = 0; i <= MAX_DEPTH; ++i) */
		/*fprintf(stderr, "F1: %d\t%d\n",i, derts1_fit[i]);*/
	/*for ( i = 0; i < MAX_DEPTH; ++i) */
		/*fprintf(stderr, "F2: %d\t%d\n",i, derts2[i]);*/
	//calculate peaks and valley in a linear fashion
	/*for (i = 3; i < MAX_DEPTH; ++i) */
		/*if ((derts1[i]^derts1[i+1]) <= 0 && (derts1[i+1] < 0 || derts1[i] < 0) && (float)depth2cnt[i]/max_cnt >= min_frac) {*/
			/*locopt_t tmp = (locopt_t){i, i, depth2cnt[i]};*/
			/*kv_push(locopt_t, locopts, tmp);*/
		/*}*/
	int peakn = 0;
	for (i = 3; i < MAX_DEPTH; ++i) 
		if ((derts1_fit[i]^derts1_fit[i+1]) <= 0 && (derts1_fit[i+1] < 0 || derts1_fit[i] < 0) && (float)depth2cnt[i]/max_cnt >= min_frac) {
			int isp = 0;
			if (derts1_fit[i+1] < 0) {
				++peakn;
				isp = 1;	
			}
			locopt_t tmp = (locopt_t){i, isp, i, 0, depth2cnt[i]};
			kv_push(locopt_t, locopts, tmp);
		}
	fprintf(stderr, "[M::%s] Find %d peaks\n", __func__, peakn);
	/*fprintf(stderr, "Find_Peaks: %d\n", peakn);		*/
	/*print_loopt(locopts.a, locopts.n);	*/
	/*fprintf(stderr, "After Processing\n");*/
	//sort and output	
	if (locopts.n >= 3) {
		/*int mean = get_mean(depth2cnt);*/
		/*fprintf(stderr, "%d\n", mean);*/
		//merge near values	
		int i, j, del = 0;
		for (i = 0, j = 1; j < locopts.n; ++j) {
			/*if (locopts.a[j].idx_s < locopts.a[i].idx_e + 3 && (float)locopts.a[j].cnt/locopts.a[i].cnt >0.95) {*/
			if (locopts.a[j].idx_s <= locopts.a[i].idx_e + 3) { // or better to compare their freqs? 
			if (locopts.a[j].cnt > locopts.a[i].cnt) 
				locopts.a[i].cnt = locopts.a[j].cnt; 
				locopts.a[i].idx_e = locopts.a[j].idx_s;	
				locopts.a[j].del = 1;
			} else 
				i = j;
		}
		//peak and valley assignment  
		for ( i = 0; i < locopts.n; ++i) {
				int is = locopts.a[i].idx_s;
				int ie = locopts.a[i].idx_e;
				if (is == ie) continue;
				if (derts1_fit[is] * derts1_fit[ie + 1] <= 0) {
					if (derts1_fit[is] < 0) 
						locopts.a[i].ispeak = 0;
					else
						locopts.a[i].ispeak = 1;
				} else 
					locopts.a[i].del = 1; // is local peaks  
		}
		
		peakn = 0;
		for (i = 0; i < locopts.n; ++i) 
			if (!locopts.a[i].del && locopts.a[i].ispeak) ++peakn;
		fprintf(stderr, "[M::%s] Merge local peaks and valleys: %d peaks remain\n", __func__, peakn);
		//would not like any peak/valleys less than LOWEST_CUT
		for ( i = 0; i < locopts.n; ++i) 
			if (locopts.a[i].idx_s < LOWEST_CUT) {
				if (locopts.a[i].ispeak) --peakn;	
				locopts.a[i].del = 1;
			}
		fprintf(stderr, "[M::%s] Remove peaks and valleys less than %d: %d peaks remain\n", __func__, LOWEST_CUT, peakn);
		//remove deleted 
		/*print_loopt(locopts.a, locopts.n);*/
		for ( i = 0, j = 0; j < locopts.n; ++j) {
			if (!locopts.a[j].del)
				locopts.a[i++] = locopts.a[j];	
		}
		locopts.n = i;
		if (locopts.n >= 3) {
			qsort(locopts.a, locopts.n, sizeof(locopt_t), cmplocopt);	
			
			//more than 3 merge 
			//in case ne		
			//check peak and valleys
			/*print_loopt(locopts.a, 3);*/
			peakn = 0;
			int valley_idx = 0, ltval_idx = 0;
			for ( i = 0; i < 3; ++i) {
				if (locopts.a[i].ispeak) ++peakn;
				else 
					valley_idx = locopts.a[i].idx_s;	
			}
			for ( i = 0; i <3; ++i) {
				if (locopts.a[i].idx_s > valley_idx) ++ltval_idx; 
			}
			fprintf(stderr, "[M::%s] Use top 3 frequent read depth\n", __func__);
			if (peakn == 2 && ltval_idx == 1) {
				fprintf(stderr, "[M::%s] Found a valley in the middle of the peaks, use two-peak mode\n", __func__);
				if (locopts.a[0].idx_s > locopts.a[1].idx_s) swap(int, locopts.a[0].idx_s, locopts.a[1].idx_s);
				if (locopts.a[1].idx_s > locopts.a[2].idx_s) swap(int, locopts.a[1].idx_s, locopts.a[2].idx_s);
				if (locopts.a[0].idx_s > locopts.a[1].idx_s) swap(int, locopts.a[0].idx_s, locopts.a[1].idx_s);
				cutoffs[1] = locopts.a[1].idx_s;
				cutoffs[0] = 2 * locopts.a[0].idx_s - cutoffs[1];
				cutoffs[2] = cutoffs[1] + 1; 
				cutoffs[3] = 2 * locopts.a[2].idx_s - cutoffs[2];
				cutoffs[4] = 3 * locopts.a[2].idx_s;	
				
				free(derts1);
				free(derts1_fit);
				free(derts2);
				kv_destroy(locopts);
				return 0;
			}
		} else {
			fprintf(stderr, "[M::%s] The valley is not in proper position, use one-peak mode\n", __func__);
			locopts.n = 2;	
		} 
	} 
	//check whether this is a haploid or diploid		
					
	int isdip = 1; 
	if (locopts.n) locopts.a[0].idx_s = locopts.a[0].idx_e = max_idx, locopts.a[0].cnt = max_cnt;	
	else {
			locopt_t tmp = (locopt_t){max_idx, 1, max_idx, 0, max_cnt};
			kv_push(locopt_t, locopts, tmp);
	} 
	if (!fhord) {
		int mean = get_mean(depth2cnt);
		fprintf(stderr, "[M::%s] mean: %d, peak: %d, mean %s than peak, treat as %s assembly\n", __func__, mean, max_idx, mean <= max_idx ? "not larger" : "larger", mean <= max_idx ? "haploid":"diploid");
		if (norm_cdf(mean<=max_idx? max_idx : mean, 0.5, mean + max_idx) <=0.95) fprintf(stderr, "[W::%s] mean is not significantly different with peak, please recheck the cutoffs\n", __func__);
		//hump is on right side is hapliod covrage
		if (mean <= max_idx) isdip = 0;
	} else if (fhord == 1) 
		isdip = 0;
	 if (isdip) {
		cutoffs[0] = cutoffs[1] = cutoffs[2] = cutoffs[3] = locopts.a[0].idx_s;
		i = cutoffs[0] >> 2;
		cutoffs[0] = cutoffs[0] - i; //.75 * opt haploid is opt 
		cutoffs[1] = cutoffs[1] + i; //1.25 * opt
		i = cutoffs[2] >> 1;
		cutoffs[2] = cutoffs[2] + i; //1.5opt 
		cutoffs[3] = (cutoffs[3] << 1) + i; //2.5opt 
		cutoffs[4] = 3 * cutoffs[2];		
	} else {
		cutoffs[0] = cutoffs[1] = cutoffs[2] = cutoffs[3] = locopts.a[0].idx_s;
		i = cutoffs[0] >> 3;
		cutoffs[0] = (cutoffs[0] >> 1) - i; //0.75 opt   haploid is opt 	
		cutoffs[1] = (cutoffs[1] >> 1) + i; // 1.25 opt
		i = cutoffs[2] >> 2;
		cutoffs[2] = cutoffs[2] - i;//1.5 * opt
		cutoffs[3] = cutoffs[3] + i;//2.5 * opt
		cutoffs[4] = 3 * cutoffs[2];
	}
	free(derts1);
	free(derts1_fit);
	free(derts2);
	kv_destroy(locopts);
	return 0;
}

int main(int argc, char *argv[])
{
	opt_t opts;
	int c;
	opts.min_mc = 7;
   	opts.min_frac = .1;
	opts.low_cov = opts.dip_cov = opts.upper_cov = -1;	
	opts.stat_fn = 0;
	opts.fhord = 0;
	/*opts.isdip = 0;*/
	char *program;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	while (~(c=getopt(argc, argv, "l:m:u:f:c:d:h"))) {
		switch (c) {
			case 'f': 
				opts.min_frac = atof(optarg);
				break;
			case 'c':
				opts.min_mc = atoi(optarg);
				break;
			case 'l':
				opts.low_cov = atoi(optarg);
				break;
			case 'm':
				opts.dip_cov = atoi(optarg);
				break;
			case 'u':
				opts.upper_cov = atoi(optarg);
				break;
			case 'd':
				opts.fhord = atoi(optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s  [<options>] <STAT> ...\n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -f    FLOAT    minimum depth count fraction to maximum depth count [.1]\n");	
				fprintf(stderr, "         -l    INT      lower bound for read depth\n");	
				fprintf(stderr, "         -m    INT      transition between haploid and diploid\n");	
				fprintf(stderr, "         -u    INT      upper bound for read depth\n");	
				/*fprintf(stderr, "         -c    INT      minimum spanning hump [7]\n");*/
				fprintf(stderr, "         -d             treat as haploid assembly or diploid assembly, 1: haploid, others: diploid [0]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 1 <= argc) 
		opts.stat_fn = argv[optind];	
	 
	 
	//read 
	if (~opts.dip_cov) {
		if (~opts.upper_cov) {
			fprintf(stderr, "Set up cutoffs manually\n");	
			fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\n", ~opts.low_cov ? opts.low_cov : LOWEST_CUT, opts.dip_cov - 1, opts.dip_cov - 1, opts.dip_cov, opts.dip_cov, opts.upper_cov);	
			return 0;	
		} else {
			fprintf(stderr, "[E::%s] please set high read depth cutoff\n", __func__);
			goto help;
		}
	}
	if (!opts.stat_fn) {
		fprintf(stderr, "[E::%s] require a read depth statistics file\n", __func__);
		goto help;
	}
	int cutoffs[5];
	uint32_t *depth2cnt = read_counts(opts.stat_fn);	
	//process 
	if (!calcuts(depth2cnt, cutoffs, opts.min_mc, opts.min_frac, opts.fhord)) 
   		//output 
		fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\n", LOWEST_CUT, cutoffs[0], cutoffs[1], cutoffs[2], cutoffs[3], cutoffs[4]);	
	return 0;
}

