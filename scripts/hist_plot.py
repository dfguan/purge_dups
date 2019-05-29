#!/usr/bin/env python3
# read depth histogram plot

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, argparse


def col_hist(stat_fn, delim):
    hists = []
    #we consider the coverage histogram start with 0
    with open(stat_fn) as f:
        for ln in f:
            lnlist =ln.strip().split(delim)
            hists.append(int(lnlist[1])) 
        f.close()
    return hists

def mk_plot(hists, ttle, xm, xM, ym, yM, out_fl):
    
    if ttle is None:
        ttle = "read depth histogram"
    if xm is None:
        xm = 0
    if xM is None:
        xM = len(hists) - 2 # ignore the last read depth count
    if ym is None:
        ym = 0
    if yM is None:
        yM = 1.2 * max(hists)
    x = [t for t in range(xm, xM)]
    plt.axis([xm, xM, ym, yM])
    width = 8
    height = 6
    plt.figure(num=None, figsize=(width, height))
    plt.plot(x, hists[xm:xM], label = "l", color="blue") # 
    plt.xticks([z for z in range(xm, xM, 10)], fontsize=3)
    plt.title(ttle)
    plt.gca().xaxis.grid(True, color="black", alpha=0.2)
    # plt.grid(True, color="black", alpha=0.2)
    # plt.gca().get_legend().remove()

    plt.tight_layout() 
    plt.savefig(out_fl, dpi = 300) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='read depth histogram plot')

    parser.add_argument('-y', '--ymin', type=int, action="store", dest = "ymin", help ='set ymin')
    parser.add_argument('-x', '--xmin', type=int, action = "store", dest = "xmin", help = 'set xmin')
    parser.add_argument('-Y', '--ymax', type=int, action="store", dest = "ymax", help ='set ymax')
    parser.add_argument('-X', '--xmax', type=int, action = "store", dest = "xmax", help = 'set xmax')
    parser.add_argument('-t', '--title', type = str, action = "store", dest = "title", help = 'figure title [NULL]', default="")
    parser.add_argument('-d', '--delim', type = str, action = "store", dest = "delim", help = 'delimiter', default="\t")
    parser.add_argument('-v', '--version', action='version', version='histplot 0.0.0')
    parser.add_argument('stat_fn', type=str, action="store", help = "stat file")
    parser.add_argument('out_fn', type=str, action="store", help = "output file")
    opts = parser.parse_args()
    hists = col_hist(opts.stat_fn, opts.delim)
    mk_plot(hists, opts.title, opts.xmin, opts.xmax, opts.ymin, opts.ymax, opts.out_fn) 

