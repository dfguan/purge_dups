#!/usr/bin/env python3

from runner.manager import manager
from runner.hpc import hpc
from multiprocessing import Process, Pool
import sys, os, json
import argparse

# utilities
def checkf(p):
    return os.path.isfile(p)

def checkd(p):
    return os.path.isdir(p)

def mkdir(d):
    if not checkd(d):
        os.makedirs(d)

def getfn(p):
    return os.path.basename(p)
def getd(p):
    dirn = os.path.dirname(p)
    return dirn if dirn else "."
# def touch(p):
    # open(p, "w").close()
def get_rm_prefix(p): # a.b.c.d return a.b.c
    return ".".join(getfn(p).split(".")[0:-1])

def get_lm_prefix(p): # a.b.c.d return a
    return getfn(p).split(".")[0]

# INPUT: assembly, output_dir, spid 
     
def split_aln(man, pltfm, ref, core_lim, mem_lim, queue, out_dir, bin_dir, spid):
    mkdir(out_dir)
    
    out_fn = "{0}/{1}.split.fa".format(out_dir, get_rm_prefix(ref))
    jcmd = "{0}/split_fa {1} > {2}".format(bin_dir, ref, out_fn)
    jjn = "split_{}".format(spid)
    jout = "{0}/{1}_%J.o".format(out_dir, jjn)
    jerr = "{0}/{1}_%J.e".format(out_dir, jjn)

    j = hpc(pltfm, cmd=jcmd, core=1, mem = 5000, jn=jjn, out=jout, err=jerr)
    rtn = man.start([j])
    if not rtn:
        in_fn = out_fn
        out_fn = "{0}/{1}.split.paf".format(out_dir, get_rm_prefix(ref))
        jcmd = "minimap2 -xasm5 -DP {0} {0} > {1}".format(in_fn, out_fn)
        jjn = "self_aln_{}".format(spid)
        jout = "{0}/{1}_%J.o".format(out_dir, jjn)
        jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
        j = hpc(pltfm, cmd=jcmd, core=core_lim, mem = mem_lim, jn=jjn, out=jout, err=jerr)
        rtn = man.start([j])
    exit(rtn)

# INPUT: input_dir, output_dir, spid 
def purge_dups(man, pltfm, paf_fn, base_cov_fn, cutoff_fn, core_lim, mem_lim, queue, out_dir, bin_dir, spid):
    mkdir(out_dir)
    fn_param = ""
    if base_cov_fn != "" and cutoff_fn != "":
        fn_param = "-c {0} -T {1}".format(base_cov_fn, cutoff_fn)
    out_fn = "{}/dups.bed".format(out_dir) 
    # jcmd = "{0}/purge_dups -1 {1} {2} > {3}".format(bin_dir, fn_param, paf_fn, out_fn)
    jcmd = "{0}/purge_dups -2 {1} {2} > {3}".format(bin_dir, fn_param, paf_fn, out_fn)
    jjn = "purge_dups_{}".format(spid)
    jout = "{0}/{1}_%J.o".format(out_dir, jjn)
    jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
    j = hpc(pltfm, cmd=jcmd, core=core_lim, mem = mem_lim, jn=jjn, out=jout, err=jerr)
    rtn = man.start([j], True)
    return rtn    
#INPUT: dups assembly
def get_seqs(man, pltfm, ref, dups_fn, core_lim, mem_lim, onlyend, out_dir, bin_dir, spid):
    mkdir(out_dir) 
    out_fn = "{0}/{1}.purged.fa".format(out_dir, get_rm_prefix(ref)) 
    out_red_fn = "{0}/{1}.red.fa".format(out_dir, get_rm_prefix(ref))
    out_prefx="{0}/{1}".format(out_dir, get_rm_prefix(ref))
    jcmd = "{0}/get_seqs -p {5} {1} {2}".format(bin_dir, dups_fn, ref, out_fn, out_red_fn, out_prefx)
    if onlyend:
        jcmd = "{0}/get_seqs -e -p {5} {1} {2}".format(bin_dir, dups_fn, ref, out_fn, out_red_fn, out_prefx)
    jjn = "get_seqs_{}".format(spid)
    jout = "{0}/{1}_%J.o".format(out_dir, jjn)
    jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
    j = hpc(pltfm, cmd=jcmd, core=core_lim, mem = mem_lim, jn=jjn, out=jout, err=jerr)
    rtn = man.start([j], True)
    return rtn 

def bwa_index(p):
    [man, pltfm, ref, out_dir, spid] = p
    jcmd = "bwa index {}".format(ref)
    mem_lim = 20000
    jjn = "bwa_index_{}".format(spid)
    jout = "{0}/{1}_%J.o".format(out_dir, jjn)
    jerr = "{0}/{1}_%J.e".format(out_dir, jjn)

    j = hpc(pltfm, cmd=jcmd, core=1, mem = mem_lim, jn=jjn, out=jout, err=jerr)


    return man.start([j])
    # print ("code {}".format(code))
    # rtn.append(man.start(jobs))
    # rtn.append(code)



#INPUT: pacbio.fofn/illumina.fofn, assembly
def cal_cov(man, pltfm, ref, ispb, isdip, fofn, core_lim, mem_lim, queue, mnmp_opt, bwa_opt, skip, out_dir, bin_dir, spid, ispurged):
    mkdir(out_dir)
    if skip == 1:
        exit(1)
    if not ispb:
        # index     
        if bwa_index([man, pltfm, ref, out_dir, spid]):
            exit(1)

    jobs = []
    out_fns = []
    with open(fofn, "r") as f:
        for fl in f:
            fl_strip = fl.strip()
            if ispb:
                fn_prefix = get_rm_prefix(getfn(fl_strip))
                out_fn = "{0}/{1}.paf".format(out_dir, fn_prefix)
                out_fns.append(out_fn)
                idx_opt = "-I {}".format("4G" if os.path.getsize(ref) < 4e9 else "10G")
                if mnmp_opt != "":
                    jcmd = "minimap2 {4} {5} -t {0} {1} {2} >{3}".format(core_lim, ref, fl_strip, out_fn, idx_opt, mnmp_opt)
                else:
                    jcmd = "minimap2 {4} -x map-pb -t {0} {1} {2} >{3}".format(core_lim, ref, fl_strip, out_fn, idx_opt)
                jjn = "minimap_{}".format(fn_prefix)
                jout = "{0}/{1}_%J.o".format(out_dir, jjn)
                jerr = "{0}/{1}_%J.e".format(out_dir, jjn)
                j = hpc(pltfm, cmd=jcmd, core=core_lim, mem = mem_lim, jn=jjn, out=jout, err=jerr)
                jobs.append(j)
            else:
                [r1, r2] = fl_strip.split('\t')
                fn_prefix = get_rm_prefix(getfn(r1))
                out_fn = "{0}/{1}.bam".format(out_dir, fn_prefix)
                out_fns.append(out_fn)
                if bwa_opt != "":
                    jcmd = "bwa mem -t {0} {5} {1} {2} {3} | samtools view -b - >{4}".format(core_lim, ref, r1, r2, out_fn, bwa_opt)
                else:
                    jcmd = "bwa mem -t {0} {1} {2} {3} | samtools view -b - >{4}".format(core_lim, ref, r1, r2, out_fn)
                jjn = "bwa_mem_{}".format(fn_prefix)
                jout = "{0}/{1}_%J.o".format(out_dir, jjn)
                jerr = "{0}/{1}_%J.e".format(out_dir, jjn)

                j = hpc(pltfm, cmd=jcmd, core=core_lim + 1, mem = mem_lim, queue=queue, jn=jjn, out=jout, err=jerr)
                jobs.append(j)
        f.close()
    rtn = man.start(jobs)

    jobs = []
    if not rtn:
        if ispb:
            jcmd = "{0}/pbcstat -O {2} {1}".format(bin_dir, " ".join(out_fns), out_dir)
            jjn = "pbcstat_{}".format(spid)
            jout = "{0}/{1}.o".format(out_dir, jjn)
            jerr = "{0}/{1}.e".format(out_dir, jjn)

            j = hpc(pltfm, cmd=jcmd, core=1, mem = 10000, jn=jjn, out=jout, err=jerr)
            jobs.append(j)
        else:
            jcmd = "{0}/ngscstat -O {2} {1}".format(bin_dir, " ".join(out_fns), out_dir)
            jjn = "ngscstat_{}".format(spid)
            jout = "{0}/{1}.o".format(out_dir, jjn)
            jerr = "{0}/{1}.e".format(out_dir, jjn)

            j = hpc(pltfm, cmd=jcmd, core=1, mem = 30000, jn=jjn, out=jout, err=jerr)
            jobs.append(j)
        rtn = man.start(jobs, True)
    # if not rtn:
        # in_fn = "{}/PB.stat".format(out_dir)
        # out_prefix = "{0}.{1}".format(spid, "purged" if ispurged == 1 else "origin")
        # jcmd = "Rscript ~/plot_depthgraph.R {} {}".format(in_fn, out_prefix)
        # jjn = "depthplot_{}".format(spid)
        # jout = "{0}/{1}.o".format(out_dir, jjn)
        # jerr = "{0}/{1}.e".format(out_dir, jjn)
        # j = hpc(pltfm, cmd=jcmd, core=1, mem = 500, jn=jjn, out=jout, err=jerr)
        # man.start([j])

    if not rtn:
        in_fn = "{}/PB.stat".format(out_dir)
        out_fn = "{}/cutoffs".format(out_dir)
        jcmd = "{0}/calcuts {3} {1} > {2}".format(bin_dir, in_fn, out_fn, "")
        jjn = "calcuts_{}".format(spid)
        jout = "{0}/{1}.o".format(out_dir, jjn)
        jerr = "{0}/{1}.e".format(out_dir, jjn)

        j = hpc(pltfm, cmd=jcmd, core=1, mem = 2000, jn=jjn, out=jout, err=jerr)
        rtn = man.start([j], True)
    exit(rtn)

def run_kcm(man, pltfm, skip, spid, fasta, mem, core, kmer, reads, prefix, tmpdir):
    if skip == 1:
        return 0
    else:
        jcmd = "run_kcm {0} {1} {2} {3} {4} {5} {6} {7}".format(spid, mem, core, kmer, reads, fasta, prefix, tmpdir)
        jjn = "kcm_{}".format(spid)
        jout = "{}.o".format(jjn)
        jerr = "{}.e".format(jjn)
        j = hpc(pltfm, cmd=jcmd, core=core, mem = mem, queue = "normal", jn=jjn, out=jout, err=jerr)
        rtn = man.start([j])
        return rtn 

def run_busco(man, pltfm, skip, workdir, spid, fasta, mem, core, queue, prefix, lineage, tmpdir):
    if skip == 1:
        return 0
    else:
        os.chdir(workdir)
        jcmd = "run_busco2 {0} {1} {2} {3} {4}".format(fasta, core, prefix, lineage, tmpdir)
        jjn = "busco_{}".format(spid)
        jout = "{}.o".format(jjn)
        jerr = "{}.e".format(jjn)
        j = hpc(pltfm, cmd=jcmd, core=core, mem = mem, queue = queue, jn=jjn, out=jout, err=jerr)
        rtn = man.start([j], True)
        return rtn 

def cont(config_fn, bin_dir, spid, pltfm, _wait, _retries):

    f = open(config_fn, "r")
    config_dict = json.load(f)

    out_dir = config_dict["out_dir"]
    ref = config_dict["ref"]
    
    ref_pref = get_rm_prefix(ref)

    mkdir(out_dir)


    man = manager(wait=_wait, retries=_retries)
    
    procs = []
    cur_d = config_dict["cc"]
    out_cov_dir = "{}/coverage".format(out_dir)
    p = Process(target=cal_cov, args=(man, pltfm, ref, cur_d["ispb"], cur_d["isdip"], cur_d["fofn"], cur_d["core"], cur_d["mem"], cur_d["queue"], cur_d["mnmp_opt"], cur_d["bwa_opt"], cur_d["skip"], out_cov_dir, bin_dir, spid, 0))
    procs.append(p)
    cur_d = config_dict["sa"]
    out_sa_dir = "{}/split_aln".format(out_dir)
    p = Process(target=split_aln, args=(man, pltfm, ref, cur_d["core"], cur_d["mem"], cur_d["queue"], out_sa_dir, bin_dir, spid))
    procs.append(p)
    print ("calculate coverage and self-alignment") 
    for p in procs:
        p.start()
    for p in procs:
        p.join()
    # purge dups
    in_paf_fn = "{0}/{1}.split.paf".format(out_sa_dir, ref_pref)
    in_cov_fn = "{0}/PB.base.cov".format(out_cov_dir)
    in_cutoffs = "{0}/cutoffs".format(out_cov_dir)
    out_pd_dir = "{0}/purge_dups".format(out_dir)
    
    rtn = 1
    print ("purge duplicates") 
    pd_mem = 20000
    pd_queue = "normal"
    if "pd" in config_dict:
        pd_mem = config_dict["pd"]["mem"]
        pd_queue = config_dict["pd"]["queue"]

    if not procs[1].exitcode:
        if not procs[0].exitcode:
            rtn = purge_dups(man, pltfm, in_paf_fn, in_cov_fn, in_cutoffs, 1, pd_mem, pd_queue, out_pd_dir, bin_dir, spid)
        else:
            rtn = purge_dups(man, pltfm, in_paf_fn, "", "", 1, pd_mem, pd_queue, out_pd_dir, bin_dir, spid)
    if not rtn:
        gs_mem = 10000
        gs_onlyend = 1
        if "gs" in config_dict:
            gs_mem = config_dict["gs"]["mem"]
        if "oe" in config_dict:
            gs_onlyend=config_dict["gs"]["oe"]
        in_dups_fn = "{0}/dups.bed".format(out_pd_dir)
        out_dir = "{}/seqs".format(out_dir)
        rtn = get_seqs(man, pltfm, ref, in_dups_fn, 1, gs_mem, gs_onlyend, out_dir, bin_dir, spid) 
    
    procs = [] 
    workdir = out_dir
    if "busco" in config_dict and not rtn:
        cur_d = config_dict["busco"]
        fasta = "{}.purged.fa".format(get_rm_prefix(ref))
        
        p = Process(target=run_busco, args=(man, pltfm, cur_d["skip"], workdir, spid, fasta, cur_d["mem"], cur_d["core"], cur_d["queue"], cur_d["prefix"], cur_d["lineage"], cur_d["tmpdir"]))
        procs.append(p)
    
    if not rtn:
        purged_ref = "{0}/{1}".format(workdir, fasta)
        out_cov_dir = "{0}/cal_cov".format(out_dir)
        cur_d = config_dict["cc"]
        
        if cur_d["ispb"] == 1:
            p = Process(target=cal_cov, args=(man, pltfm, purged_ref, cur_d["ispb"], cur_d["isdip"], cur_d["fofn"], cur_d["core"], cur_d["mem"], cur_d["queue"], cur_d["mnmp_opt"], cur_d["bwa_opt"], cur_d["skip"], out_cov_dir, bin_dir, spid, 1))
            procs.append(p)

    if not rtn and "kcp" in config_dict:
        purged_ref = "{0}/{1}".format(workdir, fasta)
        cur_d = config_dict["kcp"]
        p = Process(target=run_kcm, args=(man, pltfm, cur_d["skip"], spid, purged_ref, cur_d["mem"], cur_d["core"], 21, cur_d["fofn"], cur_d["prefix"], cur_d["tmpdir"]))
        procs.append(p)

    for p in procs:
        p.start()
    for p in procs:
        p.join()
    return rtn

if __name__ == "__main__":
# this part should be implemented in superclass ? don't know how to do it.  
    parser = argparse.ArgumentParser(description='purge_dups wrapper')
    parser.add_argument('-p', '--platform', type=str, action="store", dest = "pltfm", help ='workload management platform, input bash if you want to run locally', default='lsf')
    parser.add_argument('-w', '--wait', type=int, action = "store", dest = "wait", help = '<int> seconds sleep intervals', default = 10)
    parser.add_argument('-r', '--retries', type = int, action = "store", dest = "retries", help = 'maximum number of retries', default= 2)
    parser.add_argument('--version', action='version', version='v 0.0.3')
    parser.add_argument('config', type=str, action="store", help = "configuration file")
    parser.add_argument('bin_dir', type=str, action="store", help = "directory of purge_dups executable files")
    parser.add_argument('spid', type=str, action="store", help = "species identifier")
    opts = parser.parse_args()
    sys.exit(cont(opts.config, opts.bin_dir, opts.spid, opts.pltfm, opts.wait, opts.retries))
