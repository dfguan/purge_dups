#!/usr/bin/env python3
#generate a config.*.json used for run_purge_dups.py, just work for our directory structure 
import sys, os,json
import argparse


def fl_exist(s):
    return os.path.isfile(s)

def getfn(p):
    return os.path.basename(p)

def get_abspath(p):
    return os.path.abspath(p)

def alz_fn(d, fn):
    fl_path = "{0}/{1}".format(d, fn)
    fn = os.path.splitext(fn)[0]
    fn_list = fn.split('_')
    enzyme = fn_list[2]
    tech = fn_list[1]
    return [fl_path, enzyme, tech]

def gen_config(r, d, fn, skipB):
    #write a new config file
    ref_fn = os.path.splitext(os.path.basename(r))[0] 
    busco_lineage = {'a': "tetrapoda", 'b': "aves", 'd': "embryophyta", 'e': "metazoa", 'f': "actinopterygii", 'i': "insecta", 'h': "eukaryota", 'm': "mammalia", 'q': "arthropoda", 's': "vertebrata", 'x':"metazoa"}
    used_lineage = busco_lineage[ref_fn[0]] if ref_fn[0] in busco_lineage else ""
    # out_fn = "config.{}.json".format(ref_fn) 
    out_fn = fn 
    f = open(out_fn, 'w')
    jd = {
            "cc":{"fofn":"", "isdip":1, "core":12,"mem":20000, "queue":"normal", "mnmp_opt":"", "bwa_opt":"", "ispb":1, "skip":0},
            "sa":{"core":12, "mem":10000, "queue":"normal"},
            "busco":{"core":12, "mem":20000, "queue":"long", "skip":0, "lineage":used_lineage, "prefix":ref_fn+"_purged", "tmpdir":"busco_tmp"},
            "pd":{"mem": 20000, "queue": "normal"}, 
            "gs": {"mem": 10000, "oe": 1}, 
            "kcp": {"core":12, "mem":30000, "fofn":"", "prefix":ref_fn+"_purged_kcm", "tmpdir":"kcp_tmp", "skip": 0}, 
            "ref":r, "out_dir":ref_fn
    }
    
    fofn_fn = "{}/pb.fofn".format(d)
    if fl_exist(fofn_fn):
        jd["cc"]["fofn"] = fofn_fn
    else:
        jd["cc"]["skip"] = 1
    fofn_fn = "{}/10x.fofn".format(d)
    if fl_exist(fofn_fn):
        jd["kcp"]["fofn"] = fofn_fn
    else:
        jd["kcp"]["skip"] = 1
    if skipB: 
        jd["busco"]["skip"] = 1

    json.dump(jd, f, indent = 2) 
    f.close()
def worker(ref, ref_dir, pbfofn, txfofn, ldbdir, out_fn, skipB):
    if not os.path.isfile(ref):
        print ("file {} doesn't exist".format(ref))
        return 1
    if not os.path.isdir(ref_dir):
        os.mkdir(ref_dir)
    cp_cmd = "cp {0} {1}".format(ref, ref_dir)
    rpath = "{0}/{1}".format(ref_dir, getfn(ref))
    if ref[-3:] == ".gz":
        cp_cmd = "zcat {0} > {1}/{2}".format(ref, ref_dir, getfn(ref)[:-3])
        rpath = "{0}/{1}".format(ref_dir, getfn(ref)[:-3])
    os.system(cp_cmd)
    
    if not os.path.isfile(pbfofn):
        print ("file {} doesn't exist".format(pbfofn))
        return 1
    
    if not os.path.isdir(ldbdir):
        os.mkdir(ldbdir)
    cp_cmd = "cp {0} {1}/pb.fofn".format(pbfofn, ldbdir)
    os.system(cp_cmd)

    if not os.path.isdir(ldbdir):
        os.mkdir(ldbdir)
    if txfofn:
        cp_cmd = "cp {0} {1}/10x.fofn".format(txfofn, ldbdir)
        os.system(cp_cmd)
    gen_config(get_abspath(rpath), ldbdir, out_fn, skipB)  
    return 0

def proc_ref(ref, ref_dir, pbdbdir,txdbdir, ldbdir):
    #copy ref to refdir
    if not os.path.isfile(ref):
        print ("file {} doesn't exist".format(ref))
        return 1
    if not os.path.isdir(ref_dir):
        os.mkdir(ref_dir)
    cp_cmd = "cp {0} {1}".format(ref, ref_dir)
    rpath = "{0}/{1}".format(ref_dir, getfn(ref))
    if ref[-3:] == ".gz":
        cp_cmd = "zcat {0} > {1}/{2}".format(ref, ref_dir, getfn(ref)[:-3])
        rpath = "{0}/{1}".format(ref_dir, getfn(ref)[:-3])
    os.system(cp_cmd)
    # if os.WEXITSTATUS() \ci 
    
    d = pbdbdir
    locd = ldbdir
    if not os.path.isdir(locd):
        os.mkdir(locd)

    if os.path.isdir(d):
        d = get_abspath(d) 
        find_cmd = 'find {0}/fasta -maxdepth 1 -name "*.fasta.gz" > {1}/pb.fofn'.format(d, locd)
        os.system(find_cmd)
    else:
        print ("directory {} doesn't exist".format(d))
        return  1

    d =txdbdir 
    if os.path.isdir(d):
        d = get_abspath(d) 
        find_cmd = 'ls {0}/*R*.fastq.gz | awk \'{{ if (NR%2==1) print $1\"\\t\"23; else print $1\"\\t\"0; }}\' > {1}/10x.fofn'.format(d, locd)
        # print (find_cmd)
        os.system(find_cmd)
    # d = './' + os.path.basename(r).split(".")[0] #directory
        gen_config(get_abspath(rpath), locd)  
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='generate a configuration file in json format')
    parser.add_argument('-s', '--srfofn', type=str, action="store", dest = "srf", help ='list of short reads files (one record per line, the record is a tab splitted line of abosulte file path plus trimmed bases, refer to https://github.com/dfguan/KMC) [NONE]')
    parser.add_argument('-l', '--localdir', type=str, action="store", dest = "locd", help ='local directory to keep the reference and lists of the pacbio, short reads files [.]', default=".")
    parser.add_argument('-B', '--skipB',  action = "store_true", dest = "skipB", help = 'skip running busco [False]')
    parser.add_argument('-n', '--name', type=str, action="store", dest = "fn", help ='output config file name [config.json]', default = "config.json")
    parser.add_argument('--version', action='version', version='%(prog)s 0.0.0')
    parser.add_argument('ref', type=str, action="store", help = "reference file")
    parser.add_argument('pbfofn', type=str, action="store", help = "list of pacbio file (one absolute file path per line)")
    opts = parser.parse_args()
    sys.exit(worker(opts.ref, opts.locd, opts.pbfofn,opts.srf, opts.locd, opts.fn, opts.skipB))
