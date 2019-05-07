# Purge_Dups 

purge haplotigs and overlaps in an assembly based on read depth

## Overview

purge\_dups pipeline calculates read depth of given read files, even though the pipeline still works without sequencing data, it is high recommended to provide it. Based on the read depth, purge\_dups pipeline will resolve the coverage cutoffs for haploid and diploid sequence. Meanwhile, the pipeline also splits the draft assembly into numbers of sequences, does a self-to-self alignment, chains the segmented consistent alignments. 

![purge_dups pipeline](https://github.com/dfguan/purge_dups/blob/master/purge_dupspipeline.png)

## Directory Structure

- scripts/pd\_config.py: script to generate a configuration file used by run\_purge\_dups.py.
- scripts/run\_purge\_dups.py: script to run the purge\_dups pipeline. 
- scripts/run\_busco:  script to run busco, dependency: busco.
- scripts/run\_kcm:  script to make k-mer comparison plot. 
- scripts/sub.sh: shell script to submit a farm job.
- src: purge_dups source files
 


## Installation
Run the following commands to install runner (required):

```
git clone https://github.com/dfguan/runner.git
cd runner && python3 setup.py install --user
```
Run the following commands to intall purge_dups (required):

```
https://github.com/dfguan/purge_dups.git
cd purge_dups/src && make

```

If you also want to try k-mer comparision plot, run the following commands to install the tool (optional). 

```
git clone https://github.com/dfguan/KMC.git 
cd KMC && make -j 16
```
## Usage

### Step 1. Use pd\_config.py to generate a configuration file. 

```
pd_config.py <ref> <ref_dir> <pbdb_dir> <10xdb_dir> <local_dbdir> 
   This will generate a configuration file named with config.ASSEMBLY_PREFIX.json
   
   <ref>         the location of assembly file, can be fasta or fasta.gz
   <ref_dir>     the directory of copy of the assembly file, no need to exist
   <pbdb_dir>    Pacbio data directory 
   <10xdb_dir>   10x data directory
   <local_dbdir> diretory to keep Pacbio and 10x file lists. 

```

Example:

```
./scripts/pd_config.py ~/vgp/release/insects/iHelSar1/iHelSar1.PB.asm1/iHelSar1.PB.asm1.fa.gz refs ~/vgp/build/insects/iHelSar1/PacBio/ ~/vgp/build/insects/iHelSar1/10X/ iHelSar1
```

### Step 2. Modify the configuration file manually (optional). 

configuration file is in json format, it has all the information required by run\_purge\_dups.py. Here is an example of a configuration file. 

```
{
  "kcp": {   
  	"core": 12,
  	"skip": 0, 
  	"prefix": "iHelSar1.PB.asm1_purged_kcm", 
  	"fofn": "iHelSar1ls/10x.fofn", 
  	"mem": 30000,
  	"tmpdir": "kcp_tmp” 
  	},
  "gs": { "mem": 10000},
  "out_dir": "iHelSar1.PB.asm1", 
  "pd": { 
  	"mem": 20000, 
  	"queue": "normal”
  	},
  "cc": { 
  	"core": 12, 
  	"queue": "normal",  
  	"isdip": 1, 
  	"skip": 0, 
  	"fofn": "iHelSar1ls/pb.fofn", 
  	"mem": 20000, 
  	"ispb": 1 
  }, 
  “sa”: { 
  	“core”: 12, 
  	“mem”: 10000, 
  	“queue”: “normal”
  	},
  “busco”: { 
  	“core”: 12, 
  	“lineage”: “insecta”, 
  	“prefix”: “iHelSar1.PB.asm1_purged”, 
  	“queue”: “long”,  
  	“skip”: 0, 
  	“mem”: 20000,
  	”tmpdir”: “busco_tmp”
  	},
  “ref”:”/lustre/scratch116/vr/projects/vgp/user/dg30/dg30/projects/vgp/purge_dups/190418.primary/purge_dups/refs/iHelSar1.PB.asm1.fa”
}
```

This file use several key words to define resource allocation, input files or output files, they are listed as follows.   

- **core**: CPU number
- **skip**: Bool value set to skip this job
- **prefix**: Output file prefix
- **fofn**: Sequencing files list
- **mem**: Maximum amount of RAM in MB
- **tmpdir**: Temporary directory
- **lineage**: Busco database 
- **queue**: job queue
- **ref**: assembly file path
- **out_dir**: working directory
- **ispb**: Bool value set for pacbio data, 0 for Illumina data

**Notice**: **isdip** is deprecated. 

The dictionary "kcp" keeps paramaters for run_kcm script.  
The dictionary "gs" sets parameters for get\_seqs (purge\_dups executable file), designed to bin primary contigs and haplotigs.  
The dictionary "pd" sets parameters for purge\_dups (purge\_dups executable file), designed to purge haplotigs and overlaps in an assembly.  
The dictionary **"cc"** sets parameters for **minimap2/bwa**.  
The dictionary **sa** sets parameters for minimap2.  
The dictionary busco sets parameters for run\_busco. 

### Step 3. Use run\_purge\_dups.py to run the pipeline

```
usage: run_purge_dups.py [-h] [-p PLTFM] [-w WAIT] [-r RETRIES] [--version]
                         config bin_dir spid

purge_dups wrapper

positional arguments:
  config                configuration file
  bin_dir               directory of purge_dups executable files
  spid                  species identifier

optional arguments:
  -h, --help            show this help message and exit
  -p PLTFM, --platform PLTFM
                        workload management platform
  -w WAIT, --wait WAIT  <int> seconds sleep intervals
  -r RETRIES, --retries RETRIES
                        maximum number of retries
  --version             show program's version number and exit
```

Example: 

```
python scripts/run_purge_dups.py config.iHelSar1.json src iHelSar1
```

### Other Modification

If the busco and k-mer comparison plot scripts are working, please modify them with the following instructions. 

- run\_busco: set the PATH variables in run_busco to your own path. 
- run\_kcm: set kcm_dir to your own KMC directory.


## Results

After the pipeline is finished, there will be four new directories in the working directory (set in the configuration file).  

- **coverage**: coverage cutoffs, coverage histogram and base-level coverage files
- **split_aln**: segmented assembly file and a self-alignment paf file. 
- **purge_dups**: duplicate sequence list. 
- **seqs**: purged primary contigs ending with .purge.fa and haplotigs ending with .red.fa, also K-mer comparison plot and busco results are also in this directory.  

## Contact

Wellcome to use, you can use github webpage to report an issue or email me dfguan9@gmail.com with any advice. 









