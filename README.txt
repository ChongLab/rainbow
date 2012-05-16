Rainbow v1.0.1

Description
===========
Rainbow package consists of several programs used for RAD-seq related 
clustering and de novo assembly.

Installation
============
Type 'make' to compile Rainbow package. You can copy the executables/scripts
to your specific location (e.g. a directory in your $PATH). Or you can set
the PATH environment that leads to this directory.


Usage of Rainbow package
========================
EXAMPLE: a typical use of Rainbow step by step

	rainbow cluster -1 1.fq  -2 2.fq > rbcluster.out 2> log
	rainbow div -i rbcluster.out > rbdiv.out
	rbasm -i rbdiv.out -o rbasm.out
	rainbow merge -a rbasm.out -v rbdiv.out -p 0.002 > merged.txt ## [-p] is the estimated heterogosity
	rerun_rbasm.pl merged.txt rbdiv.out rbasm.out # the output file 'final_asm.fa' contains the final contigs

----------------------------------------------------------------------------------
rainbow 1.1 -- <ruanjue@gmail.com, chongzechen@gmail.com>
Usage: rainbow <cmd> [options]
Input  File Format: paired fasta/fastq file(s)
Output File Format: <seqid:int>\t<cluster_id:int>\t<read1:string>\t<read2:string>[\t<pre_cluster_id:int>]

 cluster
  -1 <string> Input fasta/fastq file, supports multiple '-1'
  -2 <string> Input fasta/fasta file, supports multiple '-2' [null]
  -l <int>    Read length, default: 0 variable
  -m <int>    Maximum mismatches [2]
  -e <int>    Exactly matching threshold [2000]
 div
  -i <string> Input file [stdin]
  -k <int>    K_allele, min variants to create a new group [2]
  -K <int>    K_allele, divide regardless of frequency when num of variants exceed this value [50]
  -f <float>  Frequency, min variant frequency to create a new group [0.2]
 merge  ***this procedure should be run after rbasm***
  -a <string> Input rbasm output file [stdin]
  -v <string> Input rainbow divided file [stdin]
  -p <float>  maximum heterozygosity to collapse, should be specifed according to the estimated
              polymorphism of the species [0.01]
  -l <int>    Minimum overlap to collapse two contigs [100]
  -n <int>    Maximum number of contigs to execute pairwise alignment [50]

----------------------------------------------------------------------------------
rbasm: a greedy assembler to locally assemble each cluster produced by rainbow
Local assemble fragments around restriction sites
Usage: rbasm [options]
 -i <string> Input file [STDIN]
 -o <string> Output file [STDOUT]
 -l <int>    Minium length of overlap [5]
 -s <float>  Minium similiarity of overlap [0.90]
 -r <int>    Minium reads to execute assembly [5]
 -R <int>    Maxium reads to execute assembly [200]

----------------------------------------------------------------------------------
<Obsoleted> rbmergetag: a program merges divided results to evaluate clustering performance.
                  Users should omit this program when de novo assembling RAD-seq reads.
Usage: rbmergetag [options]
Options:
 -i <string>    Input file name [stdin]
 -o <string>    Output file name [stdout]
 -j <cns|merge> Job type, cns: consensus, merge: merging, [merge]
 -m <int>       Maximum mismatches to merge two groups [1]
 -h             Show this document

----------------------------------------------------------------------------------
