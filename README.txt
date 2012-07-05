Rainbow v2.0

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
	rainbow div -i rbcluster.out -o rbdiv.out
	rainbow merge -o rbasm.out -a -i rbdiv.out -N500

----------------------------------------------------------------------------------
rainbow 1.1 -- <ruanjue@gmail.com, chongzechen@gmail.com>
Usage: rainbow <cmd> [options]
Input  File Format: paired fasta/fastq file(s)
Output File Format: <seqid:int>\t<cluster_id:int>\t<read1:string>\t<read2:string>[\t<pre_cluster_id:int>]

 cluster
  -1 <string> Input fasta/fastq file, supports multiple '-1'
  -2 <string> Input fasta/fastq file, supports multiple '-2' [null]
  -l <int>    Read length, default: 0 variable
  -m <int>    Maximum mismatches [4]
  -e <int>    Exactly matching threshold [2000]
  -L          Low level of polymorphism
 div
  -i <string> Input file [stdin]
  -o <string> Output file [stdout]
  -k <int>    K_allele, min variants to create a new group [2]
  -K <int>    K_allele, divide regardless of frequency when num of variants exceed this value [50]
  -f <float>  Frequency, min variant frequency to create a new group [0.2]
 merge
  -i <string> Input rbasm output file [stdin]
  -a          output assembly or not [not]
  -o <string> Output file for merged contigs, one line per cluster [stdout]
  -N <int>    Maximum number of divided clusters to merge [300]

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
