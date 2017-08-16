### Rainbow v2.0.4

Author: Zechen Chong & Jue Ruan
Email: chongzechen@gmail.com or ruanjue@gmail.com


## Description
Rainbow package consists of several programs used for RAD-seq related 
clustering and de novo assembly.

## Installation
Type 'make' to compile Rainbow package. You can copy the executables/scripts
to your specific location (e.g. a directory in your $PATH). Or you can set
the PATH environment that leads to this directory.


## Usage of Rainbow package

***EXAMPLE***: a typical use of Rainbow step by step
```
	rainbow cluster -1 1.fq  -2 2.fq > rbcluster.out 2> log
	rainbow div -i rbcluster.out -o rbdiv.out
	rainbow merge -o rbasm.out -a -i rbdiv.out -N500
```
The final output file of 'rainbow merge -a' is based on the final merged
clusters. Each cluster has been locally assembled by 'rainbow merge -a'. For
each cluster, rainbow outputs all assembled contigs seperated by '//' for each
record:
```
E clusterID
C contigID1
L length
S sequence
N #reads
R readIDs
//
C contigID2
L length
S sequence
N #reads
R readIDs
.
.
.
```

We have also provided four simple perl scripts that can be used to extract the assembly
information:
```
 select_all_rbcontig.pl, select_best_rbcontig.pl, select_sec_rbcontig.pl, select_best_rbcontig_plus_read1.pl
```
select_all_rbcontig.pl extracts all the assembled contigs, i.g., all the
records

select_best_rbcontig.pl and select_sec_rbcontig.pl extract the longest and
the longest plus the second longest contigs for the final clusters,
respectively

select_best_rbcontig_plus_read1.pl, as select_best_rbcontig.pl, it  extracts the longest contig for each cluster. Besides, it also outputs the read1. If read1 overlaps with the contig, it joins the two as a whole. If read1 does not overlap with the contig, it pads 10 'X' to join the read1 and the contig, thus generating a long contig. 

```
rainbow 2.0.3 -- <ruanjue@gmail.com, chongzechen@gmail.com>
Usage: rainbow <cmd> [options]

 cluster
  Input  File Format: paired fasta/fastq file(s)
  Output File Format: <seqid:int>\t<cluster_id:int>\t<read1:string>\t<read2:string>
  -1 <string> Input fasta/fastq file, supports multiple '-1'
  -2 <string> Input fasta/fastq file, supports multiple '-2' [null]
  -l <int>    Read length, default: 0 variable
  -m <int>    Maximum mismatches [4]
  -e <int>    Exactly matching threshold [2000]
  -L          Low level of polymorphism
 div
  Input File Format: <seqid:int>\t<cluster_id:int>\t<read1:string>\t<read2:string>
  Output File Format: <seqid:int>\t<cluster_id:int>\t<read1:string>\t<read2:string>[\t<pre_cluster_id:int>]
  -i <string> Input rainbow cluster output file [stdin]
  -o <string> Output file [stdout]
  -k <int>    K_allele, min variants to create a new group [2]
  -K <int>    K_allele, divide regardless of frequency when num of variants exceed this value [50]
  -f <float>  Frequency, min variant frequency to create a new group [0.2]
 merge
  Input File Format: <seqid:int>\t<cluster_id:int>\t<read1:string>\t<read2:string>[\t<pre_cluster_id:int>]
  -i <string> Input rainbow div output file [stdin]
  -a          output assembly 
  -o <string> Output file [stdout]
  -N <int>    Maximum number of divided clusters to merge [300]
  -l <int>    Minimum overlap when assemble two reads (valid only when '-a' is opened) [5]
  -f <float>  Minimum fraction of similarity when assembly (valid only when '-a' is opened) [0.90]
  -r <int>    Minimum number of reads to assemble (valid only when '-a' is opened) [5]
  -R <int>    Maximum number of reads to assemble (valid only when '-a' is opened) [300]

```
rbasm: a greedy assembler to locally assemble each cluster produced by rainbow or the other
tools. This has been integrated into the merge module. Please always open '-a' option when running
'rainbow merge'.
Local assemble fragments around restriction sites
***NOTE***: the input file format should be: <seqid:int>\t<cluster_id:int>\t<read1:string>\t<read2:string>[\t<pre_cluster_id:int>]
```
Usage: rbasm [options]
 -i <string> Input file [STDIN] 
 -o <string> Output file [STDOUT]
 -l <int>    Minium length of overlap [5]
 -s <float>  Minium similiarity of overlap [0.90]
 -r <int>    Minium reads to execute assembly [5]
 -R <int>    Maxium reads to execute assembly [200]
```

## Change log:
* v2.0.1: README and usage infomation updated
* v2.0.2: 'merge' options are riched. The 'merge' assembly work can be customized like rbasm now. Thanks Ross Whetten in NCSU for advicing this.
* v2.0.3: changed the name of script 'select_best_rbcontig2.pl' to 'select_best_rbcontig_plus_read1.pl', and documented it. 
* v2.0.4: fixed a bug that rainbow cannot be compiled in Mac OS
