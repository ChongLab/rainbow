/*
 * 
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#include "rainbow.h"

const char *version = "1.1";

int usage(){
	printf(
	"rainbow %s -- <ruanjue@gmail.com, chongzechen@gmail.com>\n"
	"Usage: rainbow <cmd> [options]\n"
	"Input  File Format: paired fasta/fastq file(s)\n"
	"Output File Format: <seqid:int>\\t<cluster_id:int>\\t<read1:string>\\t<read2:string>[\\t<pre_cluster_id:int>]\n"
	"\n"
	" cluster\n"
	"  -1 <string> Input fasta/fastq file, supports multiple '-1'\n"
	"  -2 <string> Input fasta/fasta file, supports multiple '-2' [null]\n"
	"  -l <int>    Read length, default: 0 variable\n"
	//"  -r <int>    rank of input files [1]\n"
	"  -m <int>    Maximum mismatches [2]\n"
	"  -e <int>    Exactly matching threshold [2000]\n"
	" div\n"
	"  -i <string> Input file [stdin]\n"
	"  -k <int>    K_allele, min variants to create a new group [2]\n"
	"  -K <int>    K_allele, divide regardless of frequency when num of variants exceed this value [50]\n"
	"  -f <float>  Frequency, min variant frequency to create a new group [0.2]\n"
	" merge \n"
	"  -a <string> Input rbasm output file [stdin]\n"
	"  -v <string> Input rainbow divided file [stdin]\n"
	"  -p <float>  maximum heterozygosity to collapse, should be specifed according to the estimated\n"
	"              polymorphism of the species [0.02]\n"
	"  -l <int>    Minimum overlap to collapse two contigs [100]\n"
	"  -k <int>    Minimum number of kmers to define similarity between two contigs [5]\n"
	"  -o <string> Output file for merged contigs, one line per cluster [stdout]\n" 
//	"  -n <int>    Maximum number of contigs to execute pairwise alignment [50]\n"
	"\n",
	version
	);
	return 1;
}

define_list(namelist, char*);

int cluster_invoker(int argc, char **argv){
	Cluster *cluster;
	FileReader *fr1, *fr2;
	namelist *list1, *list2;
	char *infile;
	int max_mm, c, exact_limit, is_fq1, is_fq2, fix_rd_len;
	int rank = 1;
	fr2 = NULL;
	infile = NULL;
	max_mm = 2;
	exact_limit = 2000;
	fix_rd_len = 0;
	list1 = init_namelist(2);
	list2 = init_namelist(2);

	while((c = getopt(argc, argv, "h1:2:r:m:e:l:")) != -1){
		switch(c){
			case 'h': return usage();
			case '1': push_namelist(list1, optarg); break;
			case '2': push_namelist(list2, optarg); break;
			case 'r': rank = atoi(optarg); break;
			case 'l': fix_rd_len = atoi(optarg); break;
			case 'm': max_mm = atoi(optarg); break;
			case 'e': exact_limit = atoi(optarg); break;
			default: return usage();
		}
	}
	if(count_namelist(list1) == 0) return usage();
	if(count_namelist(list2) != 0 && count_namelist(list1)!=count_namelist(list2)) {
		fprintf(stderr, "file1 and file2 should be paired\n\n");
		return usage();
	}
	is_fq1 = is_fq2 = 0;
	if((fr1 = fopen_m_filereader(count_namelist(list1), as_array_namelist(list1))) == NULL){
		fprintf(stderr, " -- Cannot open input file in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		abort();
	} else {
		is_fq1 = guess_seq_file_type(fr1);
		switch (is_fq1) {
			case 1: is_fq1 = 0; break;
			case 2: is_fq1 = 1; break;
			default: fprintf(stderr, "unknown file type\n");
			abort(); 
		}
	}
	if(count_namelist(list2) != 0) {
		if((fr2 = fopen_m_filereader(count_namelist(list2), as_array_namelist(list2))) == NULL){
			fprintf(stderr, " -- Cannot open input file in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
			abort();
		} else {
			is_fq2 = guess_seq_file_type(fr2);
			switch (is_fq2) {
				case 1: is_fq2 = 0; break;
				case 2: is_fq2 = 1; break;
				default: fprintf(stderr, "unknown file type\n");
				abort(); 
			}
		}
	}
	free_namelist(list1);
	free_namelist(list2);
	
	cluster = init_cluster(max_mm, exact_limit);
	indexing_cluster(cluster, fr1, is_fq1, fix_rd_len);
	clustering(cluster, fr2, is_fq2, fix_rd_len, stdout);
	free_cluster(cluster);
	fclose_filereader(fr1);
	if(fr2) fclose_filereader(fr2);
	fprintf(stderr, "Program exit normally\n");
	return 0;
}

int div_invoker(int argc, char **argv){
	Div *div;
	FileReader *fr;
	FILE *out;
	int c, k_allele, K_allele;
	float min_freq;
	char *infile, *outfile;
	infile = NULL;
	outfile = NULL;
	k_allele = 2;
	K_allele = 50;
	min_freq = 0.2;
	while((c = getopt(argc, argv, "hi:o:k:K:f:")) != -1){
		switch(c){
			case 'h': return usage();
			case 'i': infile = optarg; break;
			case 'o': outfile = optarg; break;
			case 'k': k_allele = atoi(optarg); break;
			case 'K': K_allele = atoi(optarg); break;
			case 'f': min_freq = atof(optarg); break;
			default: return usage();
		}
	}
	if(infile){
		if((fr = fopen_filereader(infile)) == NULL){
			fprintf(stdout, " -- Cannot open %s in %s -- %s:%d --\n", infile, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	} else fr = stdin_filereader();
	if(outfile){
		if((out = fopen(outfile, "w")) == NULL){
			fprintf(stdout, " -- Cannot write %s in %s -- %s:%d --\n", outfile, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	} else out = stdout;
	div = init_div(k_allele, K_allele, min_freq);
	div_reads(div, fr, out);
	free_div(div);
	fclose_filereader(fr);
	if(outfile) fclose(out);
	return 0;
}

int merge_invoker(int argc, char **argv) {
	FileReader *asmd, *divd;
	FILE *out = NULL;
	char *asmdf = NULL, *divdf = NULL, *outfile = NULL;
	uint32_t min_kmer = 5;
	uint32_t min_overlap = 100;
	float het = 0.02; int c;
	uint32_t kmersize = 15;

	while ((c = getopt(argc, argv, "ha:v:l:p:k:o:s:")) != -1) {
		switch (c) {
			case 'h': return usage();
			case 'a': asmdf = optarg; break;
			case 'v': divdf = optarg; break;
			case 'l': min_overlap = atoi(optarg); break;
			case 'p': het = atof(optarg); break;
			case 'k': min_kmer = atoi(optarg); break;
			case 'o': outfile = optarg; break;
			case 's': kmersize = atoi(optarg); break;
			default: return usage();
		}
	}

	if (asmdf) {
		if ((asmd = fopen_filereader(asmdf)) == NULL) {
			fprintf(stdout, " -- Cannot open %s in %s -- %s:%d --\n", asmdf, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	} else asmd = stdin_filereader();

	if (divdf) {
		if ((divd = fopen_filereader(divdf)) == NULL) {
			fprintf(stdout, " -- Cannot open %s in %s -- %s:%d --\n", divdf, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	} else divd = stdin_filereader();
	if (outfile) {
		if ((out = fopen(outfile, "w")) == NULL) {
			fprintf(stdout, " -- Cannot write %s in %s -- %s:%d --\n", divdf, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	} else out = stdout;
	merge_t *merger;
	merger = init_merger(min_kmer, min_overlap, het, kmersize);
	merge_ctgs(merger, asmd, divd, out);
	free_merger(merger);
	fclose_filereader(asmd);
	fclose_filereader(divd);
	if (outfile) fclose(out);
	return 0;
}

int main(int argc, char **argv){
	if(argc < 2) return usage();
	if(strcasecmp(argv[1], "cluster") == 0){
		return cluster_invoker(argc - 1, argv + 1);
	} else if(strcasecmp(argv[1], "div") == 0){
		return div_invoker(argc - 1, argv + 1);
	} else if(strcasecmp(argv[1], "merge") == 0) {
		return merge_invoker(argc - 1, argv + 1);
	} else {
		return usage();
	}
}
