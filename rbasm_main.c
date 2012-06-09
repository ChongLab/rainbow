#include "asm_R2.h"

int main(int argc, char **argv){
	FileReader *in;
	FILE *out;
	uint32_t n_ef, min_ol, min_read, max_read;
	float min_sm;
	char *infile, *outfile;
	int c;
	infile = NULL;
	outfile = NULL;
	min_ol = 5;
	min_sm = 0.9;
	min_read = 5;
	max_read = 200;
	while((c = getopt(argc, argv, "hi:o:r:R:l:s:")) != -1){
		switch(c){
			case 'i': infile = optarg; break;
			case 'o': outfile = optarg; break;
			case 'l': min_ol = atoi(optarg); break;
			case 's': min_sm = atof(optarg); break;
			case 'r': min_read = atoi(optarg); break;
			case 'R': max_read = atoi(optarg); break;
			case 'h': return ef_usage();
		}
	}
	if(infile == NULL) in = stdin_filereader();
	else if((in = fopen_filereader(infile)) == NULL){
		fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", infile, __FUNCTION__, __FILE__, __LINE__);
		abort();
	}
	if(outfile == NULL) out = stdout;
	else if((out = fopen(outfile, "w")) == NULL){
		fprintf(stderr, " -- Cannot write %s in %s -- %s:%d --\n", outfile, __FUNCTION__, __FILE__, __LINE__);
		abort();
	}
	n_ef = asm_ef(in, out, min_ol, min_sm, min_read, max_read);
	fclose_filereader(in);
	if(outfile) fclose(out);
	return 0;
}
