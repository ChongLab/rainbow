#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <assert.h>

#define PACKAGE_VERSION "0.1.1"

uint8_t nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
static mut_t mutmsk = (mut_t)0xf000;

typedef struct {
	int l, m; /* length and maximum buffer size */
	unsigned char *s; /* sequence */
} seq_t;

typedef struct {
	int l, m; /* length and maximum buffer size */
	mut_t *s; /* sequence */
} mutseq_t;

typedef struct {
	uint64_t l, m;
	uint64_t *idx;
} idx_t;

#define INIT_SEQ(seq) (seq).s = 0; (seq).l = (seq).m = 0
#define INIT_IDX(index) (index).idx = 0; (index).l = (index).m = 0

static int SEQ_BLOCK_SIZE = 512;

void seq_set_block_size(int size)
{
	SEQ_BLOCK_SIZE = size;
}

int seq_read_fasta(FILE *fp, seq_t *seq, char *locus, char *comment)
{
	int c, l, max;
	char *p;
	
	c = 0;
	while (!feof(fp) && fgetc(fp) != '>');
	if (feof(fp)) return -1;
	p = locus;
	while (!feof(fp) && (c = fgetc(fp)) != ' ' && c != '\t' && c != '\n')
		if (c != '\r') *p++ = c;
	*p = '\0';
	if (comment) {
		p = comment;
		if (c != '\n') {
			while (!feof(fp) && ((c = fgetc(fp)) == ' ' || c == '\t'));
			if (c != '\n') {
				*p++ = c;
				while (!feof(fp) && (c = fgetc(fp)) != '\n')
					if (c != '\r') *p++ = c;
			}
		}
		*p = '\0';
	} else if (c != '\n') while (!feof(fp) && fgetc(fp) != '\n');
	l = 0; max = seq->m;
	while (!feof(fp) && (c = fgetc(fp)) != '>') {
		if (isalpha(c) || c == '-' || c == '.') {
			if (l + 1 >= max) {
				max += SEQ_BLOCK_SIZE;
				seq->s = (unsigned char*)realloc(seq->s, sizeof(char) * max);
			}
			seq->s[l++] = (unsigned char)c;
		}
	}
	if (c == '>') ungetc(c,fp);
	seq->s[l] = 0;
	seq->m = max; seq->l = l;
	return l;
}

/* Error-checking open, copied from utils.c */

#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
		abort();
	}
	return fp;
}

static double ERR_RATE = 0.02;
static double DEPTH = 10.0;
static double MUT_RATE = 0.001;
static double HOM_RATE = 0.0;
static double INDEL_FRAC = 0.1;
static double INDEL_EXTEND = 0.3;

void uc(unsigned char *s)
{
	while (*s) {
		*s = toupper(*s);
		s++;
	}
}


int strindex(idx_t *index, unsigned char *s, unsigned char *t)
{
	uint64_t i, j, k;
	uint64_t l, max;

	i = 0;
	while (t[i]) {
		if (nst_nt4_table[(int)t[i]] == 4) {
			//printf("here%c\n", *t);
			return -1;
		}
		i++;
	}
	
	l = max = 0;
	for (i = 0; s[i] != '\0'; i++) {
		for (j=i, k=0; t[k]!='\0' && s[j]==t[k]; j++, k++)
			;
		if (k > 0 && t[k] == '\0') {
			if (l + 1 >= max) {
		  		max += SEQ_BLOCK_SIZE;
				index->idx = (uint64_t*)realloc(index->idx, sizeof(uint64_t) * max);
			}
			index->idx[l++] = i;
		}
	}

	if (l) {
		index->l = l;
		index->m = max;
		index->idx[l] = -1;
		return l;
	}
	else 
		return -1;
}

/* Simple normal random number generator, copied from genran.c */

double ran_normal()
{ 
	static int iset = 0; 
	static double gset; 
	double fac, rsq, v1, v2; 
	if (iset == 0) {
		do { 
			v1 = 2.0 * drand48() - 1.0;
			v2 = 2.0 * drand48() - 1.0; 
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq); 
		gset = v1 * fac; 
		iset = 1;
		return v2 * fac;
	} else {
		iset = 0;
		return gset;
	}
}
void maq_mut_diref(const seq_t *seq, int is_hap, mutseq_t *hap1, mutseq_t *hap2)
{
	int i, deleting = 0;
	mutseq_t *ret[2];

	ret[0] = hap1; ret[1] = hap2;
	ret[0]->l = seq->l; ret[1]->l = seq->l;
	ret[0]->m = seq->m; ret[1]->m = seq->m;
	ret[0]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));
	ret[1]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));
	for (i = 0; i != seq->l; ++i) {
		int c;
		c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)seq->s[i]];
        if (deleting) {
            if (drand48() < INDEL_EXTEND) {
                if (deleting & 1) ret[0]->s[i] |= DELETE;
                if (deleting & 2) ret[1]->s[i] |= DELETE;
                continue;
            } else deleting = 0;
        }
		if (c < 4 && drand48() < MUT_RATE) { // mutation
			if (drand48() >= INDEL_FRAC) { // substitution
				double r = drand48();
				c = (c + (int)(r * 3.0 + 1)) & 3;
				if (is_hap || drand48() < HOM_RATE) { // hom
					ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
				} else { // het
					ret[drand48()<0.5?0:1]->s[i] = SUBSTITUTE|c;
				}
			} else { // indel
				if (drand48() < 0.5) { // deletion
					if (is_hap || drand48() < HOM_RATE) { // hom-del
						ret[0]->s[i] = ret[1]->s[i] = DELETE;
                        deleting = 3;
					} else { // het-del
                        deleting = drand48()<0.5?1:2;
						ret[deleting-1]->s[i] = DELETE;
					}
				} else { // insertion
                    int num_ins = 0, ins = 0;
                    do {
                        num_ins++;
                        ins = (ins << 2) | (int)(drand48() * 4.0);
                    } while(num_ins < 4 && drand48() < INDEL_EXTEND);

					if (is_hap || drand48() < HOM_RATE) { // hom-ins
						ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
					} else { // het-ins
						ret[drand48()<0.5?0:1]->s[i] = (num_ins << 12) | (ins << 4) | c;
					}
				}
			}
		}
	}
}
void maq_print_mutref(const char *name, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2)
{
	int i;
	for (i = 0; i != seq->l; ++i) {
		int c[3];
		c[0] = nst_nt4_table[(int)seq->s[i]];
		c[1] = hap1->s[i]; c[2] = hap2->s[i];
		if (c[0] >= 4) continue;
		if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
			printf("%s\t%d\t", name, i+1);
			if (c[1] == c[2]) { // hom
				if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
					printf("%c\t%c\t-\n", "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
				} else if ((c[1]&mutmsk) == DELETE) { // del
					printf("%c\t-\t-\n", "ACGTN"[c[0]]);
				} else if (((c[1] & mutmsk) >> 12) <= 5) { // ins
					printf("-\t");
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while(n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        n--;
                    }
                    printf("\t-\n");
				}  else assert(0);
			} else { // het
				if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
					printf("%c\t%c\t+\n", "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)]);
				} else if ((c[1]&mutmsk) == DELETE) {
					printf("%c\t-\t+\n", "ACGTN"[c[0]]);
				} else if ((c[2]&mutmsk) == DELETE) {
					printf("%c\t-\t+\n", "ACGTN"[c[0]]);
				} else if (((c[1] & mutmsk) >> 12) <= 4) { // ins1
					printf("-\t");
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        n--;
                    }
                    printf("\t+\n");
				} else if (((c[2] & mutmsk) >> 12) <= 5) { // ins2
					printf("-\t");
                    int n = (c[2]&mutmsk) >> 12, ins = c[2] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t+\n");
				} else assert(0);
			}
		}
	}
}


//wiki knuth method
int poisson_num_gen(double lamda)
{
	int k = 0;
	double L = exp(-lamda);
	double p = 1.0;
	
	do {
		k++;
		p = p * drand48();
	} while (p > L);

	return k-1;
}

void ezmsim_LR_core(FILE *fpout1, FILE *fpout2, FILE *fp_fa, int size_l, int size_r, unsigned char *cut, int pos)
{
	idx_t index;
	seq_t seq;
	uint64_t tot_len, dep;
	unsigned int i, k;
	int len, n_ref, j, size[2], Q, m;
	char name[256], *qstr;
	uint8_t *tmp_seq[2], c;
	uint64_t id;

	INIT_SEQ(seq);
	INIT_IDX(index);

	Q = (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;

	srand48(time(0));
	seq_set_block_size(0x1000000);
	len = size_l > size_r?size_l:size_r;
	qstr = (char*)calloc(len+1, 1);
	tmp_seq[0] = (uint8_t*)calloc(len+2, 1);
	tmp_seq[1] = (uint8_t*)calloc(len+2, 1);
	size[0] = size_l; size[1] = size_r;
	tot_len = n_ref = 0; id = 0;
	while ((len = seq_read_fasta(fp_fa, &seq, name, 0)) >= 0) {
		uc(seq.s);
		uc(cut);
		if (strindex(&index, seq.s, cut) != -1)
			printf("chromsome %s has %llu digest sites\n", name, (unsigned long long)index.l);
		
		for (i = 0; i < index.l; i++) {
			if (i - size_l <= 0 || i + size_r > (unsigned int)len)
				continue;
			dep = poisson_num_gen(DEPTH);
			for (k = 0; k < dep; k++) {
				FILE *fpo[2];
				int is_flip = 0, s[2];
				id++;

				if (drand48() < 0.5) {
					fpo[0] = fpout1; fpo[1] = fpout2;
					s[0] = size[0]; s[1] = size[1];
				} else {
					fpo[1] = fpout1; fpo[0] = fpout2;
					s[1] = size[0]; s[0] = size[1];
					is_flip = 1;
				}

				for (j = 0; j < s[0]; j++) {
					c = nst_nt4_table[(int)seq.s[index.idx[i]+pos-j-1]];
					if (c >= 4) c = 4;
					else if (drand48() < ERR_RATE) {
						c = (c + (int)(drand48()*3.0 + 1)) & 3;
					}
					tmp_seq[0][j] = c;
				}

				for (j = 0; j < s[1]; j++) {
					c = nst_nt4_table[(int)seq.s[index.idx[i]+pos+j]];
					if (c >= 4) c = 4;
					else if (drand48() < ERR_RATE) {
						c = (c + (int)(drand48()*3.0 + 1)) & 3;
					}
					tmp_seq[1][j] = c < 4?3-c:4;
				}

				for (m = 0; m < 2; m++) {
					fprintf(fpo[m], "@%s_%s_%llu_%llu/%d\n", name, cut, (unsigned long long)index.idx[i], (unsigned long long)id, m==0?is_flip+1:2-is_flip);
					for (j = 0; j < s[m]; j++) {
						qstr[j] = Q;
						fputc("ACGTN"[(int)tmp_seq[m][j]], fpo[m]);
					} 
					qstr[j] = 0;
					fprintf(fpo[m], "\n+\n%s\n", qstr);
				}
			}
			//fprintf(stderr, "%llu ", index.idx[i]);
		}
		printf("\n");
		
		tot_len += len;
		++n_ref;
	}
	fprintf(stderr, "-- %d sequences, total length: %llu\n", n_ref, (unsigned long long)tot_len);
	rewind(fp_fa);
	
	free(seq.s);
	free(index.idx);
}

int LR_usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: ezmsim (simulate enzyme cut assembly sequences)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Zechen Chong <chongzechen@gmail.com>\n\n");
	fprintf(stderr, "Usage: ezmsim LR [options] <-z enzyme> <in.ref.fa> <out.read1.fq> <out.read2.fq>\n\n");
	fprintf(stderr, "Options: -e FLOAT    base error rate [%.3f]\n", ERR_RATE);
	fprintf(stderr, "         -D FLOAT    read depth [%.1f]\n", DEPTH);
	fprintf(stderr, "         -1 INT      length of the first read [100]\n");
	fprintf(stderr, "         -2 INT      length of the second read [100]\n");
	fprintf(stderr, "         -z STRING   enzyme sequence (must be specified)\n");
	fprintf(stderr, "         -p INT      enzyme cut position [length(enzyme)/2]\n");
	fprintf(stderr, "\n");
	return 1;
}

int LR_main(int argc, char *argv[])
{
	unsigned int size_l, size_r;
	FILE *fpout1, *fpout2, *fp_fa;
	char *cut = "";

	int c, pos; 
	pos = -1;
	//dist = 250; 
	//std_dev = 20;
	size_l = size_r = 100;
	while ((c = getopt(argc, argv, "e:D:1:2:d:s:z:p:")) != -1) {
	switch (c) {
		case 'e': ERR_RATE = atof(optarg); break;
		case 'D': DEPTH = atof(optarg); break;
		case '1': size_l = atoi(optarg); break;
		case '2': size_r = atoi(optarg); break;
		case 'z': cut = strdup(optarg);
		case 'p': pos = atoi(optarg);
				  uc((unsigned char*)cut); break;
		//case 'd': dist = atoi(optarg); break;
		//case 's': std_dev = atoi(optarg); break;
		default:
				  //printf("unknown option: %c\n", optopt);
				  //return 1;
				  return LR_usage();
		}
	}
	if (strlen(cut)==0) {
		//fprintf(stderr, "parameter z (enzyme cut) must be specified\n");
		return LR_usage();
	}
	if (argc - optind < 3) {
		fprintf(stderr, "files must be specified\n");
		return LR_usage();
	}
	fp_fa = (strcmp(argv[optind+0], "-") == 0)?stdin:xopen(argv[optind+0], "r");
	fpout1 = xopen(argv[optind+1], "w");
	fpout2 = xopen(argv[optind+2], "w");

	if (pos == -1) {
		pos = strlen(cut)/2;
	}
	ezmsim_LR_core(fpout1, fpout2, fp_fa, size_l, size_r, (unsigned char*)cut, pos);
	
	fclose(fp_fa); fclose(fpout1); fclose(fpout2);
	return 0;
}

void ezmsim_EF_core(FILE *fpout1, FILE *fpout2, FILE *fp_fa, unsigned int size_l, unsigned int size_r, unsigned char *cut, int pos, int distance, int ovlp, int stp, int reverse, int is_hap)
{
	idx_t index;
	seq_t seq;
	uint64_t tot_len, dep, i, k;
	int len, n_ref, j, size[2], Q, m, n;
	char name[256], *qstr;
	uint8_t *tmp_seq[2], c;
	uint64_t id;
	int dist, overlap, step, rev;
	mutseq_t rseq[2];
	mut_t *target;

	INIT_SEQ(seq);
	INIT_IDX(index);

	Q = (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;

	srand48(time(0));
	seq_set_block_size(0x1000000);
	len = size_l > size_r?size_l:size_r;
	qstr = (char*)calloc(len+1, 1);
	tmp_seq[0] = (uint8_t*)calloc(len+2, 1);
	tmp_seq[1] = (uint8_t*)calloc(len+2, 1);
	size[0] = size_l; size[1] = size_r;
	tot_len = n_ref = 0; id = 0;
	dist = distance; overlap = ovlp; step = stp; rev = reverse;

	while ((len = seq_read_fasta(fp_fa, &seq, name, 0)) >= 0) {
		uc(seq.s);
		uc(cut);
		if (strindex(&index, seq.s, cut) != -1)
			printf("chromsome %s has %llu digest sites\n", name, (unsigned long long)index.l);
		
		maq_mut_diref(&seq, is_hap, rseq, rseq+1);
		maq_print_mutref(name, &seq, rseq, rseq+1);

		for (i = 0; i < index.l; i++) {
			if (len < (dist + overlap*step) * 2) {
				fprintf(stderr, "[ezmsim_core] skip sequence '%s' as it is shorter than %d!\n", name, (dist + overlap*step)*2);
				continue;
			}
			
			for (n = 0; n < step; n++) {
				dep = poisson_num_gen(DEPTH);
				for (k = 0; k < dep; k++) {
					FILE *fpo[2];
					int is_flip = 0, s[2], d;
					id++;

					d = dist + (int)(drand48()*overlap);

					//if (drand48() < 0.5) {
					if (!rev) {
						fpo[0] = fpout1; fpo[1] = fpout2;
						s[0] = size[0]; s[1] = size[1];
					} else {
						fpo[1] = fpout1; fpo[0] = fpout2;
						s[1] = size[0]; s[0] = size[1];
						is_flip = 1;
					}
					//generate the read sequences
					target = rseq[drand48()<0.5?0:1].s;
					int ii, begin, end;
					for (ii = index.idx[i]+pos, j = 0, begin = 0; ii < seq.l && j < s[0]; ++ii) {
						int c = target[ii];
						int mut_type = c & mutmsk;
						if (mut_type == DELETE) continue; // deletion
						if (begin == 0) {
							begin = ii;
							if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) mut_type = NOCHANGE; // skip ins at the first base
						}
						if(mut_type == NOCHANGE || mut_type == SUBSTITUTE) {
							tmp_seq[0][j++] = c&0xf;
							continue;
						}
						int n = mut_type >> 12, ins = c >> 4;
						while (n > 0) {
							tmp_seq[0][j++] = ins & 0x3;
							ins >>= 2;
							n--;
							if ((int)k == s[0]) break;
						}
						tmp_seq[0][j++] = c&0xf;
					}
					for (ii = index.idx[i]+pos+d-1, j = 0, end = 0; ii>=0 && j < s[1];--ii) {
						int c = target[ii];
						if ((c&mutmsk) == DELETE) continue; // deletion
						if (end == 0) end = i;
						tmp_seq[1][j++] = c&0xf;
						if((c&mutmsk) == NOCHANGE || (c&mutmsk) == SUBSTITUTE) continue;
						int n = (c&mutmsk) >> 12, ins = c >> 4;
						while (n > 0) {
							if (j == s[1]) break;
							tmp_seq[1][j++] = ins & 0x3;
							ins >>= 2;
							n--;
						}
					}
					for (j = 0; j < s[0]; j++) {
						c = tmp_seq[0][j];
						//c = nst_nt4_table[(int)seq.s[index.idx[i]+pos+j]];
						if (c >= 4) c = 4;
						else if (drand48() < ERR_RATE) {
							c = (c + (int)(drand48()*3.0+1)) & 3;
						}
						tmp_seq[0][j] = c;
					}

					for (j = 0; j < s[1]; j++) {
						c = tmp_seq[1][j];
						//c = nst_nt4_table[(int)seq.s[index.idx[i]+pos+d-j]];
						if (c >= 4) c = 4;
						else if (drand48() < ERR_RATE) {
							c = (c + (int)(drand48()*3.0 + 1)) & 3;
						}
						tmp_seq[1][j] = c<4?3-c:4;
					}
					
					for (m = 0; m < 2; m++) {
						fprintf(fpo[m], "@%s_%s_%llu_%llu_%d/%d\n", name, cut, (unsigned long long)index.idx[i], (unsigned long long)id, d, m==0?is_flip+1:2-is_flip);
						for (j = 0; j < s[m]; j++) {
							qstr[j] = Q;
							fputc("ACGTN"[(int)tmp_seq[m][j]], fpo[m]);
						} 
						qstr[j] = 0;
						fprintf(fpo[m], "\n+\n%s\n", qstr);
					}

				}
				/*
				dep = poisson_num_gen(DEPTH);
				for (k = 0; k < dep; k++) {
					FILE *fpo[2];
					int is_flip = 0, s[2], d;
					id++;

					d = dist + (int)(drand48()*step);

					//if (drand48() < 0.5) {
					if (!rev) {
						fpo[0] = fpout1; fpo[1] = fpout2;
						s[0] = size[0]; s[1] = size[1];
					} else {
						fpo[1] = fpout1; fpo[0] = fpout2;
						s[1] = size[0]; s[0] = size[1];
						is_flip = 1;
					}

					for (j = 0; j < s[0]; j++) {
						c = nst_nt4_table[(int)seq.s[index.idx[i]+pos-j-1]];
						if (c >= 4) c = 4;
						else if (drand48() < ERR_RATE) {
							c = (c + (int)(drand48()*3.0+1)) & 3;
						}
						tmp_seq[0][j] = c;
					}

					for (j = 0; j < s[1]; j++) {
						c = nst_nt4_table[(int)seq.s[index.idx[i]+pos-d+j]];
						if (c >= 4) c = 4;
						else if (drand48() < ERR_RATE) {
							c = (c + (int)(drand48()*3.0 + 1)) & 3;
						}
						tmp_seq[1][j] = c<4?3-c:4;
					}
					
					for (m = 0; m < 2; m++) {
						fprintf(fpo[m], "@%s_%s_%lld_%lld_%d/%d\n", name, cut, index.idx[i], id, d, m==0?is_flip+1:2-is_flip);
						for (j = 0; j < s[m]; j++) {
							qstr[j] = Q;
							fputc("ACGTN"[(int)tmp_seq[m][j]], fpo[m]);
						} 
						qstr[j] = 0;
						fprintf(fpo[m], "\n+\n%s\n", qstr);
					}
				}*/
				dist += overlap;
			}
			dist = distance;
		}
		tot_len += len;
		++n_ref;
	}
	fprintf(stderr, "-- %d sequences, total length: %llu\n", n_ref, (unsigned long long)tot_len);
	
	free(seq.s);
	free(index.idx);
}

int EF_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: ezmsim (simulate enzyme cut assembly sequences)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Zechen Chong <chongzechen@gmail.com>\n\n");
	fprintf(stderr, "Usage: ezmsim RAD [options] <-z enzyme> <in.ref.fa> <out.read1.fq> <out.read2.fq>\n\n");
	fprintf(stderr, "Options: -e FLOAT    base error rate [%.3f]\n", ERR_RATE);
	fprintf(stderr, "         -D FLOAT    read depth [%.1f]\n", DEPTH);
	fprintf(stderr, "         -1 INT      length of the first read [100]\n");
	fprintf(stderr, "         -2 INT      length of the second read [100]\n");
	fprintf(stderr, "         -z STRING   enzyme sequence (must be specified)\n");
	fprintf(stderr, "         -p INT      enzyme cut position [length(enzyme)/2]\n");
	fprintf(stderr, "         -d INT      initial insert size distance [120]\n");
	//fprintf(stderr, "         -s INT      standard deviation of insert size [20]\n");
	fprintf(stderr, "         -o INT      insert size overlap distance [50]\n");
	fprintf(stderr, "         -t INT      elongation steps of insert size [10]\n");
	fprintf(stderr, "         -h FLOAT    rate of homozygosity[%.4f]\n", HOM_RATE);
	fprintf(stderr, "         -m FLOAT    rate of mutation[%.4f]\n", MUT_RATE);
	fprintf(stderr, "         -R FLOAT    fraction of indels [%.2f]\n", INDEL_FRAC);
	fprintf(stderr, "         -X FLOAT    probability an indel is extended [%.2f]\n", INDEL_EXTEND);
	fprintf(stderr, "         -r          reverse or not [forward only]\n");
	fprintf(stderr, "         -H          haploid mode\n");
	fprintf(stderr, "\n");
	return 1;
}

int EF_main(int argc, char *argv[])
{
	int c, size_l, size_r, pos, dist, overlap, step, rev, is_hap = 0;
	FILE *fpout1, *fpout2, *fp_fa;
	char *cut = "";

	pos = -1;
	dist = 120;  //initial distance
	//std_dev = 20; 
	overlap = 50;
	size_l = size_r = 100;
	step = 10; rev = 0;
	while ((c = getopt(argc, argv, "e:D:1:2:d:s:z:p:o:t:R:rh:Hm:")) != -1) {
	switch (c) {
		case 'e': ERR_RATE = atof(optarg); break;
		case 'D': DEPTH = atof(optarg); break;
		case '1': size_l = atoi(optarg); break;
		case '2': size_r = atoi(optarg); break;
		case 'z': cut = strdup(optarg);
		case 'p': pos = atoi(optarg);
				  uc((unsigned char*)cut); break;
		case 'd': dist = atoi(optarg); break;
		//case 's': std_dev = atoi(optarg); break;
		case 'o': overlap = atoi(optarg); break;
		case 't': step = atoi(optarg); break;
		case 'r': rev = 1; break;
		case 'h': HOM_RATE = atof(optarg); break;
		case 'H': is_hap = 1; break;
		case 'm': MUT_RATE = atof(optarg); break;
		case 'R': INDEL_FRAC = atof(optarg); break;
		case 'X': INDEL_EXTEND = atof(optarg); break;
		default:
				  //printf("unknown option: %c\n", optopt);
				  //return 1;
				  return EF_usage();
		}
	}
	if (strlen(cut)==0) {
		//fprintf(stderr, "parameter z (enzyme cut) must be specified\n");
		return EF_usage();
	}
	if (argc - optind < 3) {
		fprintf(stderr, "files must be specified\n");
		return EF_usage();
	}
	fp_fa = (strcmp(argv[optind+0], "-") == 0)?stdin:xopen(argv[optind+0], "r");
	fpout1 = xopen(argv[optind+1], "w");
	fpout2 = xopen(argv[optind+2], "w");

	if (pos == -1) {
		pos = strlen(cut)/2;
	}
	ezmsim_EF_core(fpout1, fpout2, fp_fa, size_l, size_r, (unsigned char*)cut, pos, dist, overlap, step, rev, is_hap);
	
	fclose(fp_fa); fclose(fpout1); fclose(fpout2);
	return 0;
}

int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: ezmsim (simulate enzyme cut assembly sequences)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Zechen Chong <chongzechen@gmail.com>\n\n");
	fprintf(stderr, "Usage: ezmsim <LR|EF> [options]\n\n");
	fprintf(stderr, "Options: LR          simulate LR reads\n");
	fprintf(stderr, "         RAD          simulate RAD reads\n");
	fprintf(stderr, "\n");
	return 1;
}

int main (int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "LR") == 0) return LR_main(argc-1, argv+1);
	else if (strcmp(argv[1], "RAD") == 0) return EF_main(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
