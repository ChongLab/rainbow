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
#ifndef ASM_R2_H
#define ASM_R2_H

#include "string.h"
#include "vector.h"
#include "hashset.h"
#include "file_reader.h"
#include "dna.h"
#include <stdint.h>
#include <unistd.h>

#define MAX_RD_LEN	255
#define ASM_KMER_SIZE	5
#define ASM_KMER_MASK	0x3FFu

typedef struct { uint32_t rid:16, roff:16; } rp_t;

typedef struct {
	uint32_t kmer:10, rps_idx:22;
} rhash_t;

#define rhash_code(r) (r).kmer
#define rhash_eq(r1, r2) ((r1).kmer == (r2).kmer)
define_hashset(rhash, rhash_t, rhash_code, rhash_eq);

typedef struct {
	char     seq[MAX_RD_LEN+1];
	uint32_t seq_id;
	uint32_t rd_len:10, rank:10;
	uint32_t ctg_id:24, ctg_off:19, used:1;
} FRead;

typedef struct {
	uint32_t len:31, closed:1;
	String   *seq;
	Vector   *rids;
} FContig;

typedef struct {
	uint32_t l_rid:20, r_rid:20, l_ol:8, r_ol:8, n_mm:7, used:1;
} Overlap;

typedef struct {
	uint32_t ef_id;
	char     eseq[MAX_RD_LEN];
	Vector   *rds;
	Vector   *ctgs;
	Vector   *rps;
	Vector   *ols;
	rhash    *index;
	u64hash  *uniq;
	uint32_t min_ol;
	float    min_sm;
	uint32_t inc_tag;

	Vector   *pool_vec;
	Vector   *pool_ctg;
} EF;


#ifdef __CPLUSPLUS
extern "c" {
#endif


Vector* get_pool_vec(EF *ef);
void put_pool_vec(EF *ef, Vector *vec);
FContig* get_pool_ctg(EF *ef);
int cmp_ol_func(const void *e1, const void *e2);
void put_pool_ctg(EF *ef, FContig *ctg);
void add_read2ef_core(EF *ef, char *seq, uint32_t seq_id, uint32_t rd_len, uint32_t rank);
EF* init_ef(uint32_t ef_id, char *eseq, uint32_t rd_len, uint32_t min_ol, float min_sm); 
void set_inc_tag_ef(EF *ef, uint32_t inc_tag);
void add_read2ef(EF *ef, char *seq, uint32_t seq_id, uint32_t rd_len, uint32_t rank); 
void find_overlap(char *seq1, uint32_t len1, uint32_t off1, char *seq2, uint32_t len2, uint32_t off2, uint32_t *l_ol, uint32_t *r_ol, uint32_t *n_mm);
void align_reads_ef(EF *ef);
void print_alignments(EF *ef);
void asm_ef_ctgs(EF *ef);
void output_ef_ctgs(EF *ef, FILE *out);
void reset_ef(EF *ef, uint32_t ef_id, char *eseq, uint32_t rd_len, uint32_t min_ol, float min_sm);
void free_ef(EF *ef);
int ef_usage(void );
uint32_t asm_ef(FileReader *in, FILE *out, uint32_t min_ol, float min_sm, uint32_t min_read, uint32_t max_read);

#ifdef __CPLUSPLUS
}
#endif

#endif
