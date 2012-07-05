#ifndef MERGECTG_H
#define MERGECTG_H

#include <stdint.h>
#include "list.h"
#include "file_reader.h"
#include "hashset.h"
#include "string.h"
#include "stdaln.h"
#include "asm_R2.h"
#include "bloom_filter.h"

typedef struct {
	char seq[MAX_RD_LEN+1];
	uint32_t seq_id;
	uint32_t rd_len;
	uint32_t rank;
} read_t;

define_list(readv, read_t);

typedef struct {
	uint64_t kmer:62, kpos:2;
} rd_kmer_t;

#define rd_kmer_code(r) u32hashcode((r).kmer)
#define rd_kmer_eq(r1, r2) ((r1).kmer == (r2).kmer)
define_hashset(rdkhash, rd_kmer_t, rd_kmer_code, rd_kmer_eq);
define_list(idxv, rdkhash*);
//define_list(idxv, BloomFilter*);

typedef struct {
//	uint32_t id, clsid, old_clsid, sz;
	uint32_t id;
	int closed;
//	char *seq, *sec_seq;
	String *path;
	readv *rds;
	u32list *m_rds;  // merged reads index
	rdkhash *index;
//	BloomFilter *index;
	idxv *m_idx;  // merged multiple index
//	Vector *efctgs;
} contig_t; 

#define contig_code(c) u32hashcode((c).id)
#define contig_eq(c1, c2) ((c1).id == (c2).id)
define_hashset(ctgset, contig_t, contig_code, contig_eq);

typedef struct {
	uint32_t id;
	char *seq;
} contig_seq_t;

define_list(contigv, contig_t*);
define_list(contigsv, contig_t*);

typedef struct pathtree_t pathtree_t;
struct pathtree_t {
	uint32_t tid; // leaf records contig ID
	pathtree_t *left;
	pathtree_t *right;
};

typedef struct {
	uint64_t kmer, kpos;
	uint32_t id; // which contig
	int offset;  // offset w.r.t. the current contig
	int offset2; // offset of query contig
} ctg_kmer_t;

#define kmer_code(k) u64hashcode((k).kmer)
#define kmer_eq(k1, k2) ((k1).kmer == (k2).kmer)
define_hashset(ctgkhash, ctg_kmer_t, kmer_code, kmer_eq);

typedef struct {uint32_t key; uint32_t oldid; char *path;} uuchash_t;
#define uuchash_code(e) (e).key
#define uuchash_equals(e1, e2) ((e1).key == (e2).key)
define_hashset(uuchash, uuchash_t, uuchash_code, uuchash_equals);

define_list(ctgkmerv, ctg_kmer_t);

typedef struct {
	uint64_t last; //last kmer position
	int offset; // current kmer offset
} link_t;

define_search_array(bisearch, uint64_t, native_number_cmp);

typedef struct {
	contigv *ctgs;
	contigsv *cache;
	pathtree_t *tree;
	uint32_t min_kmer; // parameter: # kmers to define two similar contigs
	uint32_t min_overlap; // parameter
	float het; // parameter
	uint32_t RD_KMER_SIZE; // parameter
	uint32_t min_ol; //parameter for asm
	float min_sm; // parameter for asm
	uint32_t min_read; // parameter for asm
	uint32_t max_read; // parameter for asm
	uint32_t sim_pairs;
	uint32_t max_cluster; //parameter
	uint32_t need_asm; // parameter
	uint32_t cid; //
	EF *ef;
	int flag;  // if == 0 first use, init; else reset
} merge_t;


#ifdef __CPLUSPLUS
extern "C" {
#endif

merge_t* init_merger(uint32_t min_kmer, uint32_t min_overlap, float het, uint32_t kmersize, uint32_t max_cluster, uint32_t need_asm);
//void merge_ctgs(merge_t *merger, FileReader *asmd, FileReader *divd, FILE *out);
void merge_ctgs(merge_t *merger, FileReader *in, FILE *out);
void merge_along_tree(merge_t *merger, pathtree_t *tree);
void merge_core(merge_t *merger);
void free_index(merge_t *merger);
void free_ctg(contig_t *ctg);
void free_ctgs(merge_t *merger);
void build_tree(merge_t *merger);
void update_ctg2merge(merge_t *merger);
int is_similar_enough(merge_t *merger, contig_t *c1, contig_t *c2);
void merge_2ctg(merge_t *merger, contig_t *ctg1, contig_t *ctg2);
void update_merger(merge_t *merger, contig_t *ctg1, contig_t *ctg2);
void prefix_path(char *s1, char *s2, int n, char *pre);
void destroy_tree(pathtree_t *t);
void free_tree(merge_t *merger);
void reset_merger(merge_t *merger);
void free_merger(merge_t *merger);
void clear_ctg(contig_t *ctg);

#ifdef __CPLUSPLUS
}
#endif

#endif
