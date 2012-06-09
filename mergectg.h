#ifndef MERGECTG_H
#define MERGECTG_H

#include <stdint.h>
#include "list.h"
#include "file_reader.h"
#include "hashset.h"
#include "string.h"
#include "stdaln.h"
#include "asm_R2.h"

typedef struct {
	char seq[MAX_RD_LEN+1];
	uint32_t seq_id;
	uint32_t rd_len;
	uint32_t rank;
} read_t;

define_list(readv, read_t);

typedef struct {
//	uint32_t id, clsid, old_clsid, sz;
	uint32_t id;
	int closed;
	char *seq, *sec_seq, *path;
	readv *rds;
	Vector *efctgs;
} contig_t; 

#define contig_code(c) u32hashcode((c).id)
#define contig_eq(c1, c2) ((c1).id == (c2).id)
define_hashset(ctgset, contig_t, contig_code, contig_eq);

typedef struct {
	uint32_t id;
	char *seq;
} contig_seq_t;

define_list(contigv, contig_t);

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
//	ctgset *ctgs;
	pathtree_t *tree;
	ctgkhash *index;
	link_t *links;
	uint64_t *idv;  // to search translate pos to ids
	ctgkmerv *kmers;   // all kmers of each kmer of query contig
	ctgkmerv *aux_kmers; // translated ids
	u32list *ids; // ctg ids prepare to handle, those parse min_kmer
	uint32_t min_kmer; // parameter: # kmers to define two similar contigs
	uint32_t min_overlap; // parameter
	float het; // parameter
	uint32_t CTG_KMER_SIZE; // parameter
	uint32_t min_ol; //parameter for asm
	float min_sm; // parameter for asm
	uint32_t min_read; // parameter for asm
	uint32_t max_read; // parameter for asm
	uint32_t sim_pairs;
	EF *ef;
	int flag;  // if == 0 first use, init; else reset
} merge_t;

static inline int cmp_kmer_pos(const void *e1, const void *e2) {
	ctg_kmer_t *t1, *t2;
	t1 = (ctg_kmer_t *)e1;
	t2 = (ctg_kmer_t *)e2;

	if (t1->kpos == t2->kpos)
		return 0;
	else if (t1->kpos > t2->kpos)
		return 1;
	else
		return -1;
}

static inline int cmp_ids(const void *e1, const void *e2) {
	ctg_kmer_t *t1, *t2;
	t1 = (ctg_kmer_t *)e1;
	t2 = (ctg_kmer_t *)e2;

	if (t1->id == t2->id)
		return 0;
	else if (t1->id > t2->id)
		return 1;
	else
		return -1;
}
//
//static inline int cmp_ctg_clsids(const void *e1, const void *e2) {
//	contig_t *t1, *t2;
//	t1 = (contig_t *)e1;
//	t2 = (contig_t *)e2;
//
//	if (t1->clsid == t2->clsid)
//		return 0;
//	else if (t1->clsid > t2->clsid)
//		return 1;
//	else
//		return -1;
//
//}
//
#ifdef __CPLUSPLUS
extern "C" {
#endif

merge_t* init_merger(uint32_t min_kmer, uint32_t min_overlap, float het, uint32_t kmersize);
//void merge_ctgs(merge_t *merger, FileReader *asmd, FileReader *divd, FILE *out);
void merge_ctgs(merge_t *merger, FileReader *in, FILE *out);
void assign_best_ctg(merge_t *merger, contig_t *ctg);
void prepare_EF(merge_t *merger, FileReader *in, uint32_t lastcid);
void merge_along_tree(merge_t *merger, pathtree_t *tree);
void merge_core(merge_t *merger);
void index_ctgs(merge_t *merger);
void free_index(merge_t *merger);
void free_ctg(contig_t *ctg);
void free_ctgs(merge_t *merger);
void build_tree(merge_t *merger);
void update_ctg2merge(merge_t *merger);
void prepare_ctgs(merge_t *merger, uint32_t i, contig_t *ctg, uint64_t pos);
int is_onegroup(merge_t *merger, contig_t *c1, contig_t *c2);
int is_similar_enough(merge_t *merger, contig_t *c1, contig_t *c2);
void merge_2ctg(merge_t *merger, contig_t *ctg1, contig_t *ctg2);
void update_merger(merge_t *merger, contig_t *ctg1, contig_t *ctg2);
void gather_leaves(pathtree_t *tree, u32list *leaves);
void prefix_path(char *s1, char *s2, int n, char *pre);
void print_clusters(merge_t *merger, FILE *out);
void destroy_tree(pathtree_t *t);
void free_tree(merge_t *merger);
void reset_merger(merge_t *merger);
void free_merger(merge_t *merger);

#ifdef __CPLUSPLUS
}
#endif

#endif
