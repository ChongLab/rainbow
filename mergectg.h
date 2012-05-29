#ifndef MERGECTG_H
#define MERGECTG_H

#include <stdint.h>
#include "list.h"
#include "file_reader.h"
#include "hashset.h"
#include "string.h"

typedef struct {
	uint32_t id, clsid, old_clsid, sz;
	char *seq, *path;
} contig_t; 

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

define_list(ctgkmerv, ctg_kmer_t);

typedef struct {
	uint64_t last; //last kmer position
	int offset; // current kmer offset
} link_t;

typedef struct {
	contigv *ctgs;
	pathtree_t *tree;
	ctgkhash *index;
	link_t *links;
	uint64_t *idv;  // to search translate pos to ids
	ctgkmerv *kmers;   // all kmers of each kmer of query contig
	ctgkmerv *aux_kmers; // translated ids
	uint32_t min_kmer; // parameter: # kmers to define two similar contigs
	uint32_t min_overlap; // parameter
	float het; // parameter
	uint32_t CTG_KMER_SIZE; // parameter
} merge_t;

#ifdef __CPLUSPLUS
extern "C" {
#endif

merge_t* init_merger(uint32_t min_kmer, uint32_t min_overlap, float het, uint32_t kmersize);
void merge_ctgs(merge_t *merger, FileReader *asmd, FileReader *divd, FILE *out);
void merge_core(merge_t *merger, FILE *out);
void index_ctgs(merge_t *merger);
void free_index(merge_t *merger);
void build_tree(merge_t *merger);
void free_tree(merge_t *merger);
void reset_merger(merge_t *merger);
void free_merger(merge_t *merger);

#ifdef __CPLUSPLUS
}
#endif

#endif
