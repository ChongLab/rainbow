#ifndef __MERGECONTIG_H
#define __MERGECONTIG_H

#include <stdint.h>
#include "list.h"
#include "file_reader.h"
#include "stdaln.h"
#include "string.h"
#include "heap.h"
#include "hashset.h"
#include "rainbow.h"

#define KMER_SIZE_CTG 15

/* char -> 17 (=16+1) nucleotides */
static unsigned char aln_nt16_table[256] = {
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,16 /*'-'*/,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15, 1,14, 4, 11,15,15, 2, 13,15,15,10, 15, 5,15,15,
        15,15, 3, 6,  8,15, 7, 9,  0,12,15,15, 15,15,15,15,
        15, 1,14, 4, 11,15,15, 2, 13,15,15,10, 15, 5,15,15,
        15,15, 3, 6,  8,15, 7, 9,  0,12,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};


typedef struct {
	uint64_t kmer;
	uint64_t kpos;
} kmer_tt;

#define kmer_code(k) u64hashcode((k).kmer)
#define kmer_eq(k1, k2) ((k1).kmer == (k2).kmer)
define_hashset(kmerhash, kmer_tt, kmer_code, kmer_eq);

typedef struct {
	uint64_t pos;
	uint32_t lastoffset;
	uint32_t offset;
} kmer_pos_t;

define_list(posv, kmer_pos_t);

static inline int cmp_kmer_pos (const void *e1, const void *e2) {
	kmer_pos_t *t1, *t2;
	t1 = (kmer_pos_t *)e1;
	t2 = (kmer_pos_t *)e2;

	if (t1->pos == t2->pos)
		return 0;
	else if (t1->pos > t2->pos)
		return 1;
	else
		return -1;
}

typedef struct {
	uint32_t id;
	uint32_t offset;
	uint32_t lastoffset;
} id_tt;

define_list(idlist, id_tt);

static inline int cmp_ids(const void *e1, const void *e2) {
	id_tt *t1, *t2;
	t1 = (id_tt *)e1;
	t2 = (id_tt *)e2;

	if (t1->id == t2->id)
		return 0;
	else if (t1->id > t2->id)
		return 1;
	else
		return -1;
}

static inline void aln_str(char *s1, char *s2, int *mm, int *mn, int *score) {
	int len, len1, len2, i, j, k, s;
	len1 = strlen(s1);
	len2 = strlen(s2);
	len = len1<len2?len1:len2;
	s = 0;
	int m = 0;
	for (i = 0; i < len; i++) {
		j = aln_nt16_table[(int)s1[i]];
		k = aln_nt16_table[(int)s2[i]];
		s += aln_sm_nt[j+k*16]; // magic aln_sm_nt table ROW_NUMBER=16
		if (s1[i] != s2[i])
			m++;
	}
	*mm = m;
	*mn = len;
	*score = s;
}

typedef struct {
	uint64_t last; //last kmer position
	uint32_t offset;
} link_t;

typedef struct {
	uint32_t id;
	uint32_t cls_id;
	uint32_t old_clsid;
	uint32_t sz;   //union tree depth
	char *seq;
} Ctg;

define_list(ctglist, Ctg);

typedef struct {
	uint32_t ctgnum;
	ctglist *ctgs;
} CtgDB;

typedef struct {
	uint32_t id0;
	uint32_t id1;
	uint32_t overlap;
	float het;
	int score; 
} PWcontig;

define_list(pwctglist, PWcontig*);

typedef struct {
	pwctglist *pwctgs;
	Heap *hp;
	ctglist *ctgv;
} PWDB;

#ifdef __CPLUSPLUS
extern "C" {
#endif

CtgDB* init_ctgdb(void );
CtgDB* load_ctgdb(FileReader *fr1, FileReader *fr2);
void print_ctgdb(CtgDB *db);
void free_ctgdb(CtgDB *db);
void free_load_ctgdb(CtgDB *db);
PWDB* pw_aln_contigs(CtgDB *db, uint32_t overlap, float het);
PWDB* pw_aln_contigs_brute(CtgDB *db);
PWDB* clustering_ctg(PWDB *db, uint32_t overlap, float het);
void print_clusters(PWDB *db);
void execute_pwaln(CtgDB *db, uint32_t overlap, float het, uint32_t max_nctg);
void free_pwdb(PWDB *db);

#ifdef __CPLUSPLUS
}
#endif

#endif
