#ifndef __MERGECONTIG_H
#define __MERGECONTIG_H

#include <stdint.h>
#include "list.h"
#include "file_reader.h"
#include "stdaln.h"
#include "string.h"
#include "heap.h"
#include "hashset.h"

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
PWDB* pw_aln_contigs(CtgDB *db);
PWDB* clustering_ctg(PWDB *db, uint32_t overlap, float het);
void print_clusters(PWDB *db);
void execute_pwaln(CtgDB *db, uint32_t overlap, float het, uint32_t max_nctg);
void free_pwdb(PWDB *db);

#ifdef __CPLUSPLUS
}
#endif

#endif
