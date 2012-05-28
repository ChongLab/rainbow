/*
 * 
 * Copyright (c) 2012, Zechen Chong <chongzechen@gmail.com>
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

#include "mergectg.h"

typedef struct {uint32_t key; uint32_t oldid; char *path;} uuchash_t;
#define uuchash_code(e) (e).key
#define uuchash_equals(e1, e2) ((e1).key == (e2).key)
define_hashset(uuchash, uuchash_t, uuchash_code, uuchash_equals);

void merge_ctgs(merge_t *merger, FileReader *asmd, FileReader *divd, FILE *out) {
	uuchash_t h, *ph = NULL;
	uuchash *map = init_uuchash(1023);
	uint32_t key, oldid, i = 0; char *path = NULL, *prepath = NULL;
	contig_t *ctg;
	int exists;
	int len = 0;
	uint32_t id = 0, preid = 0, preoldid = 0;
	String *line = init_string(1);
	char *seq = NULL; 


	while (fread_table(divd) != -1) {
		key = atoi(get_col_str(divd, 1));
		oldid = atoi(get_col_str(divd, 4));
		path = get_col_str(divd, 5);
		h.key = key;
		h.oldid = oldid;
		ph = prepare_uuchash(map, h, &exists);
		if(exists){
		} else {
			h.path = strdup(path);
			*ph = h;
		}
	}

	while (fread_line(line, asmd) != -1) {
		if (line->string[0] == 'E') {
			id = atoi(line->string+2);
			h.key = id;
			h.oldid = 0;
			h.path = NULL;
			ph = get_uuchash(map, h);
			if (preoldid != ph->oldid) {
				if (preid) {  //first load last contig
					ctg = next_ref_contigv(merger->ctgs);
					ctg->id = preid;
					ctg->clsid = i++;
					ctg->old_clsid = preoldid;
					ctg->sz = 1;
					ctg->seq = seq;
					ctg->path = prepath;
					len = 0;
				}
				if (count_contigv(merger->ctgs) > 1) {
					merge_core(merger, out); //
				}
				if (count_contigv(merger->ctgs) >= 1) {
					for (i = 0; i < count_contigv(merger->ctgs); i++) { //free seq and paths
						seq = (ref_contigv(merger->ctgs, i))->seq;
						free(seq);
						seq = NULL;
					}
				}
				reset_merger(merger);
				i = 0;
				preid = id;
				preoldid = ph->oldid;
				prepath = ph->path;
			} else {
				ctg = next_ref_contigv(merger->ctgs);
				ctg->id = preid;
				ctg->clsid = i++;
				ctg->old_clsid = preoldid;
				ctg->sz = 1;
				ctg->seq = seq;
				ctg->path = prepath;
				preid = id;
				preoldid = ph->oldid;
				prepath = ph->path;
				len = 0;
				seq = NULL;
			}
		} else if (line->string[0] == 'S') {
			if (len < (int)strlen(line->string+2)) {
				len = (int)strlen(line->string+2);
				if (seq) free(seq); seq = NULL;
				seq = strdup(line->string+2);
			}
		}
	}
	if (preid) {  //first load last contig
		ctg = next_ref_contigv(merger->ctgs);
		ctg->id = preid;
		ctg->clsid = i++;
		ctg->old_clsid = preoldid;
		ctg->sz = 1;
		ctg->seq = seq;
		ctg->path = prepath;
		len = 0;
		//free(seq);
	}
	
	if (count_contigv(merger->ctgs) > 1) {
		merge_core(merger, out); // 
	}
	
	if (count_contigv(merger->ctgs) >= 1) {
		for (i = 0; i < count_contigv(merger->ctgs); i++) { //free seq and paths
			seq = (ref_contigv(merger->ctgs, i))->seq;
			free(seq);
			seq = NULL;
		}
	}
	reset_iter_uuchash(map);
	while (ph = ref_iter_uuchash(map)) {
		//fprintf(stdout, "%s\n", ph->path);
		//fflush(stdout);
		free(ph->path);
	}
	free_string(line);
	free_uuchash(map);
}

void merge_core(merge_t *merger, FILE *out) {
	contigv *ctgs = merger->ctgs;
	contig_t *ctg; uint32_t i;
	for (i = 0; i < count_contigv(ctgs); i++) {
		ctg = ref_contigv(ctgs, i);
		fprintf(stdout, "%u\t%u\t%u\t%u\t%s\t%s\n", ctg->id, ctg->clsid, ctg->old_clsid, ctg->sz, ctg->seq, ctg->path);
		fflush(stdout);
	}
	fprintf(stdout, "\n");
}

merge_t* init_merger(uint32_t min_kmer, uint32_t min_overlap, float het, uint32_t kmersize) {
	merge_t *merger;
	merger = (merge_t *)malloc(sizeof(merge_t));
	merger->ctgs = init_contigv(2);
	merger->tree = NULL;
	merger->index = init_ctgkhash(23);
	merger->links = NULL;
	merger->idv = NULL;
	merger->kmers = init_ctgkmerv(23);
	merger->aux_kmers = init_ctgkmerv(12);
	merger->min_kmer = min_kmer;
	merger->min_overlap = min_overlap;
	merger->het = het;
	merger->CTG_KMER_SIZE = kmersize;
	
	return merger;
}

void reset_merger(merge_t *merger) {
	clear_contigv(merger->ctgs);
	clear_ctgkhash(merger->index);
	clear_ctgkmerv(merger->kmers);
	clear_ctgkmerv(merger->aux_kmers);
}

void free_merger(merge_t *merger) {
	free_contigv(merger->ctgs);
	free(merger->tree);
	free_ctgkhash(merger->index);
	free(merger->links);
	free(merger->idv);
	free_ctgkmerv(merger->kmers);
	free_ctgkmerv(merger->aux_kmers);
	
	free(merger);
}
