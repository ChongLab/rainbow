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
#include "rainbow.h"

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
					index_ctgs(merger);
					build_tree(merger);
					merge_core(merger); //TODO
					print_clusters(merger, out);
					free_index(merger);
					free_tree(merger);
				}
				if (count_contigv(merger->ctgs) >= 1) {
					for (i = 0; i < count_contigv(merger->ctgs); i++) { //free seq
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
		index_ctgs(merger);
		build_tree(merger);
		merge_core(merger); // TODO
		print_clusters(merger, out);
		free_index(merger);
		free_tree(merger);
	}
	
	if (count_contigv(merger->ctgs) >= 1) {
		for (i = 0; i < count_contigv(merger->ctgs); i++) { //free seq 
			seq = (ref_contigv(merger->ctgs, i))->seq;
			free(seq);
			seq = NULL;
		}
	}
	reset_iter_uuchash(map);
	while ((ph = ref_iter_uuchash(map))) {
		free(ph->path);
	}
	free_string(line);
	free_uuchash(map);
}

void prepare_ctgs(merge_t *merger, uint32_t i, contig_t *ctg, uint64_t pos) {
	clear_ctgkmerv(merger->kmers);
	clear_ctgkmerv(merger->aux_kmers);
	clear_u32list(merger->ids);
	uint32_t j, n;
	int seqlen, idx, count, pre;
	ctg_kmer_t K, *t, *tpos;
	uint64_t kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-merger->CTG_KMER_SIZE)*2), next, bt, p;

	seqlen = strlen(ctg->seq);
	for (j = 0; j < merger->CTG_KMER_SIZE-1; j++)
		K.kmer = (K.kmer << 2) | base_bit_table[(int)ctg->seq[j]];
	for (j = 0; j <= (unsigned)seqlen-merger->CTG_KMER_SIZE; j++) {
		K.kmer = ((K.kmer << 2) | base_bit_table[(int)ctg->seq[j+merger->CTG_KMER_SIZE-1]]) & kmask;
		t = get_ctgkhash(merger->index, K);
		bt = t->kpos;
		if (bt >= pos+seqlen || bt < pos) {
			tpos = next_ref_ctgkmerv(merger->kmers);
			tpos->kmer = K.kmer;
			tpos->kpos = bt;
			tpos->id = 0;  
			tpos->offset = merger->links[bt].offset;
			tpos->offset2 = j;
		}
		while (1) { //tracing positions
			next = merger->links[bt].last;
			if (next == bt) break;
			if (next >= pos+seqlen || next < pos) { // omit the contig and those before it
				tpos = next_ref_ctgkmerv(merger->kmers);
				tpos->kmer = K.kmer;
				tpos->kpos = bt;
				tpos->id = 0; 
				tpos->offset = merger->links[bt].offset;
				tpos->offset2 = j;
			}
			bt = next;
		}
	}
	qsort(as_array_ctgkmerv(merger->kmers), count_ctgkmerv(merger->kmers), sizeof(ctg_kmer_t), cmp_kmer_pos);

	//translate positions into ids and jump overlap kmers
	n = count_contigv(merger->ctgs);
	t = ref_ctgkmerv(merger->kmers, 0);
	p = t->kpos;
	idx = bisearch(merger->idv, n, p, NULL);
	if (idx < 0) idx = -1 - idx;
	if (p <= merger->idv[idx]) {
		if (idx > (int)i) {
			tpos = next_ref_ctgkmerv(merger->aux_kmers);
			*tpos = *t;
			tpos->id  = idx;
		}
	}
	for (j = 1; j < count_ctgkmerv(merger->kmers); j++) {
		tpos = ref_ctgkmerv(merger->kmers, j);
		if (tpos->kpos - p >= merger->CTG_KMER_SIZE) {
			t = tpos;
			p = t->kpos;
		} else {
			continue;
		}

		if (p <= merger->idv[idx]) {
			if (idx > (int)i) {
				tpos = next_ref_ctgkmerv(merger->aux_kmers);
				*tpos = *t;
				tpos->id  = idx;
			}
		} else {
			idx = bisearch(merger->idv, n, p, NULL);
			if (idx < 0) idx = -1 - idx;
			if (idx > (int)i) {
				tpos = next_ref_ctgkmerv(merger->aux_kmers);
				*tpos = *t;
				tpos->id = idx;
			}
		}
	}

	qsort(as_array_ctgkmerv(merger->aux_kmers), count_ctgkmerv(merger->aux_kmers), sizeof(ctg_kmer_t), cmp_ids);

	count = 0;
	pre = -1;
	for (j = 0; j < count_ctgkmerv(merger->aux_kmers); j++) {
		tpos = ref_ctgkmerv(merger->aux_kmers, j);
		if (pre != (int)tpos->id) {
			if (count >= (int)merger->min_kmer) {  // parameter min_kmer used here
				push_u32list(merger->ids, t->id);
			}
			count = 0;
			t = tpos;
			count++;
			pre = tpos->id;
		} else {
			count++;
			t = tpos;
		}
		//printf("%u\t%llu\t%llu\t%d\t%d\t%d\n", i, tpos->kpos,tpos->kmer, tpos->offset, tpos->offset2, tpos->id);
	}
	if (count >= (int)merger->min_kmer) {// parameter min_kmer used here
		push_u32list(merger->ids, t->id);
	} /*
	for (j = 0; j < count_u32list(merger->ids); j++) {
		printf("%d ", get_u32list(merger->ids, j));
	}
	printf("\n"); */
	clear_ctgkmerv(merger->kmers);
	clear_ctgkmerv(merger->aux_kmers);
	return;
}

int is_onegroup(merge_t *merger, contig_t *c1, contig_t *c2) {
	uint32_t i, j;
	contigv *ctgs = merger->ctgs;
	for (i = c1->clsid; i != (ref_contigv(ctgs, i))->clsid; i = (ref_contigv(ctgs, i))->clsid) ;
	for (j = c2->clsid; j != (ref_contigv(ctgs, j))->clsid; j = (ref_contigv(ctgs, j))->clsid) ;
	
	if (i == j) return 1;
	return 0;
}

int is_similar_enough(merge_t *merger, contig_t *c1, contig_t *c2) {
	int mm, mn, aln_len, i;
	AlnAln *aa;
	AlnParam ap = {10, 2, 2, aln_sm_nt, 16, 75};

	mn = mm = 0;
	aa = aln_stdaln(c1->seq, c2->seq, &ap, 0, 1);
	
	aln_len = strlen(aa->out1);
	for (i = 0; i <aln_len; i++) {
		if (aa->out1[i] == '-' || aa->out2[i] == '-')
			continue;
		if (aa->out1[i] != aa->out2[i])
			mm++;
		mn++;
	}
	
	aln_free_AlnAln(aa);
	if (mn > (int)merger->min_overlap && ((float)mm/mn - merger->het) <= 0)
		return 1;

	return 0;
}

void merge_ctg(merge_t *merger, contig_t *c1, contig_t *c2) {
	uint32_t i, j;
	contigv *ctgs = merger->ctgs;
	for (i = c1->clsid; i != (ref_contigv(ctgs, i))->clsid; i = (ref_contigv(ctgs, i))->clsid)
		ref_contigv(ctgs, i)->clsid = ref_contigv(ctgs, ref_contigv(ctgs, i)->clsid)->clsid;
	for (j = c2->clsid; j != (ref_contigv(ctgs, j))->clsid; j = (ref_contigv(ctgs, j))->clsid) ;
		ref_contigv(ctgs, j)->clsid = ref_contigv(ctgs, ref_contigv(ctgs, j)->clsid)->clsid;
	
	if (i == j) return;
	
	if (ref_contigv(ctgs, i)->sz < ref_contigv(ctgs, j)->sz) {
		ref_contigv(ctgs, i)->clsid = j;
		ref_contigv(ctgs, j)->sz += ref_contigv(ctgs, i)->sz;
	} else {
		ref_contigv(ctgs, j)->clsid = i;
		ref_contigv(ctgs, i)->sz += ref_contigv(ctgs, j)->sz;
	}

	return;
}

void prefix_path(char *s1, char *s2, int n, char *pre) {
	int i;
	for (i = 0; i < n; i++) {
		if (s1[i] == s2[i]) {
			pre[i] = s1[i];
		} else {
			break;
		}
	}
	pre[i] = '\0';

	return;
}

void gather_leaves(pathtree_t *tree, u32list *leaves) {
	if (tree == NULL)
		return;
	
	if (tree->left == NULL && tree->right == NULL) {
		push_u32list(leaves, tree->tid);
	} else {
		gather_leaves(tree->left, leaves);
		gather_leaves(tree->right, leaves);
	}

}

void update_merger(merge_t *merger, contig_t *c1, contig_t *c2) {
	u32list *subtree;
	pathtree_t *t; int i, n;
	char *prefix;
	n = strlen(c1->path)>=strlen(c2->path)?strlen(c2->path):strlen(c1->path);
	prefix = (char*)malloc(sizeof(char)*(n+1));
	memset(prefix, 0, n+1);
	prefix_path(c1->path, c2->path, n, prefix);
	if (strcmp(prefix, "") == 0) {
		free(prefix);
		return;
	}
	subtree = init_u32list(2);
	n = strlen(prefix);
	t = merger->tree;
	for (i = 0; i < n; i++) {
		if (prefix[i] == '0') {
			t = t->left;
		} else if (prefix[i] == '1') {
			t = t->right;
		}
		if (t == NULL) {
			free(prefix);
			free_u32list(subtree);
			return;
		}
	}
	gather_leaves(t, subtree);
	n = count_u32list(subtree);
	for (i = 0; i < n; i++) {
		merge_ctg(merger, c1, ref_contigv(merger->ctgs, get_u32list(subtree, i)));
	}
	free(prefix);
	free_u32list(subtree);
	return;
}

void merge_core(merge_t *merger) {
	contigv *ctgs = merger->ctgs;
	contig_t *ctg1, *ctg2; uint32_t i, j, ii;
	int seqlen;
	uint64_t pos = 0;
	for (i = 0; i < count_contigv(ctgs)-1; i++) {
		ctg1 = ref_contigv(ctgs, i);
		seqlen = strlen(ctg1->seq);
		if (seqlen < (int)merger->CTG_KMER_SIZE) continue;
		prepare_ctgs(merger, i, ctg1, pos);
		for (ii = 0; ii < count_u32list(merger->ids); ii++) {
			j = get_u32list(merger->ids, ii);
			ctg2 = ref_contigv(ctgs, j);
			if (is_onegroup(merger, ctg1, ctg2)) {
				continue;
			} else {
				if (is_similar_enough(merger, ctg1, ctg2)) {
					merge_ctg(merger, ctg1, ctg2);
					update_merger(merger, ctg1, ctg2);
				}
			}
		}
		//fprintf(stdout, "id\tclsid\told_clsid\tsz\tseq\tpath\n");
		//fprintf(stdout, "%u\t%u\t%u\t%u\t%s\t%s\n", ctg1->id, ctg1->clsid, ctg1->old_clsid, ctg1->sz, ctg1->seq, ctg1->path);
		//fflush(stdout);
		pos += seqlen;
	}
	//fprintf(stdout, "\n");
	fflush(stdout);
}

void build_tree(merge_t *merger) {
	uint32_t i, n;
	n = count_contigv(merger->ctgs);
	char *path; int len = 0, j;
	contig_t *ctg;
	pathtree_t *t;
	if (merger->tree == NULL) {
		merger->tree = (pathtree_t *) malloc(sizeof(pathtree_t));
		merger->tree->left = NULL;
		merger->tree->right = NULL;
		merger->tree->tid = 0;
	}
	t = merger->tree;
	for (i = 0; i < n; i++) {
		ctg = ref_contigv(merger->ctgs, i);
		path = ctg->path;
		len = strlen(path);
		for (j = 0; j < len; j++) {
			if (path[j] == '0') {
				if (t->left == NULL) {
					t->left = (pathtree_t *) malloc(sizeof(pathtree_t));
					t->left->left = NULL;
					t->left->right = NULL;
					t->left->tid = 0;
					if (j == len-1)
						t->left->tid = i;
				}
				t = t->left;
			} else { // '1'
				if (t->right == NULL) {
					t->right = (pathtree_t *) malloc(sizeof(pathtree_t));
					t->right->left = NULL;
					t->right->right = NULL;
					t->right->tid = 0;
					if (j == len-1)
						t->right->tid = i;
				}
				t = t->right;
			}
		}
		t = merger->tree;
	}
	return;
}
void print_clusters(merge_t *merger, FILE *out) {
	int last_cid = -1;
	int line_num = 0;
	uint32_t i, n;
	contig_t *ctg;
	contigv *ctgs = merger->ctgs;
	
	qsort(as_array_contigv(ctgs), count_contigv(ctgs), sizeof(contig_t), cmp_ctg_clsids);
	n = count_contigv(ctgs);
	for (i = 0; i < n; i++) {
		ctg = ref_contigv(ctgs, i);
		if (last_cid != (int)ctg->clsid && line_num > 0)
			fprintf(out, "\n");
		last_cid = ctg->clsid;
		line_num++;
		fprintf(out, "%d ", ctg->id);
	}
	fprintf(out, "\n");
	fflush(out);
	return;

}

void destroy_tree(pathtree_t *t) {
	if (t != NULL) {
		destroy_tree(t->left);
		destroy_tree(t->right);
		free(t);
		t = NULL;
	}
}

void free_tree(merge_t *merger) {
	destroy_tree(merger->tree);
	merger->tree = NULL;
	return;
}

void index_ctgs(merge_t *merger) {
	uint64_t kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-merger->CTG_KMER_SIZE)*2), pos;
	uint32_t n, i, j; int seqlen, exists;
	contig_t *ctg;
	merger->index = init_ctgkhash(23);
	ctg_kmer_t K, *t;
	n = count_contigv(merger->ctgs);
	merger->idv = (uint64_t *) malloc(n * sizeof(uint64_t));
	link_t *tmp;

	K.kmer = 0;
	K.kpos = 0;
	K.id = 0;
	K.offset = -1;
	K.offset2 = -1;
	pos = 0;
	for (i = 0; i < n; i++) {
		ctg = ref_contigv(merger->ctgs, i);
		seqlen = strlen(ctg->seq);
		if (seqlen < (int)merger->CTG_KMER_SIZE) continue;
		tmp = (link_t *)realloc(merger->links, (pos+seqlen)*sizeof(link_t));
		if (tmp == NULL) {
			free(merger->links);
			fprintf(stderr, "Memory allocation error!!!\n");
			abort();
		}
		merger->links = tmp;
		merger->idv[i] = pos+seqlen-1;
		
		for (j = 0; j < merger->CTG_KMER_SIZE-1; j++)
			K.kmer = (K.kmer << 2) | base_bit_table[(int)ctg->seq[j]];
		for (j = 0; j <= (unsigned)seqlen-merger->CTG_KMER_SIZE; j++) {
			K.kmer = ((K.kmer << 2) | base_bit_table[(int)ctg->seq[j+merger->CTG_KMER_SIZE-1]]) & kmask;
			t = prepare_ctgkhash(merger->index, K, &exists);
			if (exists) {
				merger->links[pos+j].last = t->kpos;
				merger->links[pos+j].offset = j;
			} else {
				t->kmer = K.kmer;
				merger->links[pos+j].last = pos+j;
				merger->links[pos+j].offset = j;
			}
			t->kpos = pos+j;
		}
		pos += seqlen;
	}
}

merge_t* init_merger(uint32_t min_kmer, uint32_t min_overlap, float het, uint32_t kmersize) {
	merge_t *merger;
	merger = (merge_t *)malloc(sizeof(merge_t));
	merger->ctgs = init_contigv(2);
	merger->tree = NULL;
	merger->index = NULL;
	merger->links = NULL;
	merger->idv = NULL;
	merger->kmers = init_ctgkmerv(23);
	merger->aux_kmers = init_ctgkmerv(12);
	merger->ids = init_u32list(2);
	merger->min_kmer = min_kmer;
	merger->min_overlap = min_overlap;
	merger->het = het;
	merger->CTG_KMER_SIZE = kmersize;
	
	return merger;
}

void reset_merger(merge_t *merger) {
	clear_contigv(merger->ctgs);
	//clear_ctgkmerv(merger->kmers);
	//clear_ctgkmerv(merger->aux_kmers);
}

void free_index(merge_t *merger) {
	free_ctgkhash(merger->index);
	free(merger->links);
	free(merger->idv);
	merger->links = NULL;
	merger->idv = NULL;
}

void free_merger(merge_t *merger) {
	free_contigv(merger->ctgs);
	//free(merger->tree);
	free_ctgkmerv(merger->kmers);
	free_ctgkmerv(merger->aux_kmers);
	free_u32list(merger->ids);
	
	free(merger);
}
