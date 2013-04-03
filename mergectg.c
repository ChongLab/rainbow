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
#include "vector.h"

contig_t* lend_ctgv_merger(merge_t *merger) {
	contig_t *ctg;
	if (!pop_contigsv(merger->cache, &ctg)) {
//		ctg = next_ref_contigv(merger->ctgs);
		ctg = malloc(sizeof(contig_t));
		ctg->path = init_string(10);
		ctg->rds = init_readv(50);
		ctg->index = init_rdkhash(24);
		ctg->m_idx = init_idxv(2);
		ctg->m_rds = init_u32list(2);
	}
	return ctg;
}

void return_ctgv_merger(merge_t *merger, contig_t *ctg) {
	clear_ctg(ctg);
	push_contigsv(merger->cache, ctg);
}

void prepare_reads(merge_t *merger, FileReader *in, uint32_t lastcid) {
	int n_col;
	uint32_t cid, ef_id, eid, eflen, id = 0;
	ef_id = 0; 
	contig_t *ctg = NULL;
	read_t *rd;
	char *efstr, *path;

	while ((n_col = fread_table(in)) != -1) {
		if (n_col == 0) continue;
		cid = atoi(get_col_str(in, 4));
		if (cid != lastcid) {
			froll_back(in);
			return;
		}
		eid = atoi(get_col_str(in, 1));
		if (eid != ef_id) {
			ef_id = eid;
			efstr = get_col_str(in, 2);
			eflen = get_col_len(in, 2);
			reverse_dna(efstr, eflen);
			ctg = lend_ctgv_merger(merger); 
			push_contigv(merger->ctgs, ctg);
			ctg->id = ef_id;
			ctg->closed = 0;
			path = get_col_str(in, 5);
			append_string(ctg->path, path, strlen(path));
//			ctg->path = strdup(get_col_str(in, 5));
//			ctg->efctgs = init_vec(sizeof(FContig*), 6);
			push_u32list(ctg->m_rds, id++);
			rd = next_ref_readv(ctg->rds);
			rd->seq_id = atol(get_col_str(in, 0));
			rd->rd_len = eflen;
			memmove(rd->seq, efstr, rd->rd_len);
			rd->seq[rd->rd_len]= '\0';
			rd->rank = 1;

		}
		rd = next_ref_readv(ctg->rds);
		rd->rank = 1;
		rd->seq_id = atol(get_col_str(in, 0));
		rd->rd_len = get_col_len(in, 3);
		memmove(rd->seq, get_col_str(in, 3), rd->rd_len);
		rd->seq[rd->rd_len]= '\0';
	}
}

void print_asm(merge_t *merger, FILE *out) {
	uint32_t i, j, k;
	contig_t *ctg, *ctg2;
	read_t *rd;

	for (k = 0; k < merger->ctgs->size; k++) {
		ctg = get_contigv(merger->ctgs, k);
		if (merger->flag && !ctg->closed && ctg->rds->size >= merger->min_read && ctg->rds->size <= merger->max_read) {  // have used ef
			rd = ref_readv(ctg->rds, 0);
			reset_ef(merger->ef, ctg->id, rd->seq, rd->rd_len, merger->min_ol, merger->min_sm);
			for (i = 1; i < ctg->rds->size; i++) {
				rd = ref_readv(ctg->rds, i);
				add_read2ef(merger->ef, rd->seq, rd->seq_id, rd->rd_len, rd->rank);
			} 
			for (j = 1; j < ctg->m_rds->size; j++) {
				ctg2 = get_contigv(merger->ctgs, get_u32list(ctg->m_rds, j));
				for (i = 0; i < ctg2->rds->size; i++) {
					rd = ref_readv(ctg2->rds, i);
					add_read2ef(merger->ef, rd->seq, rd->seq_id, rd->rd_len, rd->rank);
				}
			}
			align_reads_ef(merger->ef);
			asm_ef_ctgs(merger->ef);
			output_ef_ctgs(merger->ef, out);
			
		}
		if (!merger->flag && !ctg->closed && ctg->rds->size >= merger->min_read && ctg->rds->size <= merger->max_read) {
			merger->flag = 1;
			rd = ref_readv(ctg->rds, 0);
			merger->ef = init_ef(ctg->id, rd->seq, rd->rd_len, merger->min_ol, merger->min_sm);
			for (i = 1; i < ctg->rds->size; i++) {
				rd = ref_readv(ctg->rds, i);
				add_read2ef(merger->ef, rd->seq, rd->seq_id, rd->rd_len, rd->rank);
			}
			for (j = 1; j < ctg->m_rds->size; j++) {
				ctg2 = get_contigv(merger->ctgs, get_u32list(ctg->m_rds, j));
				for (i = 0; i < ctg2->rds->size; i++) {
					rd = ref_readv(ctg2->rds, i);
					add_read2ef(merger->ef, rd->seq, rd->seq_id, rd->rd_len, rd->rank);
				}
			}
			align_reads_ef(merger->ef);
			asm_ef_ctgs(merger->ef);
			output_ef_ctgs(merger->ef, out);

		}
	}
}

void print_asm2(merge_t *merger, FILE *out) {
	uint32_t i, j, k;
	contig_t *ctg, *ctg2;
	read_t *rd;

	for (k = 0; k < merger->ctgs->size; k++) {
		ctg = get_contigv(merger->ctgs, k);
		if (merger->flag && !ctg->closed && ctg->rds->size >= merger->min_read && ctg->rds->size <= merger->max_read) {  // have used ef
			for (i = 1; i < ctg->rds->size; i++) {
				rd = ref_readv(ctg->rds, i);
				fprintf(out, "%u\t%u\n", merger->cid, rd->seq_id);
				fflush(out);
			} 
			for (j = 1; j < ctg->m_rds->size; j++) {
				ctg2 = get_contigv(merger->ctgs, get_u32list(ctg->m_rds, j));
				for (i = 0; i < ctg2->rds->size; i++) {
					rd = ref_readv(ctg2->rds, i);
					fprintf(out, "%u\t%u\n", merger->cid, rd->seq_id);
					fflush(out);
				}
			}
			merger->cid++;
		}
		if (!merger->flag && !ctg->closed && ctg->rds->size >= merger->min_read && ctg->rds->size <= merger->max_read) {
			merger->flag = 1;
			rd = ref_readv(ctg->rds, 0);
			for (i = 1; i < ctg->rds->size; i++) {
				rd = ref_readv(ctg->rds, i);
				fprintf(out, "%u\t%u\n", merger->cid, rd->seq_id);
				fflush(out);
			}
			for (j = 1; j < ctg->m_rds->size; j++) {
				ctg2 = get_contigv(merger->ctgs, get_u32list(ctg->m_rds, j));
				for (i = 0; i < ctg2->rds->size; i++) {
					rd = ref_readv(ctg2->rds, i);
					fprintf(out, "%u\t%u\n", merger->cid, rd->seq_id);
					fflush(out);
				}
			}
			merger->cid++;
		}
	}
}

void index_rds(merge_t *merger, contig_t *ctg) {
	uint64_t kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-merger->RD_KMER_SIZE)*2), pos;
	read_t *rd; uint32_t i, j, len; int exists;
	rd_kmer_t K, *t;

	K.kmer = 0;
	K.kpos = 0;
	pos = 0;
	for (i = 1; i < ctg->rds->size; i++) { 
		rd = ref_readv(ctg->rds, i);
		len = rd->rd_len;
		if (len < merger->RD_KMER_SIZE) continue;
		for (j = 0; j < merger->RD_KMER_SIZE-1; j++)
			K.kmer = (K.kmer << 2) | base_bit_table[(int)rd->seq[j]];
		for (j = 0; j <= (unsigned)len-merger->RD_KMER_SIZE; j++) {
			K.kmer = ((K.kmer << 2) | base_bit_table[(int)rd->seq[j+merger->RD_KMER_SIZE-1]]) & kmask;
			t = prepare_rdkhash(ctg->index, K, &exists);
			if (exists) {
			} else {
				t->kmer = K.kmer;
				t->kpos = 1;
			}
		}
		pos += len;
	}
}

void update_ctg2merge(merge_t *merger) {
	uint32_t j;
	contig_t *ctg;

	for (j = 0; j < merger->ctgs->size; j++) {
		ctg = get_contigv(merger->ctgs, j);
		index_rds(merger, ctg);
		push_idxv(ctg->m_idx, ctg->index);
	}
}
/*
static void print_leaf(merge_t *merger, pathtree_t *tree, FILE *out) {
	contig_t *ctg;
	if (tree->tid) {
		ctg = get_contigv(merger->ctgs, tree->tid-1);
		fprintf(out, ">%d\n", ctg->id);
		fflush(out);
	} else {
		print_leaf(merger, tree->left, out);
		print_leaf(merger, tree->right, out); 
	}
}

static void printnode(int c, int h) {
	int i;
	for (i = 0; i < h; i++)
		printf(" ");
	printf("%d\n", c);
}

static void show(pathtree_t *tree, int h) {
	if (tree == NULL) 
	{ printnode('*', h); return;}
	show(tree->right, h+1);
	printnode(tree->tid, h);
	show(tree->left, h+1);
}
*/

void merge_leaves(merge_t *merger, uint32_t id1, uint32_t id2) {
	char *prefix;
	uint32_t i; int n, n1, n2;
	contig_t *c1, *c2, *c;
	
	c1 = get_contigv(merger->ctgs, id1);
	c2 = get_contigv(merger->ctgs, id2);
	c = id1<id2?c1:c2;
	
	if (c == c1) {
		c2->closed = 1;
	} else {
		c1->closed = 1;
	}
	n1 = strlen(c1->path->string);
	n2 = strlen(c2->path->string);
	n = n1>=n2?n1:n2;
	prefix = (char*)malloc(sizeof(char)*(n+1));
	memset(prefix, 0, n+1);
	prefix_path(c1->path->string, c2->path->string, n, prefix);
	
	if (c2->closed) {
		for (i = 0; i < c2->m_rds->size; i++)
			push_u32list(c->m_rds, get_u32list(c2->m_rds, i));
	}
	else {
		for (i = 0; i < c1->m_rds->size; i++)
			push_u32list(c->m_rds, get_u32list(c1->m_rds, i));
	}
	clear_string(c->path);
	append_string(c->path, prefix, n);
	if (c == c1) {
		for (i = 0; i < c2->m_idx->size; i++)
			push_idxv(c->m_idx, get_idxv(c2->m_idx, i));
	} else {
		for (i = 0; i < c1->m_idx->size; i++)
			push_idxv(c->m_idx, get_idxv(c1->m_idx, i));
	}
	//	align_reads_ef(merger->ef);
//	asm_ef_ctgs(merger->ef);
//	output_ef_ctgs(merger->ef, stderr);
//	assign_best_ctg(merger, c);
	free(prefix);
} 

int is_similar_enough(merge_t *merger, contig_t *c1, contig_t *c2) {
	uint64_t kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-merger->RD_KMER_SIZE)*2);
	read_t *rd; uint32_t i, j, k, len, m, n, n1, n2, cnt; contig_t *cdb, *cq, *c; int found;
	rd_kmer_t K, *t;
	
	n1 = n2 = 0;
	for (m = 0; m < c1->m_rds->size; m++) {
		c = get_contigv(merger->ctgs, get_u32list(c1->m_rds, m));
		n1 += c->rds->size;
	}
	for (m = 0; m < c2->m_rds->size; m++) {
		c = get_contigv(merger->ctgs, get_u32list(c2->m_rds, m));
		n2 += c->rds->size;
	}
	if (n1 >= n2) {
		cdb = c1;
		cq = c2;
		n = n2-1;
	} else {
		cdb = c2;
		cq = c1;
		n = n1-1;
	}

	cnt = 0;
	for (m = 0; m < cq->m_rds->size; m++) {
		c = get_contigv(merger->ctgs, get_u32list(cq->m_rds,m));
		K.kmer = 0;
		K.kpos = 0;
		found = 0;
		if (m > 0)
			i = 0;
		else
			i = 1;
		for (; i < c->rds->size; i++) { 
			rd = ref_readv(c->rds, i);
			len = rd->rd_len;
			if (len < merger->RD_KMER_SIZE) continue;
			for (j = 0; j < merger->RD_KMER_SIZE-1; j++)
				K.kmer = (K.kmer << 2) | base_bit_table[(int)rd->seq[j]];
			for (j = 0; j <= (unsigned)len-merger->RD_KMER_SIZE-3; j++) {
				K.kmer = ((K.kmer << 2) | base_bit_table[(int)rd->seq[j+merger->RD_KMER_SIZE-1]]) & kmask;
				j++;
				K.kmer = ((K.kmer << 2) | base_bit_table[(int)rd->seq[j+merger->RD_KMER_SIZE-1]]) & kmask;
				j++;
				K.kmer = ((K.kmer << 2) | base_bit_table[(int)rd->seq[j+merger->RD_KMER_SIZE-1]]) & kmask;
				j++;
				K.kmer = ((K.kmer << 2) | base_bit_table[(int)rd->seq[j+merger->RD_KMER_SIZE-1]]) & kmask;
				for (k = 0; k < cdb->m_idx->size; k++) {
					t = get_rdkhash(get_idxv(cdb->m_idx, k), K);
					if (t) {
						cnt++;
						found = 1;
						break;
					}
				}
				if (found)
					break;
			}
			found = 0;
		}
	}
//	fprintf(stderr, "cnt=%d n=%d div=%f, n1=%d, n2=%d\n", cnt, n, (float)cnt/n, n1, n2);
	if ((float)cnt/n - merger->het >= 0)
		return 1;

	return 0;
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

void merge_core(merge_t *merger) {
	contigv *ctgs = merger->ctgs;
	contig_t *ctg1, *ctg2; uint32_t i, j;
	
	for (i = 0; i < ctgs->size-1; i++) {
		ctg1 = get_contigv(ctgs, i);
		if (ctg1->closed) continue;
		for (j = i+1; j < ctgs->size; j++) {
			if (ctg1->closed) break;
			ctg2 = get_contigv(ctgs, j);
			if (ctg2->closed || ctg2->rds->size <= 5 ) continue;
			if (is_similar_enough(merger, ctg1, ctg2)) {
				merge_leaves(merger, i, j);
				break;
			}
		}
	}
}

void build_tree(merge_t *merger) {
	uint32_t i, n;
	n = count_contigv(merger->ctgs);
//	reset_iter_ctgset(merger->ctgs);
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
		ctg = get_contigv(merger->ctgs, i);	
		path = ctg->path->string;
		len = strlen(path);
		for (j = 0; j < len; j++) {
			if (path[j] == '0') {
				if (t->left == NULL) {
					t->left = (pathtree_t *) malloc(sizeof(pathtree_t));
					t->left->left = NULL;
					t->left->right = NULL;
					t->left->tid = 0;
					if (j == len-1)
						t->left->tid = i+1;
				}
				t = t->left;
			} else { // '1'
				if (t->right == NULL) {
					t->right = (pathtree_t *) malloc(sizeof(pathtree_t));
					t->right->left = NULL;
					t->right->right = NULL;
					t->right->tid = 0;
					if (j == len-1)
						t->right->tid = i+1;
				}
				t = t->right;
			}
		}
		t = merger->tree;
	}
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

merge_t* init_merger(uint32_t min_kmer, uint32_t min_overlap, float het, uint32_t kmersize, uint32_t max_cluster, uint32_t need_asm, float min_sm, uint32_t min_read, uint32_t max_read) {
	merge_t *merger;
//	uint32_t skmer;
	merger = (merge_t *)malloc(sizeof(merge_t));
	merger->ctgs = init_contigv(2);
	merger->cache = init_contigsv(2);
//	merger->ctgs = init_ctgset(2);
	merger->tree = NULL;
	merger->min_kmer = min_kmer;
	merger->min_overlap = min_overlap;
	merger->het = het;
	merger->RD_KMER_SIZE = kmersize;
//	merger->min_ol = 5;              
	merger->min_ol = min_overlap;              
	merger->min_sm = min_sm;
	merger->min_read = min_read;
	merger->max_read = max_read;
	merger->need_asm = need_asm;
	merger->max_cluster = max_cluster;
	merger->sim_pairs = 0;
	merger->ef = NULL;
	merger->flag = 0;
	merger->cid = 0;
//	skmer = kmersize;
	return merger;
}

void free_ctg(contig_t *ctg) {
	free_string(ctg->path);
	free_u32list(ctg->m_rds);
	free_rdkhash(ctg->index);
	free_idxv(ctg->m_idx);
	free_readv(ctg->rds);
}

void clear_ctg(contig_t *ctg) {
	clear_string(ctg->path);
	clear_u32list(ctg->m_rds);
	clear_rdkhash(ctg->index);
	clear_idxv(ctg->m_idx);
	if (ctg->rds->size) clear_readv(ctg->rds);
}

void free_ctgs(merge_t *merger) {
	contig_t *ctg;
	uint32_t i;

	for (i = 0; i < merger->cache->size; i++) {
		ctg = get_contigsv(merger->cache, i);
		free_ctg(ctg);
		free(ctg);
	}
}

void put_cache_ctgs(merge_t *merger, contig_t *ctg) {
	clear_ctg(ctg);
	push_contigsv(merger->cache, ctg);
}

void reset_merger(merge_t *merger) {
	uint32_t i;
	
	for (i = 0; i < merger->ctgs->size; i++) {
		put_cache_ctgs(merger, get_contigv(merger->ctgs, i));
	}

//	for (i = 0; i < count_contigsv(merger->cache); i++) {
//		return_ctgv_merger(merger, get_contigsv(merger->cache, i));
//	}
//	clear_contigsv(merger->cache);
	clear_contigv(merger->ctgs);
//	for(i=0;i<vec_size(merger->ef->ctgs);i++){ put_pool_ctg(merger->ef, gget_vec(merger->ef->ctgs, i, FContig*)); }
//	clear_vec(merger->ef->ctgs);
	/*
	for (i = 0; i < vec_size(merger->ef->pool_ctg); i++) {
		ctg = gget_vec(merger->ef->pool_ctg, i, FContig*);
		free_vec(ctg->rids);
		free_string(ctg->seq);
		free(ctg);
	}*/
//	clear_vec(merger->ef->pool_ctg);
	//clear_ctgkmerv(merger->kmers);
	
	//clear_ctgkmerv(merger->aux_kmers);
}

void free_merger(merge_t *merger) {
	uint32_t i;

	for (i = 0; i < merger->ctgs->size; i++) {
		put_cache_ctgs(merger, get_contigv(merger->ctgs, i));
	}

	free_ctgs(merger);
	free_contigsv(merger->cache);
	free_contigv(merger->ctgs);
//	free_ctgset(merger->ctgs);
	if (merger->ef)
		free_ef(merger->ef);
	
	free(merger);
}

void merge_along_tree(merge_t *merger, pathtree_t *tree) {
	contig_t *c1, *c2;
	if (tree->tid || (tree->left == NULL && tree->right == NULL))
		return ;
	if (tree->left->left == NULL && tree->left->right == NULL && tree->right->left == NULL && tree->right->right == NULL) {
		c1 = get_contigv(merger->ctgs, tree->left->tid-1);
		c2 = get_contigv(merger->ctgs, tree->right->tid-1);
		if (is_similar_enough(merger, c1, c2)) {
			merger->sim_pairs++;
			merge_leaves(merger, tree->left->tid-1, tree->right->tid-1);
			tree->tid = c1->id<c2->id?c1->id:c2->id;
//			fprintf(stderr, "alongtree %d %d\n", ref_contigv(merger->ctgs, tree->left->tid-1)->id, ref_contigv(merger->ctgs, tree->right->tid-1)->id);
		}

	}
	if (tree->left->left && tree->left->right)
		merge_along_tree(merger, tree->left);
	if (tree->right->left && tree->right->right)
		merge_along_tree(merger, tree->right);
}

void merge_ctgs(merge_t *merger, FileReader *in, FILE *out) {
	uint32_t lastcid, cid;
	int n_col;
	lastcid = 0;
	while((n_col = fread_table(in)) != -1){
		if(n_col == 0) continue;
		cid = atoi(get_col_str(in, 4));
		if (cid != lastcid) {
			if (lastcid) { 
				build_tree(merger);
				update_ctg2merge(merger);
				if (merger->ctgs->size <= merger->max_cluster) {
					do {
						merger->sim_pairs = 0;
						merge_along_tree(merger, merger->tree); 
					} while (merger->sim_pairs);
					if (merger->ctgs->size>=3){ 
	//				if (merger->ctgs->size>=4 && merger->ctgs->size<=200){ 
	//					index_ctgs(merger);
						merge_core(merger);
					}
					if (merger->need_asm)
						print_asm(merger, out);
					else
						print_asm2(merger, out);

				} else {
					merge_along_tree(merger, merger->tree); 
					if (merger->need_asm)
						print_asm(merger, out);
					else
						print_asm2(merger, out);
				}
				free_tree(merger);
//				free_ctgs(merger);
				reset_merger(merger);
			}
			lastcid = cid;
			froll_back(in);
			prepare_reads(merger, in, lastcid);
		} else {
			prepare_reads(merger, in, lastcid);
		}
	}
	if (lastcid) {
		build_tree(merger);
		update_ctg2merge(merger);
		if (merger->ctgs->size <= merger->max_cluster) {
			do {
				merger->sim_pairs = 0;
				merge_along_tree(merger, merger->tree);
			} while (merger->sim_pairs);
			if (merger->ctgs->size>=3){ 
	//		if (merger->ctgs->size>=4 && merger->ctgs->size<=200){ 
	//			index_ctgs(merger);
				merge_core(merger);
			}
			if (merger->need_asm)
				print_asm(merger, out);
			else
				print_asm2(merger, out);
		} else {
			merge_along_tree(merger, merger->tree); 
			if (merger->need_asm)
				print_asm(merger, out);
			else
				print_asm2(merger, out);
		}
		free_tree(merger);
//		free_ctgs(merger);
		reset_merger(merger);
		//free_merger(merger);
	}
}
