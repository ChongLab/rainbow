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

void prepare_reads(merge_t *merger, FileReader *in, uint32_t lastcid) {
	int n_col;
	uint32_t cid, ef_id, eid, eflen;
	ef_id = 0; 
	contig_t *ctg = NULL;
	read_t *rd;
	char *efstr;

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
			ctg = next_ref_contigv(merger->ctgs);
			ctg->id = ef_id;
			ctg->closed = 0;
			ctg->path = strdup(get_col_str(in, 5));
			ctg->efctgs = init_vec(sizeof(FContig*), 6);
			ctg->seq = NULL;
			ctg->sec_seq = NULL;
			ctg->rds = init_readv(24);
			rd = ref_next_readv(ctg->rds);
			rd->seq_id = atol(get_col_str(in, 0));
			rd->rd_len = eflen;
			memcpy(rd->seq, efstr, rd->rd_len);
			rd->seq[rd->rd_len]= '\0';
			rd->rank = 1;

			/*
			ctgh.id = ef_id;
			ctg = prepare_ctgset(merger->ctgs, ctgh, &exists);
			if (exists) {
			} else {
				ctgh.seq = NULL;
				ctgh.path = strdup(get_col_str(in, 5));
				ctgh.rds = init_readv(56);
				*ctg = ctgh;
				rd = ref_next_readv(ctg->rds);
				rd->seq_id = atol(get_col_str(in, 0));
				rd->rd_len = eflen;
				memcpy(rd->seq, efstr, rd->rd_len);
				rd->seq[rd->rd_len]= '\0';
				rd->rank = 1;
			}*/
		}
		rd = next_ref_readv(ctg->rds);
		rd->rank = 1;
		rd->seq_id = atol(get_col_str(in, 0));
		rd->rd_len = get_col_len(in, 3);
		memcpy(rd->seq, get_col_str(in, 3), rd->rd_len);
//		if (vec_size(ctg->ef->rds) <= merger->max_read) {
//			add_read2ef(ctg->ef, get_col_str(in, 3), seqid, get_col_len(in, 3), rank);
//		}
	}
}
void free_ctg(contig_t *ctg) {
	if (ctg->seq) free(ctg->seq);
	if (ctg->sec_seq) free(ctg->sec_seq);
	free(ctg->path);
	ctg->seq = ctg->sec_seq = ctg->path = NULL;
	if (ctg->rds->size) free_readv(ctg->rds);
}
void free_ctgs(merge_t *merger) {
	contig_t *ctg;
	uint32_t i;

	for (i = 0; i < merger->ctgs->size; i++) {
		ctg = ref_contigv(merger->ctgs, i);
//		put_pool_ctg(merger->ef, gget_vec(ctg->efctgs, i, FContig*));
		free_vec(ctg->efctgs);
		free_ctg(ctg);
	}
}

void print_asm(merge_t *merger, FILE *out) {
	uint32_t i, k;
//	uint32_t i, j, k, cid;
	contig_t *ctg; //FContig *c;
	read_t *rd;

//	reset_iter_ctgset(merger->ctgs);
	for (k = 0; k < merger->ctgs->size; k++) {
		ctg = ref_contigv(merger->ctgs, k);
//	while ((ctg = ref_iter_ctgset(merger->ctgs))) {
		if (merger->flag && !ctg->closed && ctg->rds->size >= merger->min_read && ctg->rds->size <= merger->max_read) {  // have used ef
			
			rd = ref_readv(ctg->rds, 0);
			reset_ef(merger->ef, ctg->id, rd->seq, rd->rd_len, merger->min_ol, merger->min_sm);
			for (i = 1; i < ctg->rds->size; i++) {
				rd = ref_readv(ctg->rds, i);
				add_read2ef(merger->ef, rd->seq, rd->seq_id, rd->rd_len, rd->rank);
			}
			align_reads_ef(merger->ef);
			asm_ef_ctgs(merger->ef);
			output_ef_ctgs(merger->ef, out);
			
			/*
			cid = 0;
			fprintf(out, "E %u\n", ctg->id);
			for(i=0;i<vec_size(ctg->efctgs);i++){
				c = gget_vec(ctg->efctgs, i, FContig*);
				if(c->closed) continue;
				fprintf(out, "C %u\n", cid);
				cid ++;
				fprintf(out, "L %d\n", (int)strlen(c->seq->string));
				fprintf(out, "S %s\n", c->seq->string);
				fprintf(out, "N %u\n", (uint32_t)ctg->rds->size);
				fprintf(out, "R");
				for(j=0;j<ctg->rds->size;j++){
					rd = ref_readv(ctg->rds, j);
					fprintf(out, " %u:%u", rd->seq_id, rd->rank);
				}
				fprintf(out, "\n//\n");
				fflush(out);
			} */
			/*
			if (!ctg->closed && ctg->rds->size >= merger->min_read && ctg->rds->size <= merger->max_read) {
				fprintf(out, ">%d\n%s\n", ctg->id, ctg->seq);
				fflush(out); 
			}
			*/
		}
	}

}

void update_ctg2merge(merge_t *merger) {
	uint32_t i, j;
	contig_t *ctg;
	read_t *rd;

//	reset_iter_ctgset(merger->ctgs);
//	while ((ctg = ref_iter_ctgset(merger->ctgs))) {
	for (j = 0; j < merger->ctgs->size; j++) {
		ctg = ref_contigv(merger->ctgs, j);
		if (merger->flag) { // have used ef
			rd = ref_readv(ctg->rds, 0);
			reset_ef(merger->ef, ctg->id, rd->seq, rd->rd_len, merger->min_ol, merger->min_sm);
			for (i = 1; i < ctg->rds->size; i++) {
				rd = ref_readv(ctg->rds, i);
				add_read2ef(merger->ef, rd->seq, rd->seq_id, rd->rd_len, rd->rank);
			}
			align_reads_ef(merger->ef);
			asm_ef_ctgs(merger->ef);
			assign_best_ctg(merger, ctg);
		} else { // first use ef
			merger->flag = 1;
			rd = ref_readv(ctg->rds, 0);
			merger->ef = init_ef(ctg->id, rd->seq, rd->rd_len, merger->min_ol, merger->min_sm);
			for (i = 1; i < ctg->rds->size; i++) {
				rd = ref_readv(ctg->rds, i);
				add_read2ef(merger->ef, rd->seq, rd->seq_id, rd->rd_len, rd->rank);
			}
			align_reads_ef(merger->ef);
			asm_ef_ctgs(merger->ef);
			assign_best_ctg(merger, ctg);
		}
	}
}

void assign_best_ctg(merge_t *merger, contig_t *ctg) {
	uint32_t i; int len = 0, sec_len = 0;
	FContig *c = NULL;
	

//	for(i=0;i<vec_size(ctg->efctgs);i++){ put_pool_ctg(merger->ef, gget_vec(ctg->efctgs, i, FContig*)); }
	clear_vec(ctg->efctgs);

	for(i=0;i<vec_size(merger->ef->ctgs);i++) {gpush_vec(ctg->efctgs, gget_vec(merger->ef->ctgs, i, FContig*), FContig*);}

	for (i = 0; i < vec_size(merger->ef->ctgs); i++) {
		c = gget_vec(merger->ef->ctgs, i, FContig*);
		if (c->closed) continue;
		if (len < c->seq->size) {
			sec_len = len;
			ctg->sec_seq = ctg->seq;
			len = c->seq->size;
			ctg->seq = c->seq->string;
		}
		if (sec_len < c->seq->size && len != c->seq->size) {
			sec_len = c->seq->size;
			ctg->sec_seq = c->seq->string;
		}
	}
	ctg->seq = strdup(ctg->seq);
	if (ctg->sec_seq) ctg->sec_seq = strdup(ctg->sec_seq);

}

static void print_leaf(merge_t *merger, pathtree_t *tree, FILE *out) {
	contig_t *ctg;
	if (tree->tid) {
		ctg = ref_contigv(merger->ctgs, tree->tid-1);
		fprintf(out, ">%d\n", ctg->id);
		fprintf(out, "%s\n", ctg->seq);
		fflush(out);
	} else {
		print_leaf(merger, tree->left, out);
		print_leaf(merger, tree->right, out); 
	}
	/*
	if (tree->left == NULL && tree->right == NULL) {
		c.id = tree->tid;
		ctg = get_ctgset(merger->ctgs, c);
		fprintf(out, ">%d\n", ctg->id);
		fprintf(out, "%s\n", ctg->seq);
		fflush(out);
	}
	print_leaf(merger, tree->left, out);
	print_leaf(merger, tree->right, out);*/
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

void merge_leaves(merge_t *merger, uint32_t id1, uint32_t id2) {
	read_t *rd;
	char *prefix;
	uint32_t i; int n;
	contig_t *c1, *c2, *c;
	c1 = ref_contigv(merger->ctgs, id1);
	c2 = ref_contigv(merger->ctgs, id2);
	c = id1<id2?c1:c2;
	if (c == c1) {
		c2->closed = 1;
	} else {
		c1->closed = 1;
	}
	n = strlen(c1->path)>=strlen(c2->path)?strlen(c2->path):strlen(c1->path);
	prefix = (char*)malloc(sizeof(char)*(n+1));
	memset(prefix, 0, n+1);
	prefix_path(c1->path, c2->path, n, prefix);
	rd = ref_readv(c1->rds, 0);
//	reset_ef(merger->ef, c1->id, rd->seq, rd->rd_len, merger->min_ol, merger->min_sm);
	for (i = 1; i < c1->rds->size; i++) {
		rd = ref_readv(c1->rds, i);
//		add_read2ef(merger->ef, rd->seq, rd->seq_id, rd->rd_len, rd->rank);
		if (c == c2) push_readv(c->rds, *rd);
	}
	for (i = 0; i < c2->rds->size; i++) {
		rd = ref_readv(c2->rds, i);
//		add_read2ef(merger->ef, rd->seq, rd->seq_id, rd->rd_len, rd->rank);
		if (c == c1) push_readv(c->rds, *rd);
	}
	free(c->path);
	c->path = NULL;
	c->path = strdup(prefix);
//	if (c->seq) free(c->seq);
//	if (c->sec_seq) free(c->sec_seq);
//	c->seq = c->sec_seq = NULL;
	if (strlen(c1->seq) > strlen(c2->seq)) {
		if (c==c2) { 
			if (c->seq) free(c->seq);
			c->seq = strdup(c1->seq);
		}
	} else {
		if (c==c1) {
			if (c->seq) free(c->seq);
			c->seq = strdup(c2->seq);
		}
	}

	if (c1->sec_seq && !c2->sec_seq) {
		if (c == c2) {
			c->sec_seq = strdup(c1->sec_seq);
		}
	}

	if (!c1->sec_seq && c2->sec_seq) {
		if (c == c1) {
			c->sec_seq = strdup(c2->sec_seq);
		}
	}

	if (c1->sec_seq && c2->sec_seq) {
		if (strlen(c1->sec_seq) > strlen(c2->sec_seq)) {
			if (c == c2) {
				if (c->sec_seq) free(c->sec_seq);
				c->sec_seq = strdup(c1->sec_seq);
			}
		} else {
			if (c == c1) {
				if (c->sec_seq) free(c->sec_seq);
				c->sec_seq = strdup(c2->sec_seq);
			}
		} 
	}
	//	align_reads_ef(merger->ef);
//	asm_ef_ctgs(merger->ef);
//	output_ef_ctgs(merger->ef, stderr);
//	assign_best_ctg(merger, c);
	free(prefix);
} 

void merge_along_tree(merge_t *merger, pathtree_t *tree) {
	contig_t *c1, *c2;
	if (tree->tid || (tree->left == NULL && tree->right == NULL))
		return ;
	if (tree->left->left == NULL && tree->left->right == NULL && tree->right->left == NULL && tree->right->right == NULL) {
		c1 = ref_contigv(merger->ctgs, tree->left->tid-1);
		c2 = ref_contigv(merger->ctgs, tree->right->tid-1);
		if (is_similar_enough(merger, c1, c2)) {
			merger->sim_pairs++;
			merge_leaves(merger, tree->left->tid-1, tree->right->tid-1);
			tree->tid = c1->id<c2->id?c1->id:c2->id;
//			tree->tid = tree->left->tid;
//			free(tree->left);
//			free(tree->right);
//			tree->left = tree->right = NULL;
//			fprintf(stderr, "alongtree %d %d\n", ref_contigv(merger->ctgs, tree->left->tid-1)->id, ref_contigv(merger->ctgs, tree->right->tid-1)->id);
//			return;
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
			if (lastcid) { // TODO merger here
				build_tree(merger);
				update_ctg2merge(merger);
				do {
					merger->sim_pairs = 0;
					merge_along_tree(merger, merger->tree); 
				} while (merger->sim_pairs);

//				printf("count=%d\n", count_ctgset(merger->ctgs));
				if (merger->ctgs->size>=4){ 
					index_ctgs(merger);
					merge_core(merger);
				}// print_leaf(merger, merger->tree, out);
				print_asm(merger, out);
				free_tree(merger);
				free_ctgs(merger);
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
		do {
			merger->sim_pairs = 0;
			merge_along_tree(merger, merger->tree);
		} while (merger->sim_pairs);
//		printf("count=%d\n", count_ctgset(merger->ctgs));
//		if (count_ctgset(merger->ctgs)>1) print_leaf(merger, merger->tree, out);
		if (merger->ctgs->size>=4){ 
			index_ctgs(merger);
			merge_core(merger);
		}// print_leaf(merger, merger->tree, out);
		print_asm(merger, out);
		free_tree(merger);
		free_ctgs(merger);
		reset_merger(merger);
		//free_merger(merger);
	}
}
void prepare_ctgs(merge_t *merger, uint32_t i, contig_t *ctg, uint64_t pos) {
	clear_ctgkmerv(merger->kmers);
	clear_ctgkmerv(merger->aux_kmers);
	clear_u32list(merger->ids);
	uint32_t j, n;
	int seqlen, idx, count, pre;
	ctg_kmer_t K, *t, *tpos;
	uint64_t kmask = 0xFFFFFFFFFFFFFFFFLLU >> ((32-merger->CTG_KMER_SIZE)*2), next, bt, p;

	K.kmer = 0;
	K.kpos = 0;
	K.id = 0;
	K.offset = -1;
	K.offset2 = -1;
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

//int is_onegroup(merge_t *merger, contig_t *c1, contig_t *c2) {
//	uint32_t i, j;
//	contigv *ctgs = merger->ctgs;
//	for (i = c1->clsid; i != (ref_contigv(ctgs, i))->clsid; i = (ref_contigv(ctgs, i))->clsid) ;
//	for (j = c2->clsid; j != (ref_contigv(ctgs, j))->clsid; j = (ref_contigv(ctgs, j))->clsid) ;
//	
//	if (i == j) return 1;
//	return 0;
//}

int is_similar_enough(merge_t *merger, contig_t *c1, contig_t *c2) {
	int mm, mn, aln_len, i, ii;
	AlnAln *aa = NULL;
	AlnParam ap = {10, 2, 2, aln_sm_nt, 16, 75};

	mn = mm = 0;
	for (ii = 0; ii < 4; ii++) {
		switch (ii) {
			case 0: aa = aln_stdaln(c1->seq, c2->seq, &ap, 0, 1); break;
			case 1: if (c2->sec_seq) aa = aln_stdaln(c1->seq, c2->sec_seq, &ap, 0, 1); break;
			case 2: if (c1->sec_seq) aa = aln_stdaln(c1->sec_seq, c2->seq, &ap, 0, 1); break;
			case 3: if (c1->sec_seq && c2->sec_seq) aa = aln_stdaln(c1->sec_seq, c2->sec_seq, &ap, 0, 1); break;

		}
		if (!aa) continue;
		aln_len = strlen(aa->out1);
		for (i = 0; i <aln_len; i++) {
			if (aa->out1[i] == '-' || aa->out2[i] == '-')
				continue;
			if (aa->out1[i] != aa->out2[i])
				mm++;
			mn++;
		}
		
//		fprintf(stderr, "%d %d\n%s\n%s\n", c1->id, c2->id, aa->out1, aa->out2);
//		fflush(stderr);
		aln_free_AlnAln(aa); aa = NULL;
		if (mn > (int)merger->min_overlap && ((float)mm/mn - merger->het) <= 0) {
			return 1; 
		}
	}
	return 0;
}

//void merge_2ctg(merge_t *merger, contig_t *c1, contig_t *c2) {
//	uint32_t i, j;
//	contigv *ctgs = merger->ctgs;
//	for (i = c1->clsid; i != (ref_contigv(ctgs, i))->clsid; i = (ref_contigv(ctgs, i))->clsid)
//		ref_contigv(ctgs, i)->clsid = ref_contigv(ctgs, ref_contigv(ctgs, i)->clsid)->clsid;
//	for (j = c2->clsid; j != (ref_contigv(ctgs, j))->clsid; j = (ref_contigv(ctgs, j))->clsid) ;
//		ref_contigv(ctgs, j)->clsid = ref_contigv(ctgs, ref_contigv(ctgs, j)->clsid)->clsid;
//	
//	if (i == j) return;
//	
//	if (ref_contigv(ctgs, i)->sz < ref_contigv(ctgs, j)->sz) {
//		ref_contigv(ctgs, i)->clsid = j;
//		ref_contigv(ctgs, j)->sz += ref_contigv(ctgs, i)->sz;
//	} else {
//		ref_contigv(ctgs, j)->clsid = i;
//		ref_contigv(ctgs, i)->sz += ref_contigv(ctgs, j)->sz;
//	}
//
//	return;
//}

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

//void gather_leaves(pathtree_t *tree, u32list *leaves) {
//	if (tree == NULL)
//		return;
//	
//	if (tree->left == NULL && tree->right == NULL) {
//		push_u32list(leaves, tree->tid);
//	} else {
//		gather_leaves(tree->left, leaves);
//		gather_leaves(tree->right, leaves);
//	}
//
//}
//
//void update_merger(merge_t *merger, contig_t *c1, contig_t *c2) {
//	u32list *subtree;
//	pathtree_t *t; int i, n;
//	char *prefix;
//	n = strlen(c1->path)>=strlen(c2->path)?strlen(c2->path):strlen(c1->path);
//	prefix = (char*)malloc(sizeof(char)*(n+1));
//	memset(prefix, 0, n+1);
//	prefix_path(c1->path, c2->path, n, prefix);
//	if (strcmp(prefix, "") == 0) {
//		free(prefix);
//		return;
//	}
//	subtree = init_u32list(2);
//	n = strlen(prefix);
//	t = merger->tree;
//	for (i = 0; i < n; i++) {
//		if (prefix[i] == '0') {
//			t = t->left;
//		} else if (prefix[i] == '1') {
//			t = t->right;
//		}
//		if (t == NULL) {
//			free(prefix);
//			free_u32list(subtree);
//			return;
//		}
//	}
//	gather_leaves(t, subtree);
//	n = count_u32list(subtree);
//	for (i = 0; i < n; i++) {
//		merge_2ctg(merger, c1, ref_contigv(merger->ctgs, get_u32list(subtree, i)));
//	}
//	free(prefix);
//	free_u32list(subtree);
//	return;
//}

void merge_core(merge_t *merger) {
	contigv *ctgs = merger->ctgs;
	contig_t *ctg1, *ctg2; uint32_t i, j, ii;
	int seqlen;
	uint64_t pos = 0;
	for (i = 0; i < ctgs->size-1; i++) {
		ctg1 = ref_contigv(ctgs, i);
		if (ctg1->closed) continue;
		seqlen = strlen(ctg1->seq);
		if (seqlen < (int)merger->CTG_KMER_SIZE) continue;
		prepare_ctgs(merger, i, ctg1, pos);
		for (ii = 0; ii < merger->ids->size; ii++) {
			j = get_u32list(merger->ids, ii);
			ctg2 = ref_contigv(ctgs, j);
			if (ctg2->closed) continue;
//			if (is_onegroup(merger, ctg1, ctg2)) {
//				continue;
//			} else {
			if (is_similar_enough(merger, ctg1, ctg2)) {
				merge_leaves(merger, i, j);
//				fprintf(stderr, "finalleaves %d %d\n", ctg1->id, ctg2->id);
//				merge_2ctg(merger, ctg1, ctg2);
//				update_merger(merger, ctg1, ctg2);
			}
//			}
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
//	while ((ctg = ref_iter_ctgset(merger->ctgs)) != NULL) {
		ctg = ref_contigv(merger->ctgs, i);	
//		i = ctg->id;
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
//void print_clusters(merge_t *merger, FILE *out) {
//	int last_cid = -1;
//	int line_num = 0;
//	uint32_t i, n;
//	contig_t *ctg;
//	contigv *ctgs = merger->ctgs;
//	
//	qsort(as_array_contigv(ctgs), count_contigv(ctgs), sizeof(contig_t), cmp_ctg_clsids);
//	n = count_contigv(ctgs);
//	for (i = 0; i < n; i++) {
//		ctg = ref_contigv(ctgs, i);
//		if (last_cid != (int)ctg->clsid && line_num > 0)
//			fprintf(out, "\n");
//		last_cid = ctg->clsid;
//		line_num++;
//		fprintf(out, "%d ", ctg->id);
//	}
//	fprintf(out, "\n");
//	fflush(out);
//	return;
//
//}

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
//	merger->index = init_ctgkhash(23);
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
//	merger->ctgs = init_ctgset(2);
	merger->tree = NULL;
	merger->index = init_ctgkhash(23);
	merger->links = NULL;
	merger->idv = NULL;
	merger->kmers = init_ctgkmerv(23);
	merger->aux_kmers = init_ctgkmerv(12);
	merger->ids = init_u32list(2);
	merger->min_kmer = min_kmer;
	merger->min_overlap = min_overlap;
	merger->het = het;
	merger->CTG_KMER_SIZE = kmersize;
	merger->min_ol = 5;              // TODO: parameter added
	merger->min_sm = 0.90;
	merger->min_read = 5;
	merger->max_read = 300;
	merger->sim_pairs = 0;
	merger->ef = NULL;
	merger->flag = 0;
	return merger;
}

void reset_merger(merge_t *merger) {
	FContig *ctg; uint32_t i;
	clear_contigv(merger->ctgs);
	clear_ctgkhash(merger->index);
	for(i=0;i<vec_size(merger->ef->ctgs);i++){ put_pool_ctg(merger->ef, gget_vec(merger->ef->ctgs, i, FContig*)); }
	clear_vec(merger->ef->ctgs);
	for (i = 0; i < vec_size(merger->ef->pool_ctg); i++) {
		ctg = gget_vec(merger->ef->pool_ctg, i, FContig*);
		free_vec(ctg->rids);
		free_string(ctg->seq);
		free(ctg);
	}
	clear_vec(merger->ef->pool_ctg);
	if (merger->links) free(merger->links);
	if (merger->idv) free(merger->idv);
	merger->links = NULL;
	merger->idv = NULL;
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
/*
void free_ctg(contig_t *ctg) {
//	free(ctg->seq);
	free(ctg->path);
//	free_ef(ctg->ef);
}
*/
void free_merger(merge_t *merger) {
	free_contigv(merger->ctgs);
	free_ctgkhash(merger->index);
//	free_ctgset(merger->ctgs);
	if (merger->ef)
		free_ef(merger->ef);
	free_ctgkmerv(merger->kmers);
	free_ctgkmerv(merger->aux_kmers);
	free_u32list(merger->ids);
	
	free(merger);
}
