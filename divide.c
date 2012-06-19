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
 
#include "rainbow.h"

u32list* lend_ulist_div(Div *div){
	u32list *list;
	if(!pop_u32slist(div->cache, &list)) list = init_u32list(4);
	return list;
}

void return_ulist_div(Div *div, u32list *list){ if(list){ clear_u32list(list); push_u32slist(div->cache, list); } }

typedef struct {
	uint32_t cnt;
	uint32_t base;
} BaseCnt;

static inline uint32_t C_N_2(uint32_t n){
	if(n == 0) return 0;
	else return n * (n - 1) / 2;
}
uint32_t _call_key_col(Div *div, uint32_t gid){
	ReadInfo *rd;
	u32list *grp;
	uint32_t i, j, col, row, key, c, max_non, tol, base;
	BaseCnt cnts[4];
	key = div->n_col;
	base = 0;
	grp = get_u32slist(div->grps, gid);
	max_non = 0;
	for(col=0;col<div->n_col;col++){
		for(i=0;i<4;i++){ cnts[i].base = i; cnts[i].cnt = 0; }
		for(row=0;row<count_u32list(grp);row++){
			rd = ref_rilist(div->rds, get_u32list(grp, row));
			if(rd->seqlen1 <= col) continue;
//			c  = base_bit_table[(int)get_u8list(div->seqs, rd->seqoff + col)];
			c = div->seqs->buffer[rd->seqoff + col];
			cnts[c&0x03].cnt ++;
		}
		tol = cnts[0].cnt + cnts[1].cnt + cnts[2].cnt + cnts[3].cnt;
		if(tol == 0) break;
		for(i=0;i<2;i++){
			for(j=3;j>i;j--){
				if(cnts[j].cnt > cnts[j-1].cnt){
					swap_tmp(cnts[j].cnt, cnts[j-1].cnt, c);
					swap_tmp(cnts[j].base, cnts[j-1].base, c);
				}
			}
		}
		if(cnts[1].cnt < div->k_allele) continue;
		if(cnts[1].cnt < div->K_allele && cnts[1].cnt < div->min_freq * tol) continue;
		if(cnts[1].cnt > max_non){
			max_non = cnts[1].cnt;
			key = col;
			base = cnts[1].base;
		}
	}
	return (key << 2) | base;
}

uint32_t call_key_col(Div *div, uint32_t gid){
	ReadInfo *rd;
	u32list *grp;
	uint32_t i, j, k, col, row, key, c, max_non, tol, base, s1, s2;
	BaseCnt cnts[4];
	col_base_t *cb;
	uint64_t MM1, MM2;
	uint32_t n_p1, n_p2, idx;
	double min_mm, mm1, mm2;
	key = div->n_col;
	base = 0;
	grp = get_u32slist(div->grps, gid);
	max_non = 0;
	clear_cbv(div->cbs);
	for(col=0;col<div->n_col;col++){
		for(i=0;i<4;i++){ cnts[i].base = i; cnts[i].cnt = 0; }
		for(row=0;row<count_u32list(grp);row++){
			rd = ref_rilist(div->rds, get_u32list(grp, row));
			if(rd->seqlen1 <= col) continue;
			c = div->seqs->buffer[rd->seqoff + col];
			cnts[c&0x03].cnt ++;
		}
		tol = cnts[0].cnt + cnts[1].cnt + cnts[2].cnt + cnts[3].cnt;
		if(tol == 0) break;
		for(i=0;i<2;i++){
			for(j=3;j>i;j--){
				if(cnts[j].cnt > cnts[j-1].cnt){
					swap_tmp(cnts[j].cnt, cnts[j-1].cnt, c);
					swap_tmp(cnts[j].base, cnts[j-1].base, c);
				}
			}
		}
		if(cnts[1].cnt < div->k_allele) continue;
		if(cnts[1].cnt < div->K_allele && cnts[1].cnt < div->min_freq * tol) continue;
		cb = next_ref_cbv(div->cbs);
		cb->col  = col;
		cb->base = cnts[1].base;
		cb->cnt  = cnts[1].cnt;
	}
	if(div->cbs->size == 1){
		key = ref_cbv(div->cbs, 0)->col;
		base = ref_cbv(div->cbs, 0)->base;
	}
	if(div->cbs->size > 1){
		encap_u32list(div->ps1, div->n_col * 4);
		encap_u32list(div->ps2, div->n_col * 4);
		min_mm = 10000000;
		for(i=0;i<div->cbs->size;i++){
			cb = ref_cbv(div->cbs, i);
			n_p1 = cb->cnt;
			n_p2 = grp->size - cb->cnt;
			memset(div->ps1->buffer, 0, div->n_col * 4 * 4);
			memset(div->ps2->buffer, 0, div->n_col * 4 * 4); 
			for(row=0;row<grp->size;row++){
				rd = ref_rilist(div->rds, get_u32list(grp, row));
				if(rd->seqlen1 <= cb->col) idx = 1;
				else idx = (div->seqs->buffer[rd->seqoff + cb->col] != cb->base);
				if(idx){
					for(j=0;j<div->n_col;j++){
						div->ps2->buffer[div->seqs->buffer[rd->seqoff + j] + 4 * j] ++;
					}
				} else {
					for(j=0;j<div->n_col;j++){
						div->ps1->buffer[div->seqs->buffer[rd->seqoff + j] + 4 * j] ++;
					}
				}
			}
			MM1 = MM2 = 0; 
			for(j=0;j<div->n_col;j++){
				if(j == i) continue;
				s1 = C_N_2(n_p1);
				s2 = C_N_2(n_p2);
				for (k = 0; k < 4; k++) {
					s1 -= C_N_2(div->ps1->buffer[k + 4 * j]);
					s2 -= C_N_2(div->ps2->buffer[k + 4 * j]);
				}
				MM1 += s1;
				MM2 += s2;
				//fprintf(stdout, " -- %u %u in %s -- %s:%d --\n", s1, s2, __FUNCTION__, __FILE__, __LINE__);
			}
			//fprintf(stdout, "col%d mm1 %lld\n", cb->col, MM1);
			//fprintf(stdout, "col%d mm2 %lld\n", cb->col, MM2);
			mm1 = ((long double)MM1) / (n_p1*(n_p1-1)/2);
			mm2 = ((long double)MM2) / (n_p2*(n_p2-1)/2);
			if(mm1 < mm2) mm1 = mm2;
			//fprintf(stdout, "gid%u col%d %f\n", gid, cb->col, mm1);
			if(mm1 - min_mm < 0.00000000001){
				min_mm = mm1;
				key = cb->col;
				base = cb->base;
			}
		}
	}
	return (key << 2) | base;
}

void dividing_core(Div *div, uint32_t gid, int dep){
	ReadInfo *rd;
	u32list *grp, *sub;
//	uint64_t mark0;
	uint32_t i, j, col, rid, gids[2], b;
	col = call_key_col(div, gid);
	b = col & 0x03;
	col >>= 2;
	if(col >= div->n_col || div->rds->size < div->K_allele || dep > 255){
		push_u32list(div->gids, gid);
		push_u32list(div->deps, dep);
		return;
	}
	for(i=0;i<2;i++){
		gids[i] = count_u32slist(div->grps);
		sub = lend_ulist_div(div);
		push_u32slist(div->grps, sub);
	}
	grp = get_u32slist(div->grps, gid);
	/*
	char str[257];
	for(i=0;(int)i<dep;i++){
		mark0 = get_u64list(div->markers[i/64], gid);
		str[i] = '0' + ((mark0 >> (i%64))& 0x01);
	}
	str[i] = '\0';
	fprintf(stderr, "%s\t%d\t%c\n", str, col, "ACGT"[b]);
	for (j = 0; j < 4; j++) {
		push_u64list(div->markers[j], get_u64list(div->markers[j], gid));
		push_u64list(div->markers[j], get_u64list(div->markers[j], gid));
	}
	if (dep <= 255) {
		set_u64list(div->markers[dep/64], gid+1, get_u64list(div->markers[j], gid) | (1LLU << (dep%64)));
	}
	*/
	
//	if (dep <= 255) {
	for (j = 0; (int)j < dep/64; j++) {
		push_u64list(div->markers[j], get_u64list(div->markers[j], gid));
		push_u64list(div->markers[j], get_u64list(div->markers[j], gid));
	}
		push_u64list(div->markers[j], get_u64list(div->markers[j], gid));
		push_u64list(div->markers[j], get_u64list(div->markers[j], gid) | (1LLU << (dep%64)));
	j++;
	for (; j < 4; j++) {
		push_u64list(div->markers[j], 0);
		push_u64list(div->markers[j], 0);
	}
//	}

	for(i=0;i<count_u32list(grp);i++){
		rid = get_u32list(grp, i);
		rd = ref_rilist(div->rds, rid);
		if(rd->seqlen1 >= col && div->seqs->buffer[rd->seqoff + col] == b){
			push_u32list(get_u32slist(div->grps, gids[1]), rid);
		} else {
			push_u32list(get_u32slist(div->grps, gids[0]), rid);
		}
	}
	for(i=0;i<2;i++){
//		if(count_u32list(get_u32slist(div->grps, gids[i])) > 2 * div->k_allele) dividing_core(div, gids[i], dep + 1);
		dividing_core(div, gids[i], dep + 1);
	}
}

void dividing(Div *div, uint32_t old_gid, FILE *out){
	ReadInfo *rd;
	u32list *grp;
	uint64_t marker;
	uint32_t i, j, k, gid, dep;
	char route[257];
	String *seq1, *seq2;
	seq1 = init_string(1024);
	seq2 = init_string(1024);
	for (i = 0; i < 4; i++) {
		clear_u64list(div->markers[i]);
		push_u64list(div->markers[i], 0); 
	}
	dividing_core(div, 0, 0);
	for(i=0;i<count_u32list(div->gids);i++){
		grp = get_u32slist(div->grps, get_u32list(div->gids, i));
		dep = get_u32list(div->deps, i);
		if (dep>255) dep = 255;
//		marker1 = get_u64list(div->markers1, get_u32list(div->gids, i));
//		marker2 = get_u64list(div->markers2, get_u32list(div->gids, i));
		for(j=0;j<dep;j++){
			marker = get_u64list(div->markers[j/64], get_u32list(div->gids, i)); 
			route[j] = '0' + ((marker >> (j%64)) & 0x01);
		}
		route[dep] = 0;
		gid = ++div->gidoff;
		for(j=0;j<count_u32list(grp);j++){
			rd = ref_rilist(div->rds, get_u32list(grp, j));
			for(k=0;k<rd->seqlen1;k++) seq1->string[k] = bit_base_table[div->seqs->buffer[rd->seqoff + k]];
			seq1->string[k] = 0;
			for(k=0;k<rd->seqlen2;k++) seq2->string[k] = bit_base_table[div->seqs->buffer[rd->seqoff + rd->seqlen1 + k]];
			seq2->string[k] = 0;
			fprintf(out, "%u\t%u\t%s\t%s\t%u\t%s\n",
				rd->seqid, gid, seq1->string, seq2->string, old_gid, route);
		}
	}
	fflush(out);
	old_gid = old_gid;
	free_string(seq1);
	free_string(seq2);
}

Div* init_div(uint32_t k_allele, uint32_t K_allele, float min_freq){
	Div *div; int i;
	div = malloc(sizeof(Div));
	div->gidoff   = 0;
	div->n_col    = 0;
	div->k_allele = k_allele;
	div->K_allele = K_allele;
	div->min_freq = min_freq;
	div->rds   = init_rilist(128);
	div->seqs  = init_u8list(128 * 80);
	div->grps  = init_u32slist(64);
	for (i = 0; i < 4; i++) {
		div->markers[i] = init_u64list(64);
	}
	div->deps = init_u32list(64);
	div->cache = init_u32slist(64);
	div->gids  = init_u32list(8);
	div->cbs = init_cbv(12);
	div->ps1 = init_u32list(32);
	div->ps2 = init_u32list(32); 
	return div;
}

void reset_div(Div *div){
	uint32_t i;
	clear_rilist(div->rds);
	clear_u8list(div->seqs);
	for(i=0;i<count_u32slist(div->grps);i++){
		return_ulist_div(div, get_u32slist(div->grps, i));
	}
	clear_u32slist(div->grps);
	clear_u32list(div->gids);
	for (i = 0; i < 4; i++) {
		clear_u64list(div->markers[i]);
	}
	clear_u32list(div->deps);
	div->n_col = 0;
}

void free_div(Div *div){
	uint32_t i;
	reset_div(div);
	free_rilist(div->rds);
	free_u8list(div->seqs);
	free_u32slist(div->grps);
	for(i=0;i<count_u32slist(div->cache);i++){
		free_u32list(get_u32slist(div->cache, i));
	}
	free_u32list(div->ps1);
	free_u32list(div->ps2);
	free_cbv(div->cbs);
	free_u32slist(div->cache);
	free_u32list(div->gids);
	for (i = 0; i < 4; i++) {
		free_u64list(div->markers[i]);
	}
	free_u32list(div->deps);
	free(div);
}

uint32_t div_reads(Div *div, FileReader *fr, FILE *out){
	ReadInfo *rd;
	uint32_t seqid, rank, gid, last_gid, rid, ret;
	char *seq1, *seq2;
	int i;
	last_gid = 0;
	ret = 0;
	while(fread_table(fr) != -1){
		seqid = atoll(get_col_str(fr, 0));
		rank  = 1;
		gid   = atoll(get_col_str(fr, 1));
		seq1  = get_col_str(fr, 2);
		seq2  = get_col_str(fr, 3);
		if(gid != last_gid){
			ret ++;
			if(last_gid) dividing(div, last_gid, out);
			last_gid = gid;
			reset_div(div);
			push_u32slist(div->grps, lend_ulist_div(div));
		}
		if(get_col_len(fr, 2) > (int)div->n_col) div->n_col = get_col_len(fr, 2);
		rid = count_rilist(div->rds);
		rd  = next_ref_rilist(div->rds);
		rd->seqid   = seqid;
		rd->rank    = rank;
		rd->seqoff  = div->seqs->size;
		rd->seqlen1 = get_col_len(fr, 2);
		rd->seqlen2 = get_col_len(fr, 3);
		for(i=0;i<rd->seqlen1;i++) push_u8list(div->seqs, base_bit_table[(int)seq1[i]]);
		for(i=0;i<rd->seqlen2;i++) push_u8list(div->seqs, base_bit_table[(int)seq2[i]]);
		push_u32list(get_u32slist(div->grps, 0), rid);
	}
	if(last_gid) dividing(div, last_gid, out);
	return ret;
}

