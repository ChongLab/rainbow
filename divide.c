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

uint32_t call_key_col(Div *div, uint32_t gid){
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
			c  = base_bit_table[(int)get_u8list(div->seqs, rd->seqoff + col)];
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

void dividing_core(Div *div, uint32_t gid, int dep){
	ReadInfo *rd;
	u32list *grp, *sub;
	uint64_t mark0;
	uint32_t i, col, rid, gids[2], b;
	col = call_key_col(div, gid);
	b = col & 0x03;
	col >>= 2;
	if(col >= div->n_col){
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
	mark0 = get_u64list(div->markers, gid);
	if(dep < 64){
		push_u64list(div->markers, mark0);
		push_u64list(div->markers, mark0 | (1LLU << dep));
	} else {
		push_u64list(div->markers, mark0);
		push_u64list(div->markers, mark0);
	}
	for(i=0;i<count_u32list(grp);i++){
		rid = get_u32list(grp, i);
		rd = ref_rilist(div->rds, rid);
		if(rd->seqlen1 >= col && base_bit_table[(int)get_u8list(div->seqs, rd->seqoff + col)] == b){
			push_u32list(get_u32slist(div->grps, gids[1]), rid);
		} else {
			push_u32list(get_u32slist(div->grps, gids[0]), rid);
		}
	}
	for(i=0;i<2;i++){
		if(count_u32list(get_u32slist(div->grps, gids[i])) > 2 * div->k_allele) dividing_core(div, gids[i], dep + 1);
	}
}

void dividing(Div *div, uint32_t old_gid, FILE *out){
	ReadInfo *rd;
	u32list *grp;
	uint64_t marker;
	uint32_t i, j, gid, dep;
	char route[65];
	clear_u64list(div->markers);
	push_u64list(div->markers, 0);
	dividing_core(div, 0, 0);
	for(i=0;i<count_u32list(div->gids);i++){
		grp = get_u32slist(div->grps, get_u32list(div->gids, i));
		dep = get_u32list(div->deps, i);
		marker = get_u64list(div->markers, get_u32list(div->gids, i));
		for(j=0;j<dep;j++){
			route[j] = '0' + ((marker >> j) & 0x01);
		}
		route[dep] = 0;
		gid = ++div->gidoff;
		for(j=0;j<count_u32list(grp);j++){
			rd = ref_rilist(div->rds, get_u32list(grp, j));
			fprintf(out, "%u\t%u\t%s\t%s\t%u\t%s\n",
				rd->seqid, gid,
				as_array_u8list(div->seqs) + rd->seqoff, 
				as_array_u8list(div->seqs) + rd->seqoff + rd->seqlen1 + 1, old_gid, route);
		}
	}
	fflush(out);
	old_gid = old_gid;
}

Div* init_div(uint32_t k_allele, uint32_t K_allele, float min_freq){
	Div *div;
	div = malloc(sizeof(Div));
	div->gidoff   = 0;
	div->n_col    = 0;
	div->k_allele = k_allele;
	div->K_allele = K_allele;
	div->min_freq = min_freq;
	div->rds   = init_rilist(128);
	div->seqs  = init_u8list(128 * 80);
	div->grps  = init_u32slist(64);
	div->markers = init_u64list(64);
	div->deps = init_u32list(64);
	div->cache = init_u32slist(64);
	div->gids  = init_u32list(8);
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
	clear_u64list(div->markers);
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
	free_u32slist(div->cache);
	free_u32list(div->gids);
	free_u64list(div->markers);
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
		rd->seqoff  = count_u8list(div->seqs);
		rd->seqlen1 = get_col_len(fr, 2);
		rd->seqlen2 = get_col_len(fr, 3);
		for(i=0;i<get_col_len(fr, 2);i++) push_u8list(div->seqs, seq1[i]);
		push_u8list(div->seqs, 0);
		for(i=0;i<get_col_len(fr, 3);i++) push_u8list(div->seqs, seq2[i]);
		push_u8list(div->seqs, 0);
		push_u32list(get_u32slist(div->grps, 0), rid);
	}
	if(last_gid) dividing(div, last_gid, out);
	return ret;
}

