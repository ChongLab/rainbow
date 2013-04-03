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
 
#include "asm_R2.h"

Vector* get_pool_vec(EF *ef){
	Vector *vec;
	vec = NULL;
	if(vec_size(ef->pool_vec)){
		gpop_vec(ef->pool_vec, vec, Vector*);
	} else {
		vec = init_vec(sizeof(rp_t), 64);
	}
	return vec;
}

void put_pool_vec(EF *ef, Vector *vec){
	clear_vec(vec);
	gpush_vec(ef->pool_vec, vec, Vector*);
}

FContig* get_pool_ctg(EF *ef){
	FContig* ctg;
	ctg = NULL;
	if(vec_size(ef->pool_ctg)){ // czc modified here to hold all ctgs
		gpop_vec(ef->pool_ctg, ctg, FContig*);
	} else {
		ctg = malloc(sizeof(FContig));
		ctg->rids   = init_vec(sizeof(uint32_t), 6);
		ctg->seq    = init_string(1024);
	}
	return ctg;
}

int cmp_ol_func(const void *e1, const void *e2){
	Overlap *o1, *o2;
	o1 = (Overlap*)e1;
	o2 = (Overlap*)e2;
	if(o1->l_ol < o2->l_ol) return 1;
	else if(o1->l_ol > o2->l_ol) return -1;
	else return 0;
}

void put_pool_ctg(EF *ef, FContig *ctg){
	clear_vec(ctg->rids);
	clear_string(ctg->seq);
	gpush_vec(ef->pool_ctg, ctg, FContig*);
}

void add_read2ef_core(EF *ef, char *seq, uint32_t seq_id, uint32_t rd_len, uint32_t rank){
	FRead *rd;
	FContig *ctg;
	rp_t *p;
	rhash_t RH, *rh;
	Vector *vs;
	uint32_t rid, cid, i, kmer;
	int exists;
	if(rd_len > MAX_RD_LEN) rd_len = MAX_RD_LEN;
	if(rd_len == 0) return;
	rid = vec_size(ef->rds);
	rd  = get_next_vec_ref(ef->rds);
	cid = vec_size(ef->ctgs);
	ctg = get_pool_ctg(ef);
	gpush_vec(ef->ctgs, ctg, FContig*);
	rd->seq_id  = seq_id;
	rd->rank    = rank;
	rd->rd_len  = rd_len;
	rd->ctg_id  = cid;
	rd->ctg_off = 0;
	rd->used    = 0;
	ctg->len    = rd_len;
	ctg->closed = 0;
	gpush_vec(ctg->rids, rid, uint32_t);
	append_string(ctg->seq, seq, rd_len);
	memcpy(rd->seq, seq, rd_len);
	rd->seq[rd_len] = 0;
	kmer = 0;
	RH.rps_idx = 0;
	for(i=0;i<rd_len;i++){
		kmer = (((kmer << 2) | base_bit_table[(int)rd->seq[i]])) & ASM_KMER_MASK;
		if(i + 1 < ASM_KMER_SIZE) continue;
		RH.kmer = kmer;
		rh = prepare_rhash(ef->index, RH, &exists);
		if(exists){
			vs = gget_vec(ef->rps, rh->rps_idx, Vector*);
		} else {
			rh->kmer  = kmer;
			rh->rps_idx = vec_size(ef->rps);
			vs = get_pool_vec(ef);
			gpush_vec(ef->rps, vs, Vector*);
		}
		p = get_next_vec_ref(vs);
		p->rid  = rid;
		p->roff = i;
	}
}

EF* init_ef(uint32_t ef_id, char *eseq, uint32_t rd_len, uint32_t min_ol, float min_sm){
	EF *ef;
	ef = malloc(sizeof(EF));
	ef->ef_id  = ef_id;
	ef->min_ol = min_ol;
	ef->min_sm = min_sm;
	ef->inc_tag = 1;
	ef->rds    = init_vec(sizeof(FRead), 64);
	ef->ols    = init_vec(sizeof(Overlap), 64);
	ef->rps    = init_vec(sizeof(Vector*), 64);
	ef->ctgs   = init_vec(sizeof(FContig*), 6);
	ef->index  = init_rhash(1023);
	ef->uniq   = init_u64hash(1023);
	ef->pool_vec = init_vec(sizeof(Vector*), 64);
	ef->pool_ctg = init_vec(sizeof(FContig*), 64);
	memcpy(ef->eseq, eseq, rd_len);
	ef->eseq[rd_len] = 0;
	add_read2ef_core(ef, eseq, ef_id, rd_len, 0);
	return ef;
}

void set_inc_tag_ef(EF *ef, uint32_t inc){
	ef->inc_tag = inc;
}

void add_read2ef(EF *ef, char *seq, uint32_t seq_id, uint32_t rd_len, uint32_t rank){ add_read2ef_core(ef, seq, seq_id, rd_len, (rank == 0)? 1 : rank); }


void find_overlap(char *seq1, uint32_t len1, uint32_t off1, char *seq2, uint32_t len2, uint32_t off2, uint32_t *l_ol, uint32_t *r_ol, uint32_t *n_mm){
	uint32_t i, l, r, ol, mm;
	l = (off1 <= off2)? off1 : off2;
	r = (len1 - off1 <= len2 - off2)? (len1 - off1) : (len2 - off2);
	*r_ol = *l_ol = ol = l + r;
	mm = 0;
	for(i=0;i<ol;i++){
		if(seq1[off1-l+i] != seq2[off2-l+i]) mm ++;
	}
	*n_mm = mm;
}

void align_reads_ef(EF *ef){
	FRead *rd1, *rd2;
	Vector *rp;
	Overlap *ol;
	rp_t   *p1, *p2;
	uint32_t i, j, k, r_ol, l_ol, n_mm;
	uint64_t *uniq, aln_id;
	int exists;
	for(i=0;i<vec_size(ef->rps);i++){
		rp = gget_vec(ef->rps, i, Vector*);
		for(j=0;j<vec_size(rp);j++){
			p1  = get_vec_ref(rp, j);
			rd1 = get_vec_ref(ef->rds, p1->rid);
			for(k=j+1;k<vec_size(rp);k++){
				p2  = get_vec_ref(rp, k);
				if(p2->rid > p1->rid){
					aln_id = (p1->rid << 16) | p2->rid;
				} else {
					aln_id = (p2->rid << 16) | p1->rid;
				}
				if(p1->roff >= p2->roff){
					aln_id |= (((uint64_t)p1->roff - p2->roff) << 32);
				} else {
					aln_id |= (((uint64_t)p2->roff - p1->roff) << 48);
				}
				uniq = prepare_u64hash(ef->uniq, aln_id, &exists);
				if(exists) continue;
				*uniq = aln_id;
				rd2 = get_vec_ref(ef->rds, p2->rid);
				find_overlap(rd1->seq, rd1->rd_len, p1->roff, rd2->seq, rd2->rd_len, p2->roff, &l_ol, &r_ol, &n_mm);
				if(r_ol < ef->min_ol && l_ol < ef->min_ol) continue;
				if(n_mm > (uint32_t)((1 - ef->min_sm) * r_ol) || n_mm > (uint32_t)((1 - ef->min_sm) * l_ol)) continue;
				ol = get_next_vec_ref(ef->ols);
				ol->used = 0;
				if(p1->roff >= p2->roff){
					ol->l_rid = p1->rid;
					ol->r_rid = p2->rid;
					ol->l_ol  = l_ol;
					ol->r_ol  = r_ol;
					ol->n_mm  = n_mm;
				} else {
					ol->l_rid = p2->rid;
					ol->r_rid = p1->rid;
					ol->l_ol  = r_ol;
					ol->r_ol  = l_ol;
					ol->n_mm  = n_mm;
				}
			}
		}
	}
	qsort_vec(ef->ols, cmp_ol_func);
}

void print_alignments(EF *ef){
	uint32_t i, j;
	Overlap *ol;
	FRead *r1, *r2;
	for(i=0;i<vec_size(ef->ols);i++){
		ol = get_vec_ref(ef->ols, i);
		r1 = get_vec_ref(ef->rds, ol->l_rid);
		r2 = get_vec_ref(ef->rds, ol->r_rid);
		printf("%u <-> %u = %u:%u\n", ol->l_rid, ol->r_rid, ol->l_ol, ol->n_mm);
		printf("%s\n", r1->seq);
		for(j=0;(int)j<r1->rd_len-ol->l_ol;j++) printf(" ");
		printf("%s\n", r2->seq);
	}
}

void asm_ef_ctgs(EF *ef){
	Overlap *ol;
	FRead *rd1, *rd2, *rd;
	FContig *ctg1, *ctg2;
	uint32_t i, off1, off2, l_ol, r_ol, n_mm, offset;
	uint32_t j, rid, rank_type;
	for(rank_type=0;rank_type<4;rank_type++){
		for(i=0;i<vec_size(ef->ols);i++){
			ol = get_vec_ref(ef->ols, i);
			if(ol->used) continue;
			if(ef->inc_tag == 0 && (ol->l_rid == 0 || ol->r_rid == 0)) continue;
			rd1 = get_vec_ref(ef->rds, ol->l_rid);
			rd2 = get_vec_ref(ef->rds, ol->r_rid);
			if(rank_type == 0){
				if(rd1->rank != rd2->rank) continue;
			} else if(rank_type == 1){
				if(rd1->rank + 1 != rd2->rank) continue;
			} else if(rank_type == 2){
				if(rd1->rank > rd2->rank) continue;
			}
			ctg1 = gget_vec(ef->ctgs, rd1->ctg_id, FContig*);
			ctg2 = gget_vec(ef->ctgs, rd2->ctg_id, FContig*);
			if(ctg1 == ctg2){ ol->used = 1; continue; }
			off1 = rd1->ctg_off + rd1->rd_len - ol->l_ol;
			off2 = rd2->ctg_off;
			find_overlap(ctg1->seq->string, ctg1->len, off1, ctg2->seq->string, ctg2->len, off2, &l_ol, &r_ol, &n_mm);
			if(l_ol < ef->min_ol && r_ol < ef->min_ol){ continue; }
			if(n_mm > (uint32_t)(l_ol * (1 - ef->min_sm)) || n_mm > (uint32_t)(r_ol * (1 - ef->min_sm))){ continue; }
			ol->used = 1;
			if(off1 >= off2){
				ctg2->closed = 1;
				offset = off1 - off2;
				for(j=0;j<vec_size(ctg2->rids);j++){
					rid = gget_vec(ctg2->rids, j, uint32_t);
					rd  = get_vec_ref(ef->rds, rid);
					gpush_vec(ctg1->rids, rid, uint32_t);
					rd->ctg_id  = rd1->ctg_id;
					rd->ctg_off = rd->ctg_off + offset;
				}
				if(offset + ctg2->len > ctg1->len){
					append_string(ctg1->seq, ctg2->seq->string + (ctg1->len - offset), ctg2->len - (ctg1->len - offset));
					ctg1->len = offset + ctg2->len;
				}
			} else {
				// ABCDEG
				//continue;
				ctg1->closed = 1;
				offset = off2 - off1;
				for(j=0;j<vec_size(ctg1->rids);j++){
					rid = gget_vec(ctg1->rids, j, uint32_t);
					rd  = get_vec_ref(ef->rds, rid);
					gpush_vec(ctg2->rids, rid, uint32_t);
					rd->ctg_id  = rd2->ctg_id;
					rd->ctg_off = rd->ctg_off + offset;
				}
				if(offset + ctg1->len > ctg2->len){
					append_string(ctg2->seq, ctg1->seq->string + ctg2->len - offset, ctg1->len - (ctg2->len - offset));
					ctg2->len = offset + ctg1->len;
				}
			}
		}
	}
}

void output_ef_ctgs(EF *ef, FILE *out){
	uint32_t i, j, cid;
	FContig *ctg;
	FRead *rd;
	cid = 0;
	fprintf(out, "E %u\n", ef->ef_id);
	for(i=0;i<vec_size(ef->ctgs);i++){
		ctg = gget_vec(ef->ctgs, i, FContig*);
		if(ctg->closed) continue;
		fprintf(out, "C %u\n", cid);
		cid ++;
		fprintf(out, "L %d\n", ctg->seq->size);
		fprintf(out, "S %s\n", ctg->seq->string);
		fprintf(out, "N %u\n", (uint32_t)vec_size(ctg->rids));
		fprintf(out, "R");
		for(j=0;j<vec_size(ctg->rids);j++){
			rd = get_vec_ref(ef->rds, gget_vec(ctg->rids, j, uint32_t));
			fprintf(out, " %u:%u:%u", rd->seq_id, rd->ctg_off, rd->rank);
		}
		fprintf(out, "\n//\n");
		fflush(out);
	}
}

void reset_ef(EF *ef, uint32_t ef_id, char *eseq, uint32_t rd_len, uint32_t min_ol, float min_sm){
	uint32_t i;
	ef->ef_id = ef_id;
	ef->min_ol = min_ol;
	ef->min_sm = min_sm;
	clear_vec(ef->rds);
	clear_vec(ef->ols);
	for(i=0;i<vec_size(ef->rps);i++){ put_pool_vec(ef, gget_vec(ef->rps, i, Vector*)); }
	clear_vec(ef->rps);
	for(i=0;i<vec_size(ef->ctgs);i++){ put_pool_ctg(ef, gget_vec(ef->ctgs, i, FContig*)); }
	clear_vec(ef->ctgs);
	clear_rhash(ef->index);
	clear_u64hash(ef->uniq);
	ef->eseq[rd_len] = 0;
	add_read2ef_core(ef, eseq, ef_id, rd_len, 0);
}

void free_ef(EF *ef){
	FContig *ctg;
	uint32_t i;
	free_vec(ef->rds);
	for(i=0;i<vec_size(ef->ctgs);i++){ put_pool_ctg(ef, gget_vec(ef->ctgs, i, FContig*)); }
	free_vec(ef->ctgs);
	for(i=0;i<vec_size(ef->rps);i++){ put_pool_vec(ef, gget_vec(ef->rps, i, Vector*)); }
	free_vec(ef->rps);
	free_vec(ef->ols);
	free_rhash(ef->index);
	free_u64hash(ef->uniq);
	for(i=0;i<vec_size(ef->pool_vec);i++){ free_vec(gget_vec(ef->pool_vec, i, Vector*)); }
	for(i=0;i<vec_size(ef->pool_ctg);i++){
		ctg = gget_vec(ef->pool_ctg, i, FContig*);
		free_vec(ctg->rids);
		free_string(ctg->seq);
		free(ctg);
	}
	free_vec(ef->pool_vec);
	free_vec(ef->pool_ctg);
	free(ef);
	ef = NULL;
}

uint32_t asm_ef(FileReader *in, FILE *out, uint32_t min_ol, float min_sm, uint32_t min_read, uint32_t max_read){
	EF *ef;
	uint32_t ret, ef_id, eid, rank, seqid;
	int n_col;
	ef = NULL;
	ret = 0;
	ef_id = 0;
	while((n_col = fread_table(in)) != -1){
		if(n_col == 0) continue;
		eid = atoi(get_col_str(in, 1));
		if(eid != ef_id){
			ef_id = eid;
			ret ++;
			reverse_dna(get_col_str(in, 2), get_col_len(in, 2));
			if(ef){
				if(vec_size(ef->rds) >= min_read){  //magic number 5
					align_reads_ef(ef);
					//print_alignments(ef);
					asm_ef_ctgs(ef);
					output_ef_ctgs(ef, out);
				}
				reset_ef(ef, ef_id, get_col_str(in, 2), get_col_len(in, 2), min_ol, min_sm);
			} else {
				ef = init_ef(ef_id, get_col_str(in, 2), get_col_len(in, 2), min_ol, min_sm);
			}
		}
		//rank  = atoi(get_col_str(in, 1));
		rank  = 1;
		seqid = atol(get_col_str(in, 0));
		if (vec_size(ef->rds) <= max_read) {  //magic number 200
			add_read2ef(ef, get_col_str(in, 3), seqid, get_col_len(in, 3), rank);
		}
	}
	if(ef && vec_size(ef->rds) >= min_read){
		align_reads_ef(ef);
		//print_alignments(ef);
		asm_ef_ctgs(ef);
		output_ef_ctgs(ef, out);
		free_ef(ef);
	}
	return ret;
}

int ef_usage(){
	printf(
"Local assemble fragments around restriction sites\n"
"Usage: ef [options]\n"
" -i <string> Input file [STDIN]\n"
" -o <string> Output file [STDOUT]\n"
" -l <int>    Minium length of overlap [5]\n"
" -s <float>  Minium similiarity of overlap [0.90]\n"
" -r <int>    Minium reads to execute assembly [5]\n"
" -R <int>    Maxium reads to execute assembly [200]\n"
"\n"
	);
	return 1;
}
