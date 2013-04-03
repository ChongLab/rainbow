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
 
#include "file_reader.h"
#include "dna.h"
#include "list.h"
#include "sort.h"
#include <unistd.h>

typedef struct {
	uint32_t rid, cid, gid, rank, off1, len1, off2, len2;
} REC;

define_list(recv, REC);

typedef struct {
	uint32_t gid, off, len;
	uint32_t cns_off, cns_len;
} Block;

define_list(blockv, Block);

void consensus(recv *divs, String *seqs, uint32_t beg, uint32_t end, uint32_t *cns_off, uint32_t *cns_len){
	REC *r;
	uint32_t i, j, len, acgtn[5], ref;
	len = 0;
	for(i=beg;i<end;i++){
		r = ref_recv(divs, i);
		if(r->len1 > len) len = r->len1;
	}
	*cns_off = seqs->size;
	*cns_len = len;
	for(i=0;i<len;i++){
		acgtn[0] = 0;
		acgtn[1] = 0;
		acgtn[2] = 0;
		acgtn[3] = 0;
		acgtn[4] = 0;
		for(j=beg;j<end;j++){
			r = ref_recv(divs, j);
			if(r->len1 <= i) continue;
			acgtn[base_bit_table[(int)seqs->string[r->off1 + i]]] ++;
		}
		ref = 0;
		for(j=1;j<4;j++){
			if(acgtn[j] > acgtn[ref]) ref = j;
		}
		add_char_string(seqs, bit_base_table[ref]);
	}
	add_char_string(seqs, '\0');
}

uint32_t cal_mm(String *seqs, Block *b1, Block *b2){
	uint32_t mm, i, len;
	mm = 0;
	len = b1->cns_len;
	if(len > b2->cns_len) len = b2->cns_len;
	for(i=0;i<len;i++){
		if(seqs->string[b1->cns_off + i] != seqs->string[b2->cns_off + i]) mm ++;
	}
	if(len < b1->cns_len) mm += b1->cns_len - len;
	else if(len < b2->cns_len) mm += b2->cns_len - len;
	return mm;
}

void merge_core(recv *divs, String *seqs, uint32_t max_mm, int task, blockv *blocks, FILE *out){
	Block *b;
	REC *r;
	uint32_t i, j, gid, beg, mm;
	beg = 0;
	gid = 0;
	clear_blockv(blocks);
	for(i=0;;i++){
		if(i < count_recv(divs) && ref_recv(divs, i)->gid == gid) continue;
		if(i > beg){
			b = next_ref_blockv(blocks);
			b->gid = gid;
			b->off = beg;
			b->len = i - beg;
			consensus(divs, seqs, beg, i, &b->cns_off, &b->cns_len);
			if(task == 1){
				fprintf(out, "%u\t%s\n", gid, seqs->string + b->cns_off);
			}
		}
		if(i == count_recv(divs)) break;
		beg = i;
		gid = ref_recv(divs, i)->gid;
	}
	if(task == 1) return;
	for(i=0;i+1<count_blockv(blocks);i++){
		for(j=i+1;j<count_blockv(blocks);j++){
			mm = cal_mm(seqs, ref_blockv(blocks, i), ref_blockv(blocks, j));
			if(mm <= max_mm) ref_blockv(blocks, j)->gid = ref_blockv(blocks, i)->gid;
		}
	}
	sort_array(blocks->buffer, blocks->size, Block, ((a.gid == b.gid)? 0 : ((a.gid < b.gid)? -1 : 1)));
	for(i=0;i<blocks->size;i++){
		b = ref_blockv(blocks, i);
		for(j=0;j<b->len;j++){
			r = ref_recv(divs, b->off + j);
			fprintf(out, "%u\t%u\t%s\t%s\t%u\n", r->rid, b->gid, seqs->string + r->off1, seqs->string + r->off2, r->cid);
		}
	}
}

int usage(){
	printf(
			"Usage: rbmergetag [options]\n"
			"Options:\n"
			" -i <string>    Input file name [stdin]\n"
			" -o <string>    Output file name [stdout]\n"
			" -j <cns|merge> Job type, cns: consensus, merge: merging, [merge]\n"
			" -m <int>       Maximum mismatches to merge two groups [1]\n"
			" -h             Show this document\n"
		  );
	return 1;
}

int main(int argc, char **argv){
	FileReader *fr;
	recv *divs;
	REC *r;
	blockv *blocks;
	FILE *out;
	String *seqs;
	char *inf, *ouf;
	uint32_t cid, max_mm;
	int n, c, task;
	max_mm = 1;
	task = 2;
	inf = NULL;
	ouf = NULL;
	while((c = getopt(argc, argv, "hi:o:j:m:")) != -1){
		switch(c){
			case 'i': inf = optarg; break;
			case 'o': ouf = optarg; break;
			case 'j': task = (strcasecmp(optarg, "cns") == 0)? 1 : 2; break;
			case 'm': max_mm = atoi(optarg); break;
			default: return usage();
		}
	}
	if(inf == NULL){ fr = stdin_filereader(); }
	else if((fr = fopen_filereader(inf)) == NULL){
		fprintf(stderr, "Cannot read '%s'\n", inf);
		return 1;
	}
	if(ouf == NULL){ out = stdout; }
	else if((out = fopen(ouf, "w")) == NULL){
		fprintf(stderr, "Cannot write'%s'\n", ouf);
		return 1;
	}
	divs = init_recv(1024);
	seqs = init_string(1024);
	blocks = init_blockv(12);
	cid = 0;
	while(1){
		n = fread_table(fr);
		if(n == -1 || (uint32_t)atoll(get_col_str(fr, 2)) != cid){
			if(count_recv(divs)){ merge_core(divs, seqs, max_mm, task, blocks, out); }
			clear_string(seqs);
			clear_recv(divs);
			if(n == -1) break;
			cid = atoll(get_col_str(fr, 2));
		}
		{
			r = next_ref_recv(divs);
			r->rid  = atoll(get_col_str(fr, 0));
			r->rank = 1;
			r->cid  = atoll(get_col_str(fr, 1));
			r->gid  = atoll(get_col_str(fr, 2));
			r->off1 = seqs->size;
			r->len1 = get_col_len(fr, 3);
			append_string(seqs, get_col_str(fr, 3), get_col_len(fr, 3));
			add_char_string(seqs, '\0');
			r->off2 = seqs->size;
			r->len2 = get_col_len(fr, 4);
			append_string(seqs, get_col_str(fr, 4), get_col_len(fr, 4));
			add_char_string(seqs, '\0');
		}
	}
	free_recv(divs);
	free_string(seqs);
	free_blockv(blocks);
	fclose_filereader(fr);
	if(ouf) fclose(out);
	return 0;
}
