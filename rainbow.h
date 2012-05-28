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
 
#ifndef __RAINBOW_RJ_H
#define __RAINBOW_RJ_H

#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include "bitvec.h"
#include "hashset.h"
#include "list.h"
#include "sort.h"
#include "dna.h"
#include "file_reader.h"
#include "string.h"
//#include "mergecontig.h"
#include "mergectg.h"


#define KMER_SIZE	8
#define KMER_NUM	4

typedef struct {
	uint64_t kmer:32, seqid:32;
} kmer_t;

#define kmer_hashcode(k) u32hashcode((k).kmer)
#define kmer_equals(k1, k2) ((k1).kmer == (k2).kmer)
define_hashset(khash, kmer_t, kmer_hashcode, kmer_equals);

typedef struct {
	uint64_t *seqs;
	uint32_t n_rd;
	uint8_t  rd_len, max_rd_len;
	u64list *seqoffs;
	u8list  *seqlens;
} SeqDB;

typedef struct {
	uint32_t bt;
	uint32_t len;
	uint64_t seq[8];
} SBT;

define_list(sbtv, SBT);

typedef struct {
	SeqDB    *sdb, *sdb2;
	uint64_t seq1[10], seq2[10];
	uint32_t gidoff;
	uint32_t max_seqid;
	uint32_t max_pair_len;
	uint32_t max_mm;
	uint32_t exact_limit;
	uint32_t idxs[2];
	khash *index;
	u32list *links;
	BitVec  *flags;
	//uuhash *gid_map;
	u32list *gid_map;
	u32list *gids;
	u32list *bts;
	sbtv    *sbts;
} Cluster;

Cluster* init_cluster(uint32_t max_mm, uint32_t exact_limit);
void indexing_cluster(Cluster *cluster, FileReader *fr1, int is_fq, int fix_rd_len);
void clustering(Cluster *cluster, FileReader *fr2, int is_fq, int fix_rd_len, FILE *out);
void free_cluster(Cluster *cluster);

typedef struct {
	uint32_t seqid, seqoff, seqlen1:10, seqlen2:10, rank:6, revsed:6;
} ReadInfo;

define_list(rilist, ReadInfo);

define_list(u32slist, u32list*);

typedef struct {
	uint32_t gidoff;
	rilist *rds;
	u8list *seqs;
	u32slist *grps, *cache;
	u64list *markers;
	u32list *deps;
	u32list *gids;
	uint32_t n_col;
	uint32_t k_allele, K_allele;
	float min_freq;
} Div;

Div* init_div(uint32_t k_allele, uint32_t K_allele, float min_freq);
uint32_t div_reads(Div *div, FileReader *fr, FILE *out);
void reset_div(Div *div);
void free_div(Div *div);


#endif
