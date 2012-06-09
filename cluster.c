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

SeqDB* load_seqdb(FileReader *fr, int is_fq, int fix_rd_len){
	SeqDB *sdb;
	Sequence *seq;
	u64list *seqs;
	uint64_t offset;
	uint8_t len;
	sdb = malloc(sizeof(SeqDB));
	sdb->n_rd = 0;
	sdb->rd_len = 0;
	sdb->max_rd_len = 0;
	if(!fix_rd_len){
		sdb->seqoffs = init_u64list(1024);
		sdb->seqlens = init_u8list(1024);
	} else {
		sdb->rd_len = fix_rd_len;
		sdb->max_rd_len = fix_rd_len;
		sdb->seqoffs = NULL;
		sdb->seqlens = NULL;
	}
	seq = NULL;
	seqs = init_u64list(1024);
	offset = 0;
	while(is_fq? fread_fastq_adv(&seq, fr, 5) : fread_fasta_adv(&seq, fr, 1)){
		sdb->n_rd ++;
		len = seq->seq.size;
		if(fix_rd_len){
			if(len < sdb->rd_len){
				continue;
			} else if(len > sdb->rd_len){
				len = sdb->rd_len;
			}
		} else {
			if(sdb->max_rd_len < len) sdb->max_rd_len = len;
			push_u64list(sdb->seqoffs, offset);
			push_u8list(sdb->seqlens, len);
		}
		encap_u64list(seqs, (offset + len + 31) / 32 + 2);
		seq2bits(seqs->buffer, offset, seq->seq.string, len);
		offset += len;
		if((sdb->n_rd & 0xFFF) == 0){ fprintf(stderr, "\r %u k reads  ", (sdb->n_rd >> 10)); fflush(stderr); }
	}
	sdb->seqs = seqs->buffer;
	free(seqs);
	return sdb;
}

uint8_t prepare_seq_seqdb(SeqDB *sdb, uint32_t rid, uint64_t *seqs){
	uint64_t off;
	uint32_t i, j;
	uint8_t len;
	if(sdb->rd_len){
		off = sdb->rd_len * rid;
		len = sdb->rd_len;
	} else {
		off = get_u64list(sdb->seqoffs, rid);
		len = get_u8list(sdb->seqlens, rid);
	}
	j = 0;
	for(i=0;i+32<len;i+=32){
		seqs[j++] = sub32seqbits(sdb->seqs, off + i);
	}
	if(i < len){
		seqs[j] = sub32seqbits(sdb->seqs, off + i);
		seqs[j] >>= (32 - (len - i)) << 1;
	}
	return len;
}

uint8_t cal_2seq_mm_core(uint64_t *seq1, uint64_t *seq2, uint8_t len1, uint8_t len2){
	uint32_t i, len, mm;
	len = (len1 < len2)? len1 : len2;
	mm = 0;
	for(i=0;i<2&&i<len;i+=32){
		mm += count_ones_bit64(dna_xor2ones(seq1[i >> 5] ^ seq2[i >> 5]));
	}
	return mm;
}

uint32_t linking_core(Cluster *cluster, uint32_t seqid, uint64_t *seq, uint32_t seqlen){
	kmer_t K, *k;
	uint32_t j, off, c;
	uint32_t link;
	int exists;
	if(seqlen < (cluster->idxs[1] + 1) * KMER_SIZE) return seqid;
	K.kmer1 = 0;
	K.kmer2 = 0;
	K.seqid = seqid;
	{
		off = cluster->idxs[0] * KMER_SIZE;
		for(j=0;j<KMER_SIZE;j++){
			c = bits2bit(seq, off + j);
			K.kmer1 = (K.kmer1 << 2) | c;
		}
	}
	{
		off = cluster->idxs[1] * KMER_SIZE;
		for(j=0;j<KMER_SIZE;j++){
			c = bits2bit(seq, off + j);
			K.kmer2 = (K.kmer2 << 2) | c;
		}
	}
	/*
	for(i=0;i<2;i++){
		off = cluster->idxs[i] * KMER_SIZE;
		for(j=0;j<KMER_SIZE;j++){
			c = bits2bit(seq, off + j);
			K.kmer1 = (K.kmer1 << 2) | c;
		}
	}
	*/
	k = prepare_khash(cluster->index, K, &exists);
	if(exists){
		link = k->seqid;
	} else {
		k->kmer1 = K.kmer1;
		k->kmer2 = K.kmer2;
		link = seqid;
	}
	k->seqid = seqid;
	return link;
}

void tracing_core(Cluster *cluster, uint32_t bt){
	uint32_t next;
	clear_u32list(cluster->bts);
	push_u32list(cluster->bts, bt);
	while(1){
		one_bitvec(cluster->flags, bt);
		next = get_u32list(cluster->links, bt);
		if(next == bt) break;
		push_u32list(cluster->bts, next);
		bt = next;
	}
}

inline int cmp_sbt(const void *e1, const void *e2){
	SBT *t1, *t2;
	uint32_t i, len;
	t1 = (SBT*)e1;
	t2 = (SBT*)e2;
	len = (t1->len < t2->len)? t2->len : t1->len;
	len = (31 + len) / 32;
	for(i=0;i<len;i++){
		if(t1->seq[i] < t2->seq[i]) return -1;
		if(t1->seq[i] > t2->seq[i]) return 1;
	}
	return 0;
}

uint32_t sorting_core(Cluster *cluster){
	SBT *sbt1, *sbt2;
	uint32_t ret, i, *gid1, *gid2;
	ret = 0;
	clear_sbtv(cluster->sbts);
	for(i=0;i<count_u32list(cluster->bts);i++){
		sbt1 = next_ref_sbtv(cluster->sbts);
		sbt1->bt  = get_u32list(cluster->bts, i);
		sbt1->len = prepare_seq_seqdb(cluster->sdb, sbt1->bt, sbt1->seq);
	}
	qsort(as_array_sbtv(cluster->sbts), count_sbtv(cluster->sbts), sizeof(SBT), cmp_sbt);
	sbt1 = ref_sbtv(cluster->sbts, 0);
	gid1 = ref_u32list(cluster->gids, sbt1->bt);
	*gid1 = get_u32list(cluster->gid_map, *gid1);
	for(i=1;i<count_sbtv(cluster->sbts);i++){
		sbt2 = ref_sbtv(cluster->sbts, i);
		gid2 = ref_u32list(cluster->gids, sbt2->bt);
		*gid2 = get_u32list(cluster->gid_map, *gid2);
		if(cmp_sbt((const void *)sbt1, (const void *)sbt2) == 0){
			ret ++;
			if(*gid1){
				if(*gid2){
					if(*gid1 < *gid2){
						set_u32list(cluster->gid_map, *gid2, *gid1);
						*gid2 = *gid1;
					} else if(*gid1 > *gid2){
						set_u32list(cluster->gid_map, *gid1, *gid2);
						*gid1 = *gid2;
					}
				} else {
					*gid2 = *gid1;
				}
			} else {
				if(*gid2){
					*gid1 = *gid2;
				} else {
					push_u32list(cluster->gid_map, ++cluster->gidoff);
					*gid1 = cluster->gidoff;
					*gid2 = *gid1;
				}
			}
		} else {
			sbt1 = sbt2;
			gid1 = gid2;
		}
	}
	return ret;
}

uint32_t alning_core(Cluster *cluster){
	uint32_t idx1, idx2, ret, m, n, mm, *gid1, *gid2;
	uint8_t len1, len2;
	ret = 0;
	for(m=0;m+1<count_u32list(cluster->bts);m++){
		idx1 = get_u32list(cluster->bts, m);
		len1 = prepare_seq_seqdb(cluster->sdb, idx1, cluster->seq1);
		gid1 = ref_u32list(cluster->gids, idx1);
		*gid1 = get_u32list(cluster->gid_map, *gid1);
		for(n=m+1;n<count_u32list(cluster->bts);n++){
			idx2 = get_u32list(cluster->bts, n);
			len2 = prepare_seq_seqdb(cluster->sdb, idx2, cluster->seq2);
			gid2 = ref_u32list(cluster->gids, idx2);
			*gid2 = get_u32list(cluster->gid_map, *gid2);
			if(*gid1 && *gid1 == *gid2) continue;
			mm = cal_2seq_mm_core(cluster->seq1, cluster->seq2, len1, len2);
			if(mm > cluster->max_mm) continue;
			ret ++;
			if(*gid1){
				if(*gid2){
					if(*gid1 < *gid2){
						set_u32list(cluster->gid_map, *gid2, *gid1);
						*gid2 = *gid1;
					} else {
						set_u32list(cluster->gid_map, *gid1, *gid2);
						*gid1 = *gid2;
					}
				} else {
					*gid2 = *gid1;
				}
			} else {
				if(*gid2){
					*gid1 = *gid2;
				} else {
					push_u32list(cluster->gid_map, ++cluster->gidoff);
					*gid1 = cluster->gidoff;
					*gid2 = *gid1;
				}
			}
		}
	}
	return ret;
}

void indexing_cluster(Cluster *cluster, FileReader *fr, int is_fq, int fix_rd_len){
	uint64_t cnt;
	uint32_t i, seqid, len1, bt, *gid, max_rd_len;
	clock_t t0, t1;
	t0 = clock();
	fprintf(stderr, "Load pair1\n"); fflush(stderr);
	cluster->sdb = load_seqdb(fr, is_fq, fix_rd_len);
	max_rd_len = cluster->sdb->max_rd_len;
	t1 = clock();
	fprintf(stderr, "\r %u reads, %.2f secs [OK]\n", cluster->sdb->n_rd, ((double)t1 - t0) / CLOCKS_PER_SEC); fflush(stderr);
	for(cluster->idxs[0]=0;cluster->idxs[0]<KMER_NUM;cluster->idxs[0]++){
		for(cluster->idxs[1]=cluster->idxs[0]+1;cluster->idxs[1]<KMER_NUM;cluster->idxs[1]++){
			fprintf(stderr, "Iterating %u/%u %u/%u\n", cluster->idxs[0], KMER_NUM, cluster->idxs[1], KMER_NUM); fflush(stderr);
			if(max_rd_len && (cluster->idxs[1] + 1) * KMER_SIZE > max_rd_len){
				fprintf(stderr, "- Skip\n"); fflush(stderr);
				continue;
			}
			clear_u32list(cluster->links);
			clear_bitvec(cluster->flags);
			cluster->index = init_khash(1023);
			t0 = clock();
			fprintf(stderr, "- Linking\n0 k"); fflush(stderr);
			for(seqid=0;seqid<cluster->sdb->n_rd;seqid++){
				len1 = prepare_seq_seqdb(cluster->sdb, seqid, cluster->seq1);
				push_u32list(cluster->links, linking_core(cluster, seqid, cluster->seq1, len1));
				if((seqid & 0xFFFU) == 0){ fprintf(stderr, "\r %u k", (seqid>>10)); fflush(stderr); }
			}
			free_khash(cluster->index);
			t1 = clock();
			fprintf(stderr, "\r %u k reads, %0.2f secs [OK]\n", (seqid >> 10), ((double)t1 - t0) / CLOCKS_PER_SEC); fflush(stderr);
			if(cluster->max_seqid == 0){
				cluster->max_seqid = seqid;
				clear_u32list(cluster->gids);
				for(i=0;i<seqid;i++){
					push_u32list(cluster->gids, 0);
				}
			}
			t0 = clock();
			cluster->gid_map = init_u32list(cluster->gidoff + 1024);
			for(i=0;i<=cluster->gidoff;i++){ push_u32list(cluster->gid_map, i); }
			fprintf(stderr, "- Aligning (%u mismatches)\n", cluster->max_mm); fflush(stderr);
			t0 = clock();
			encap_bitvec(cluster->flags, count_u32list(cluster->links));
			zeros_bitvec(cluster->flags);
			cnt = 0;
			for(i=cluster->max_seqid;i;i--){
				bt = i - 1;
				if(get_bitvec(cluster->flags, bt) == 0){
					tracing_core(cluster, bt);
					if(count_u32list(cluster->bts) == 1){
					} else if(count_u32list(cluster->bts) >= cluster->exact_limit){
						cnt += sorting_core(cluster);
					} else {
						cnt += alning_core(cluster);
					}
				}
				if((i&0xFFFF) == 0){ fprintf(stderr, "\r hits: %u", (unsigned)cnt); fflush(stderr); }
			}
			t1 = clock();
			fprintf(stderr, "\r hits: %u, %0.2f secs", (unsigned)cnt, ((double)t1 - t0) / CLOCKS_PER_SEC);
			fprintf(stderr, " [OK]\n"); fflush(stderr);
			fprintf(stderr, "- Translating group ids "); fflush(stderr);
			t0 = clock();
			cnt = 0;
			for(i=0;i<count_u32list(cluster->gids);i++){
				gid = ref_u32list(cluster->gids, i);
				if(get_u32list(cluster->gid_map, *gid) != *gid){ cnt ++;  *gid = get_u32list(cluster->gid_map, *gid); }
			}
			free_u32list(cluster->gid_map);
			t1 = clock();
			fprintf(stderr, " %llu in %0.2f secs [OK]\n", (unsigned long long)cnt, ((double)t1 - t0) / CLOCKS_PER_SEC); fflush(stderr);
		}
	}
}

void clustering(Cluster *cluster, FileReader *fr2, int is_fq2, int fix_rd_len, FILE *out){
	uint32_t i, seqid, gid;
	char seq1[256], seq2[256];
	u32list *rids;
	clock_t t0, t1;
	t0 = clock();
	fprintf(stderr, "sorting groups ... "); fflush(stderr);
	rids = init_u32list(cluster->sdb->n_rd);
	for(i=0;i<cluster->sdb->n_rd;i++) push_u32list(rids, i);
	sort_array(rids->buffer, rids->size, uint32_t, (((int64_t)cluster->gids->buffer[a]) - ((int64_t)cluster->gids->buffer[b])));
	t1 = clock();
	fprintf(stderr, " %.2f secs [OK]\n", ((double)t1 - t0) / CLOCKS_PER_SEC); fflush(stderr);
	if(fr2){
		t0 = clock();
		fprintf(stderr, "Load pair2\n"); fflush(stderr);
		cluster->sdb2 = load_seqdb(fr2, is_fq2, fix_rd_len);
		t1 = clock();
		fprintf(stderr, "\r %u reads, %.2f secs [OK]\n", cluster->sdb2->n_rd, ((double)t1 - t0) / CLOCKS_PER_SEC); fflush(stderr);
		if(cluster->sdb->n_rd != cluster->sdb2->n_rd){
			fprintf(stderr, " Pair2 didn't match Pair2 -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
		}
	}
	for(i=0;i<cluster->sdb->n_rd;i++){
		seqid = get_u32list(rids, i);
		gid = get_u32list(cluster->gids, seqid);
		if(gid == 0) continue;
		if(cluster->sdb->rd_len){
			bits2seq(seq1, cluster->sdb->seqs, seqid * cluster->sdb->rd_len, cluster->sdb->rd_len);
		} else {
			bits2seq(seq1, cluster->sdb->seqs, get_u64list(cluster->sdb->seqoffs, seqid), get_u8list(cluster->sdb->seqlens, seqid));
		}
		if(cluster->sdb2){
			if(cluster->sdb2->rd_len){
				bits2seq(seq2, cluster->sdb2->seqs, seqid * cluster->sdb2->rd_len, cluster->sdb2->rd_len);
			} else {
				bits2seq(seq2, cluster->sdb2->seqs, get_u64list(cluster->sdb2->seqoffs, seqid), get_u8list(cluster->sdb2->seqlens, seqid));
			}
		} else {
			seq2[0] = 'N';
			seq2[1] = '\0';
		}
		fprintf(out, "%u\t%u\t%s\t%s\n", seqid, gid, seq1, seq2);
		if((i & 0xFFFU) == 0){ fprintf(stderr, "\r output %u k seq    ", (unsigned)(i >> 10)); fflush(stderr); }
	}
	fprintf(stderr, "\r output %u k seq, %.2f secs [OK]\n", (unsigned)(i >> 10), ((double)t1 - t0) / CLOCKS_PER_SEC); fflush(stderr);
	free_u32list(rids);
}

Cluster* init_cluster(uint32_t max_mm, uint32_t exact_limit){
	Cluster *cluster;
	cluster = malloc(sizeof(Cluster));
	cluster->sdb  = NULL;
	cluster->sdb2 = NULL;
	cluster->gidoff  = 0;
	cluster->max_seqid = 0;
	cluster->max_mm  = max_mm;
	cluster->exact_limit = exact_limit;
	cluster->max_pair_len = 2 * KMER_SIZE;
	cluster->idxs[0] = 0;
	cluster->idxs[1] = 0;
	cluster->index = NULL;
	cluster->flags = init_bitvec(1024);
	cluster->links = init_u32list(1024);
	cluster->bts   = init_u32list(64);
	cluster->sbts  = init_sbtv(1024);
	cluster->gids  = init_u32list(1024);
	cluster->gid_map = NULL;
	return cluster;
}

void free_cluster(Cluster *cluster){
	if(cluster->sdb){
		free(cluster->sdb->seqs);
		if(cluster->sdb->rd_len == 0){
			free_u64list(cluster->sdb->seqoffs);
			free_u8list(cluster->sdb->seqlens);
		}
		free(cluster->sdb);
	}
	if(cluster->sdb2){
		free(cluster->sdb2->seqs);
		if(cluster->sdb2->rd_len == 0){
			free_u64list(cluster->sdb2->seqoffs);
			free_u8list(cluster->sdb2->seqlens);
		}
		free(cluster->sdb2);
	}
	free_bitvec(cluster->flags);
	free_u32list(cluster->links);
	free_u32list(cluster->gids);
	free_u32list(cluster->bts);
	free_sbtv(cluster->sbts);
	free(cluster);
}
