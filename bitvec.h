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
 
#ifndef __BIT_VEC_RJ_H
#define __BIT_VEC_RJ_H

#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define get_bit8(bits, idx) ((((bits)[(idx) >> 3]) >> ((idx) & 0x07)) & 0x01)
#define get_bit16(bits, idx) ((((bits)[(idx) >> 4]) >> ((idx) & 0x0F)) & 0x01)
#define get_bit32(bits, idx) ((((bits)[(idx) >> 5]) >> ((idx) & 0x1F)) & 0x01)
#define get_bit64(bits, idx) ((((bits)[(idx) >> 6]) >> ((idx) & 0x3F)) & 0x01)

#define get_2bit8(bits, idx) ((((bits)[(idx) >> 2]) >> (((idx) & 0x03) << 1)) & 0x03)
#define get_2bit16(bits, idx) ((((bits)[(idx) >> 3]) >> (((idx) & 0x07) << 1)) & 0x03)
#define get_2bit32(bits, idx) ((((bits)[(idx) >> 4]) >> (((idx) & 0x0F) << 1)) & 0x03)
#define get_2bit64(bits, idx) ((((bits)[(idx) >> 5]) >> (((idx) & 0x1F) << 1)) & 0x03)

typedef struct {
	uint64_t *bits;
	uint64_t n_bit;
	uint64_t n_cap;
	uint64_t *sums;
	uint64_t iter_idx;
} BitVec;

static inline BitVec* init_bitvec(uint64_t n_bit){
	BitVec *bitv;
	if(n_bit == 0) n_bit = 64 * 8;
	bitv = (BitVec*)malloc(sizeof(BitVec));
	bitv->n_bit = 0;
	bitv->n_cap = (((n_bit + 63) / 64) + 7) / 8 * 64 * 8;
	bitv->bits  = (uint64_t*)malloc(bitv->n_cap / 8);
	memset(bitv->bits, 0, bitv->n_cap / 8);
	bitv->sums = NULL;
	return bitv;
}

static inline void clear_bitvec(BitVec *bitv){ bitv->n_bit = 0; }

static inline void zeros_bitvec(BitVec *bitv){ memset(bitv->bits, 0, bitv->n_cap / 8); }

static inline void ones_bitvec(BitVec *bitv){ memset(bitv->bits, 0xFFU, bitv->n_cap / 8); }

static inline void flip_bitvec(BitVec *bitv, uint64_t idx){ bitv->bits[idx>>6] ^= 1LLU << (idx&0x3FU); }

static inline void one_bitvec(BitVec *bitv, uint64_t idx){ bitv->bits[idx>>6] |= 1LLU << (idx&0x3FU); }

static inline void zero_bitvec(BitVec *bitv, uint64_t idx){ bitv->bits[idx>>6] &= ~(1LLU << (idx&0x3FU)); }

static inline uint64_t get_bitvec(BitVec *bitv, uint64_t idx){ return (bitv->bits[idx>>6] >> (idx&0x3FU)) & 0x01LLU; }

static inline void encap_bitvec(BitVec *bitv, uint64_t num){
	if(bitv->n_bit + num < bitv->n_cap) return;
	while(bitv->n_bit + num >= bitv->n_cap){
		if(bitv->n_cap < 1024 * 1024 * 8){
			bitv->n_cap <<= 1;
		} else bitv->n_cap += 1024 * 1024 * 8;
	}
	bitv->bits = (uint64_t*)realloc(bitv->bits, bitv->n_cap / 8);
	memset(((void*)bitv->bits) + bitv->n_bit / 8, 0, (bitv->n_cap - bitv->n_bit) / 8);
}

static inline void one2bitvec(BitVec *bitv){ encap_bitvec(bitv, 1); one_bitvec(bitv, bitv->n_bit); bitv->n_bit ++; }

static inline void zero2bitvec(BitVec *bitv){ encap_bitvec(bitv, 1); zero_bitvec(bitv, bitv->n_bit); bitv->n_bit ++; }

static inline uint32_t count_ones_bit32(uint32_t v){
	v = v - ((v >> 1) & 0x55555555U);                        // reuse input as temporary
	v = (v & 0x33333333U) + ((v >> 2) & 0x33333333U);        // temp
	return (((v + (v >> 4)) & 0xF0F0F0FU) * 0x1010101U) >> 24; // count
}

#define ONES_STEP_4 0x1111111111111111ULL
#define ONES_STEP_8 0x0101010101010101ULL

static inline int count_ones_bit64(const uint64_t x){
	register uint64_t byte_sums = x - ((x & 0xa * ONES_STEP_4) >> 1);
	byte_sums = (byte_sums & 3 * ONES_STEP_4) + ((byte_sums >> 2) & 3 * ONES_STEP_4);
	byte_sums = (byte_sums + (byte_sums >> 4)) & 0x0f * ONES_STEP_8;
	return byte_sums * ONES_STEP_8 >> 56;
}

static inline void index_bitvec(BitVec *bitv){
	uint64_t i, s, t;
	if(bitv->sums) free(bitv->sums);
	bitv->sums = (uint64_t*)malloc((bitv->n_cap / 64 / 8 * 2 + 1) * 8);
	memset(bitv->sums, 0, bitv->n_cap / 64 / 8 * 2 * 8);
	t = 0;
	for(i=0;i<bitv->n_cap;i+=64*8){
		bitv->sums[((i>>6) >> 3) << 1] = t;
		s = 0;
		s += count_ones_bit64(bitv->bits[(i>>6)+0]);
		bitv->sums[(((i>>6) >> 3) << 1)+1] |= s << 0;
		s += count_ones_bit64(bitv->bits[(i>>6)+1]);
		bitv->sums[(((i>>6) >> 3) << 1)+1] |= s << 9;
		s += count_ones_bit64(bitv->bits[(i>>6)+2]);
		bitv->sums[(((i>>6) >> 3) << 1)+1] |= s << 18;
		s += count_ones_bit64(bitv->bits[(i>>6)+3]);
		bitv->sums[(((i>>6) >> 3) << 1)+1] |= s << 27;
		s += count_ones_bit64(bitv->bits[(i>>6)+4]);
		bitv->sums[(((i>>6) >> 3) << 1)+1] |= s << 36;
		s += count_ones_bit64(bitv->bits[(i>>6)+5]);
		bitv->sums[(((i>>6) >> 3) << 1)+1] |= s << 45;
		s += count_ones_bit64(bitv->bits[(i>>6)+6]);
		bitv->sums[(((i>>6) >> 3) << 1)+1] |= s << 54;
		s += count_ones_bit64(bitv->bits[(i>>6)+7]);
		t += s;
	}
	bitv->sums[((i>>6) >> 3) << 1] = t;
}

static inline uint64_t rank_bitvec(BitVec *bitv, uint64_t idx){
	uint64_t p, s, sum;
	p = (idx>>6)>>3;
	s = (idx >> 6) & 0x07U;
	sum = bitv->sums[p<<1];
	if(s) sum += (bitv->sums[(p<<1)+1] >> (9 * (s - 1))) & 0x1FFU;
	if(idx & 0x3FU) sum += count_ones_bit64(bitv->bits[idx>>6]<<(64-(idx&0x3FU)));
	return sum;
}

static inline void begin_iter_bitvec(BitVec *bitv){ bitv->iter_idx = 0; }

static inline uint64_t iter_bitvec(BitVec *bitv){
	while(bitv->iter_idx < bitv->n_cap){
		if((bitv->iter_idx & 0x1FFU) == 0 && bitv->sums[(((bitv->iter_idx>>6)>>3)<<1)] == bitv->sums[((((bitv->iter_idx>>6)>>3)+1)<<1)]){
			bitv->iter_idx += 64 * 8;
			continue;
		}
		if((bitv->bits[bitv->iter_idx>>6] >> (bitv->iter_idx & 0x3FU)) == 0){
			bitv->iter_idx = ((bitv->iter_idx >> 6) + 1) << 6;
			continue;
		}
		if((bitv->bits[bitv->iter_idx>>6] >> (bitv->iter_idx&0x3FU) & 0x01U)){
			bitv->iter_idx ++;
			return bitv->iter_idx - 1;
		} else {
			bitv->iter_idx ++;
		}
	}
	return 0xFFFFFFFFFFFFFFFFLLU;
}

static inline void free_bitvec(BitVec *bitv){
	free(bitv->bits);
	if(bitv->sums) free(bitv->sums);
	free(bitv);
}

#endif
