CC=gcc
CFLAGS= -W -O2 -Wall -Wno-self-assign -Wno-unused-function
DFLAGS= -D_FILE_OFFSET_BITS=64
GLIBS=-lm
GENERIC_SRC= string.h bitvec.h file_reader.h hashset.h sort.h list.h dna.h heap.h stdaln.h vector.h

.SUFFIXES:.c .o

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

#all: rainbow rbasm rbmergetag ezmsim 
all: rainbow ezmsim rbasm

rainbow: main.o  divide.o file_reader.o asm_R2.o mergectg.o cluster.o
	$(CC) $(CFLAGS) -o $@ $^ $(GLIBS) 

rbasm: asm_R2.o rbasm_main.o file_reader.o
	$(CC) $(CFLAGS) -o $@ $^ $(GLIBS)

ezmsim: ezmsim.o
	$(CC) $(CFLAGS) -o $@ $^ $(GLIBS)

asm_R2.o: asm_R2.c asm_R2.h string.h vector.h hashset.h file_reader.h \
  dna.h
cluster.o: cluster.c rainbow.h bitvec.h hashset.h list.h sort.h dna.h \
  file_reader.h string.h vector.h mergectg.h stdaln.h asm_R2.h \
  bloom_filter.h
divide.o: divide.c rainbow.h bitvec.h hashset.h list.h sort.h dna.h \
  file_reader.h string.h vector.h mergectg.h stdaln.h asm_R2.h \
  bloom_filter.h
ezmsim.o: ezmsim.c
file_reader.o: file_reader.c file_reader.h string.h vector.h
main.o: main.c rainbow.h bitvec.h hashset.h list.h sort.h dna.h \
  file_reader.h string.h vector.h mergectg.h stdaln.h asm_R2.h \
  bloom_filter.h
mergectg.o: mergectg.c mergectg.h list.h sort.h file_reader.h string.h \
  vector.h hashset.h stdaln.h asm_R2.h dna.h bloom_filter.h bitvec.h \
  rainbow.h
rbasm_main.o: rbasm_main.c asm_R2.h string.h vector.h hashset.h \
  file_reader.h dna.h
stdaln.o: stdaln.c stdaln.h

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out rainbow rbasm ezmsim rbmergetag *.exe

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out
