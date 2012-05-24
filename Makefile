CC=gcc
CFLAGS= -W -O2 -Wall -finline-functions -D_FILE_OFFSET_BITS=64
GLIBS=-lm
GENERIC_SRC= string.h bitvec.h file_reader.h hashset.h sort.h list.h dna.h heap.h stdaln.h

all: rainbow rbasm rbmergetag ezmsim

rainbow: $(GENERIC_SRC) file_reader.c rainbow.h mergecontig.h cluster.c divide.c stdaln.c mergecontig.c main.c
	$(CC) $(CFLAGS) $(GLIBS) -o $@ $^

rbasm: $(GENERIC_SRC) file_reader.c asm_R2.c
	$(CC) $(CFLAGS) $(GLIBS) -o $@ $^

rbmergetag: $(GENERIC_SRC) file_reader.c mergetag.c
	$(CC) $(CFLAGS) $(GLIBS) -o $@ $^

ezmsim: ezmsim.c
	$(CC) $(CFLAGS) $(GLIBS) -o $@ $^

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out rainbow rbasm rbmergetag *.exe

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out
