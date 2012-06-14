CC=gcc
CFLAGS= -O2 -W -Wall -finline-functions -D_FILE_OFFSET_BITS=64 
GLIBS=-lm
GENERIC_SRC= string.h bitvec.h file_reader.h hashset.h sort.h list.h dna.h heap.h stdaln.h vector.h

#all: rainbow rbasm rbmergetag ezmsim 
all: rainbow rbmergetag ezmsim rbasm

rainbow: $(GENERIC_SRC) file_reader.c rainbow.h asm_R2.c asm_R2.h vector.h mergectg.h cluster.c divide.c stdaln.c mergectg.c main.c
	$(CC) $(CFLAGS) $(GLIBS) -o $@ $^

rbasm: $(GENERIC_SRC) file_reader.c asm_R2.c rbasm_main.c
	$(CC) $(CFLAGS) $(GLIBS) -o $@ $^

rbmergetag: $(GENERIC_SRC) file_reader.c mergetag.c
	$(CC) $(CFLAGS) $(GLIBS) -o $@ $^

ezmsim: ezmsim.c
	$(CC) $(CFLAGS) $(GLIBS) -o $@ $^

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out rainbow rbasm rbmergetag *.exe

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out
