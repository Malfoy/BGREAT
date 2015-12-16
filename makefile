# CC=/usr/bin/g++
CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread
LDFLAGS=-pthread


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4   -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=bgreat  getLargeUnitigs

all: $(EXEC)




getLargeUnitigs: getLargeUnitigs.o
	$(CC) -o $@ $^ $(LDFLAGS)

aligner.o: aligner.cpp aligner.h
	$(CC) -o $@ -c $< $(CFLAGS)

bgreat: bgreat.o   aligner.o
	$(CC) -o $@ $^ $(LDFLAGS)

getBigUnitigs.o: getBigUnitigs.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

bgreat.o: bgreat.cpp  aligner.h
	$(CC) -o $@ -c $< $(CFLAGS)





clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
