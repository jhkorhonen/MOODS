CXXFLAGS = -Wall -O3 -fPIC -fforce-addr -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4
OBJECT_FILES = pssm_algorithms.o find_pssm_dna.o

all: find_pssm_dna libpssm.a

find_pssm_dna: libpssm.a find_pssm_dna.o
	g++ find_pssm_dna.o -L. -lpssm -o find_pssm_dna
	
libpssm.a: pssm_algorithms.o
	ar rcs libpssm.a pssm_algorithms.o
	
find_pssm_dna.o: find_pssm_dna.cpp seq_buffer.h
	g++ $(CXXFLAGS) -c find_pssm_dna.cpp -o find_pssm_dna.o

pssm_algorithms.o: pssm_algorithms.cpp pssm_algorithms.hpp seq_buffer.h
	g++ $(CXXFLAGS) -c pssm_algorithms.cpp -o pssm_algorithms.o

clean:
	rm -f  libpssm.a find_pssm_dna $(OBJECT_FILES)
