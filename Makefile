CFLAGS	= -DNDEBUG -march=native -Wno-char-subscripts -Wall -Ofast -std=c++14 -I. -I./sdsl-lite/build/include -I./htslib/htslib -I./KMC -fopenmp -fomit-frame-pointer -foptimize-strlen -faggressive-loop-optimizations -funroll-loops -funsafe-math-optimizations 
#CFLAGS	= -g -Wno-char-subscripts -Wall -O0 -std=c++14 -I. -I./sdsl-lite/build/include -I./htslib/htslib -I./KMC -fopenmp
CXXFLAGS= ${CFLAGS}
LIBS = -L./sdsl-lite/build/lib -L./sdsl-lite/build/external/libdivsufsort/lib -L./htslib/ -lhts -lz -lsdsl -ldivsufsort -ldivsufsort64 -ltbb

.PHONY: all

all: read-filter

read-filter: main.o bloomfilter.o MurmurHash3.o bloomfilter.h\
						./KMC/kmc_api/kmc_file.o ./KMC/kmc_api/kmer_api.o ./KMC/kmc_api/mmer.o
	@echo "* Linking read-filter"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(LDFLAGS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o
	rm -rf KMC/kmc_api/*.o
