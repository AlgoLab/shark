CFLAGS	= -g -Wno-char-subscripts -Wall -O3 -std=c++14 -I. -I./sdsl-lite/build/include -I./htslib/htslib -I./KMC -fopenmp
CXXFLAGS= ${CFLAGS}
LIBS = -L./sdsl-lite/build/lib -L./sdsl-lite/build/external/libdivsufsort/lib -L./htslib/ -lhts -lz -lsdsl -ldivsufsort -ldivsufsort64

.PHONY: all

all: read-filter

read-filter: main.o MurmurHash3.o bloomfilter.h\
						./KMC/kmc_api/kmc_file.o ./KMC/kmc_api/kmer_api.o ./KMC/kmc_api/mmer.o
	@echo "* Linking read-filter"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(LDFLAGS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o
	rm -rf KMC/kmc_api/*.o
