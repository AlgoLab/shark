CFLAGS	= -DNDEBUG -march=native -Wno-char-subscripts -Wall -O3 -std=c++14 -I. -I./include -fopenmp -fno-omit-frame-pointer
CXXFLAGS= ${CFLAGS}
LIBS = -L./lib -lz -lsdsl -ldivsufsort -ldivsufsort64 -ltbb

.PHONY: all clean

all: shark

shark: bloomfilter.o MurmurHash3.o main.o
	@echo "* Linking shark"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(LDFLAGS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

main.o: argument_parser.hpp bloomfilter.h BloomfilterFiller.hpp KmerBuilder.hpp FastaSplitter.hpp ReadAnalyzer.hpp ReadOutput.hpp kmer_utils.hpp

bloomfilter.o: bloomfilter.h

clean:
	rm -rf *.o
