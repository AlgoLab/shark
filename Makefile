CFLAGS	= -DNDEBUG -UDEBUG -Wno-char-subscripts -Wall -O3 -std=c++14 -I. -I./include
CXXFLAGS= ${CFLAGS}
LIBS = -L./lib -lz -lsdsl -pthread

.PHONY: all clean

all: shark

shark: main.o
	@echo "* Linking shark"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(LDFLAGS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

main.o: common.hpp argument_parser.hpp bloomfilter.h KmerBuilder.hpp FastaSplitter.hpp FastqSplitter.hpp ReadAnalyzer.hpp ReadOutput.hpp kmer_utils.hpp small_vector.hpp

clean:
	rm -rf *.o
