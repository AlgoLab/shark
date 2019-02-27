CFLAGS	= -DNDEBUG -march=native -Wno-char-subscripts -Wall -Ofast -std=c++14 -I. -I./sdsl-lite/build/include -fopenmp -fomit-frame-pointer -foptimize-strlen -faggressive-loop-optimizations -funroll-loops -funsafe-math-optimizations
CXXFLAGS= ${CFLAGS}
LIBS = -L./sdsl-lite/build/lib -L./sdsl-lite/build/external/libdivsufsort/lib -lz -lsdsl -ldivsufsort -ldivsufsort64 -ltbb

.PHONY: all

all: read-filter

read-filter: main.o bloomfilter.o MurmurHash3.o bloomfilter.h
	@echo "* Linking read-filter"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(LDFLAGS)

%.o: %.cpp
	@echo '* Compiling $<'
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o
