#ifndef BF_FILLER_HPP
#define BF_FILLER_HPP

#include "kseq.h"
#include <zlib.h>
#include <string>
#include <vector>
#include <memory>
#include "bloomfilter.h"

using namespace std;

class BloomfilterFiller {
public:
  BloomfilterFiller(BF *_bf) : bf(_bf) {}

  void operator()(vector<uint64_t> *positions) const {
    if(positions) {
      for(const auto & p : *positions) {
        bf->add_at(p % bf->_size);
      }
      delete positions;
    }
  }

private:
  BF* bf;
};
#endif
