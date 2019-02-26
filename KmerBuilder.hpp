#ifndef KMER_BUILDER_HPP
#define KMER_BUILDER_HPP

#include "kseq.h"
#include <zlib.h>
#include <string>
#include <vector>
#include <memory>
#include "bloomfilter.h"

using namespace std;

class KmerBuilder {

public:
  KmerBuilder(size_t _k) : k(_k) {}

  vector<uint64_t>* operator()(vector<pair<string, string>> *texts) const {
    if(texts) {
      vector<uint64_t>* kmer_pos = new vector<uint64_t>();
      array<uint64_t, 2> hashes;
      for(const auto & p : *texts) {
        if(p.second.size() >= k) {
          for(size_t i = 0; i < p.second.size() - k + 1; ++i) {
            string kmer = BF::_minrc(p.second.substr(i, k));
            MurmurHash3_x64_128(kmer.c_str(), kmer.size(), 0,
                                reinterpret_cast<void *>(&hashes));
            kmer_pos->push_back(hashes[0]);
          }
        }
      }
      delete texts;
      return kmer_pos;
    }
    return NULL;
  }

private:
  size_t k;
};

#endif
