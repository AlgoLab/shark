#ifndef KMER_BUILDER_HPP
#define KMER_BUILDER_HPP

#include "kseq.h"
#include <zlib.h>
#include <string>
#include <vector>
#include <memory>
#include "bloomfilter.h"
#include "kmer_utils.hpp"

using namespace std;

class KmerBuilder {

public:
  KmerBuilder(size_t _k) : k(_k) {}

  vector<uint64_t>* operator()(vector<pair<string, string>> *texts) const {
    if(texts) {
      vector<uint64_t>* kmer_pos = new vector<uint64_t>();
      array<uint64_t, 2> hashes;
      uint64_t kmer, rckmer, key;
      for(const auto & p : *texts) {
        if(p.second.size() >= k) {
          int _pos = 0;
          kmer = build_kmer(p.second, &_pos, k);
          if(kmer == (uint64_t)-1) continue;
          rckmer = revcompl(kmer, k);
          key = min(kmer, rckmer);
          MurmurHash3_x64_128(&key, sizeof(uint64_t), 0,
                              reinterpret_cast<void *>(&hashes));
          kmer_pos->push_back(hashes[0]);

          for (int pos = _pos; pos < (int)p.second.size(); ++pos) {
            uint8_t new_char = to_int[p.second[pos]];
            if(new_char == 0) { // Found a char different from A, C, G, T
              ++pos; // we skip this character then we build a new kmer
              kmer = build_kmer(p.second, &pos, k);
              if(kmer == (uint64_t)-1) break;
              rckmer = revcompl(kmer, k);
              --pos; // p must point to the ending position of the kmer, it will be incremented by the for
            } else {
              --new_char; // A is 1 but it should be 0
              kmer = lsappend(kmer, new_char, k);
              rckmer = rsprepend(rckmer, reverse_char(new_char), k);
            }
            key = min(kmer, rckmer);
            MurmurHash3_x64_128(&key, sizeof(uint64_t), 0,
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
