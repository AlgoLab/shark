/**
 * shark - Mapping-free filtering of useless RNA-Seq reads
 * Copyright (C) 2019 Tamara Ceccato, Luca Denti, Yuri Pirola, Marco Previtali
 *
 * This file is part of shark.
 *
 * shark is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * shark is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with shark; see the file LICENSE. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#ifndef KMER_BUILDER_HPP
#define KMER_BUILDER_HPP

#include "common.hpp"

#include "kseq.h"
#include <zlib.h>
#include <string>
#include <vector>
#include <memory>
#include "kmer_utils.hpp"

using namespace std;

class KmerBuilder {

public:
  KmerBuilder(size_t _k) : k(_k) {}

  void operator()(const vector<string>& texts, vector<vector<kmer_t>>& kmer_poss) const {
    kmer_t kmer, rckmer, key;
    for(const auto & p : texts) {
      kmer_poss.push_back(vector<kmer_t>());
      vector<kmer_t>& kmer_pos = kmer_poss.back();
      if(p.size() >= k) {
        int _pos = 0;
        kmer = build_kmer(p, _pos, k);
        if(kmer == (uint64_t)-1) continue;
        rckmer = revcompl(kmer, k);
        key = min(kmer, rckmer);
        kmer_pos.push_back(_get_hash(key));

        for (int pos = _pos; pos < (int)p.size(); ++pos) {
          uint8_t new_char = to_int[p[pos]];
          if(new_char == 0) { // Found a char different from A, C, G, T
            ++pos; // we skip this character then we build a new kmer
            kmer = build_kmer(p, pos, k);
            if(kmer == (uint64_t)-1) break;
            rckmer = revcompl(kmer, k);
            --pos; // p must point to the ending position of the kmer, it will be incremented by the for
          } else {
            --new_char; // A is 1 but it should be 0
            kmer = lsappend(kmer, new_char, k);
            rckmer = rsprepend(rckmer, reverse_char(new_char), k);
          }
          key = min(kmer, rckmer);
          kmer_pos.push_back(_get_hash(key));
        }
      }
    }
  }

private:
  const size_t k;
};

#endif
