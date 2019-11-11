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

#ifndef FASTQ_SPLITTER_HPP
#define FASTQ_SPLITTER_HPP

#include "common.hpp"
#include "kseq.h"
#include <string>
#include <vector>
#include <memory>

#include <tbb/pipeline.h>

using namespace std;

class FastqSplitter {
public:

  typedef vector<elem_t> output_t;

  FastqSplitter(kseq_t * const _seq1, kseq_t * const _seq2, const int _maxnum, const char _min_quality, const bool _full_mode)
    : seq1(_seq1), seq2(_seq2), maxnum(_maxnum), min_quality(_min_quality), full_mode(_full_mode), empty_el({})
  {
  }

  ~FastqSplitter() {
  }

  output_t* operator()(tbb::flow_control &fc) const {
    output_t* const fastq = new output_t();
    fastq->reserve(maxnum);
    int seq_len1, seq_len2;
    if (min_quality == 0) {
      if (seq2 == nullptr) {
        while (fastq->size() < maxnum && (seq_len1 = kseq_read(seq1)) >= 0) {
          fastq->push_back({
            seq1->seq.s,
            { { seq1->name.s, full_mode ? seq1->seq.s : "", full_mode ? seq1->qual.s : "" },
              empty_el }
          });
        }
      } else {
        while (fastq->size() < maxnum && (seq_len1 = kseq_read(seq1)) >= 0 && (seq_len2 = kseq_read(seq2)) >= 0) {
          fastq->push_back({
            string(seq1->seq.s) + "N" + string(seq2->seq.s),
            { { seq1->name.s, full_mode ? seq1->seq.s : "", full_mode ? seq1->qual.s : "" },
              { seq2->name.s, full_mode ? seq2->seq.s : "", full_mode ? seq2->qual.s : "" } }
          });
        }
      }
    } else {
      const char mq = min_quality + 33;
      if (seq2 == nullptr) {
        while (fastq->size() < maxnum && (seq_len1 = kseq_read(seq1)) >= 0) {
          fastq->push_back({
            mask_seq(seq1->seq.s, seq1->qual.s, seq1->qual.l, mq),
            { { seq1->name.s, full_mode ? seq1->seq.s : "", full_mode ? seq1->qual.s : "" },
              empty_el }
          });
        }
      } else {
        while (fastq->size() < maxnum && (seq_len1 = kseq_read(seq1)) >= 0 && (seq_len2 = kseq_read(seq2)) >= 0) {
          fastq->push_back({
            mask_seq(
              string(seq1->seq.s) + "N" + string(seq2->seq.s),
              string(seq1->qual.s) + "\33" + string(seq2->qual.s),
              mq
            ),
            { { seq1->name.s, full_mode ? seq1->seq.s : "", full_mode ? seq1->qual.s : "" },
              { seq2->name.s, full_mode ? seq2->seq.s : "", full_mode ? seq2->qual.s : "" } }
          });
        }
      }
    }
    if(fastq->size() > 0) return fastq;
    fc.stop();
    delete fastq;
    return nullptr;
  }
private:
  kseq_t * const seq1;
  kseq_t * const seq2;
  const size_t maxnum;
  const char min_quality;
  const bool full_mode;
  const sharseq_t empty_el;

  static string mask_seq(string seq, const char* const qual, const size_t l, const char min_quality) {
    for (size_t i = 0; i < l; ++i) {
      if (qual[i] < min_quality) seq[i] = seq[i] - 64;
    }
    return seq;
  }

  static string mask_seq(string seq, const string& qual, const char min_quality) {
    return mask_seq(seq, qual.c_str(), qual.length(), min_quality);
  }
};

#endif
