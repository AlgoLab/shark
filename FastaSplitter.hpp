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

#ifndef FASTA_SPLITTER_HPP
#define FASTA_SPLITTER_HPP

#include "kseq.h"
#include <zlib.h>
#include <string>
#include <vector>
#include <memory>

#include <tbb/pipeline.h>

using namespace std;

class FastaSplitter {
public:
  FastaSplitter(kseq_t * const _seq, const int _maxnum)
    : seq(_seq), maxnum(_maxnum)
  { }

  ~FastaSplitter() {
  }

  vector<pair<string, string>>* operator()(tbb::flow_control &fc) const {
    vector<pair<string, string>>* const fasta = new vector<pair<string, string>>();
    fasta->reserve(maxnum);
    int seq_len;
    while(fasta->size() < maxnum && (seq_len = kseq_read(seq)) >= 0) {
      fasta->emplace_back(seq->name.s, seq->seq.s);
    }
    if(fasta->size() > 0) return fasta;
    fc.stop();
    delete fasta;
    return NULL;
  }

private:
  kseq_t * const seq;
  const size_t maxnum;

};

#endif
