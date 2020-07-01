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
#include <string>
#include <vector>
#include <memory>
#include <mutex>

using namespace std;

class FastaSplitter {
public:
  FastaSplitter(kseq_t * const _seq, const int _maxnum, vector<string>* const _ids = nullptr)
    : seq(_seq), maxnum(_maxnum), ids(_ids)
  { }

  ~FastaSplitter() {
  }

  vector<pair<string, string>>* operator()() {
    std::lock_guard<std::mutex> lock(mtx);
    vector<pair<string, string>>* const fasta = new vector<pair<string, string>>();
    fasta->reserve(maxnum);
    int seq_len;
    while(fasta->size() < maxnum && (seq_len = kseq_read(seq)) >= 0) {
      if (ids != nullptr) ids->push_back(seq->name.s);
      fasta->emplace_back(seq->name.s, seq->seq.s);
    }
    if (!fasta->empty()) return fasta;
    delete fasta;
    return nullptr;
  }

private:
  kseq_t * const seq;
  const size_t maxnum;
  vector<string>* const ids;
  std::mutex mtx;

};

#endif
