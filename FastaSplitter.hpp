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


template <size_t N>
class FastaSplitter {
public:
  FastaSplitter(kseq_t * const _seq, std::vector<std::string>* const _ids = nullptr)
    : seq(_seq), ids(_ids), idx(0)
  { }

  ~FastaSplitter() {
  }

  void operator()(vector<string>& fasta, size_t& base_idx) {
    std::lock_guard<std::mutex> lock(mtx);
    int seq_len;
    base_idx = idx;
    while (fasta.size() < N && (seq_len = kseq_read(seq)) >= 0) {
      if (ids != nullptr) ids->push_back(seq->name.s);
      fasta.emplace_back(seq->seq.s);
      ++idx;
    }
  }

private:
  kseq_t * const seq;
  std::vector<std::string>* const ids;
  std::mutex mtx;
  size_t idx;

};

#endif
