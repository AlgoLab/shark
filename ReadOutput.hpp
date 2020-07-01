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

#ifndef READOUTPUT__HPP
#define READOUTPUT__HPP

#include <iostream>
#include <vector>
#include <mutex>

#include "common.hpp"

class ReadOutput {
public:
  ReadOutput(FILE* const _out1 = nullptr, FILE* const _out2 = nullptr)
    : out1(_out1), out2(_out2)
  { }

  void operator()(const std::vector<assoc_t>& associations) {
    std::lock_guard<std::mutex> lock(mtx);
    string previd = "";
    for(const auto & a : associations) {
      const sharseq_t& s1 = a.second.first;
      const sharseq_t& s2 = a.second.second;
      printf("%s %s\n", s1.id.c_str(), a.first.c_str());
      if (out1 != nullptr && previd != s1.id)
        fprintf(out1, "@%s\n%s\n+\n%s\n", s1.id.c_str(), s1.seq.c_str(), s1.qual.c_str());
      if (out2 != nullptr && previd != s1.id)
        fprintf(out2, "@%s\n%s\n+\n%s\n", s2.id.c_str(), s2.seq.c_str(), s2.qual.c_str());
      previd = std::move(s1.id);
    }
  }

private:
  FILE* const out1;
  FILE* const out2;
  std::mutex mtx;
};

#endif
