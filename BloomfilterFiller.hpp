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

  void operator()(vector<uint64_t> *positions) {
    {
      std::lock_guard<std::mutex> lock(mtx);
      for(const auto p : *positions) {
        bf->add_at(p % bf->_size);
      }
    }
    delete positions;
  }

private:
  BF *const bf;
  std::mutex mtx;

};
#endif
