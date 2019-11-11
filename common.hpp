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

#ifndef SHARK_COMMON_HPP
#define SHARK_COMMON_HPP

#include <string>
#include <utility>

using std::string;

struct sharseq_t {
  string id, seq, qual;
};

typedef std::pair<string, std::pair<sharseq_t, sharseq_t>> elem_t;

typedef std::pair<string, std::pair<sharseq_t, sharseq_t>> assoc_t;

#endif
