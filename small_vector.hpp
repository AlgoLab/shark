/**
 * shark - Mapping-free filtering of useless RNA-Seq reads
 * Copyright (C) 2020 Tamara Ceccato, Luca Denti, Yuri Pirola, Marco Previtali
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

#include <cstdint>
#include <vector>

struct small_vector_t {
  union {
    struct {
      uint8_t flag;
      uint8_t size;
      uint16_t arr[3];
    } s;
    std::vector<uint16_t>* l;
  } v;

  small_vector_t() {
    v.s.flag = 1u;
    v.s.size = 0;
  }

  ~small_vector_t() {
    if ((v.s.flag & 0x1) == 0) {
      delete v.l;
    }
  }

  void push_back(uint16_t x) {
    if ((v.s.flag & 0x1) != 0) {
      if (v.s.size < 3) v.s.arr[v.s.size++] = x;
      else {
        std::vector<uint16_t>* ptr = new std::vector<uint16_t>(v.s.arr, v.s.arr + 3);
        ptr->push_back(x);
        v.l = ptr;
      }
    } else {
      v.l->push_back(x);
    }
  }

  size_t size() const {
    if ((v.s.flag & 0x1) != 0) {
      return v.s.size;
    } else {
      return v.l->size();
    }
  }

  uint16_t last() const {
    if ((v.s.flag & 0x1) != 0) {
      return v.s.arr[v.s.size - 1];
    } else {
      return v.l->back();
    }
  }

  const uint16_t* begin() const {
    if ((v.s.flag & 0x1) != 0) {
      return v.s.arr;
    } else {
      return v.l->data();
    }
  }

  const uint16_t* end() const {
    if ((v.s.flag & 0x1) != 0) {
      return v.s.arr + v.s.size;
    } else {
      return v.l->data() + v.l->size();
    }
  }

};
