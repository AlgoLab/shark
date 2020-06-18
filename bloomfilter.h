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

#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include <algorithm>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <string>

#include "kmer_utils.hpp"

using namespace std;
using namespace sdsl;

class KmerBuilder;
class BloomfilterFiller;

class BF {
  friend class KmerBuilder;
  friend class BloomfilterFiller;

public:

  typedef uint64_t kmer_t;
  typedef uint64_t hash_t;
  typedef bit_vector bit_vector_t;
  typedef bit_vector_t::rank_1_type rank_t;
  typedef vector<int> index_t;
  typedef vector<uint16_t> set_index_t;
  typedef vector<uint16_t> index_kmer_t;
  typedef bit_vector_t::select_1_type select_t;

  BF(const size_t size) :
    _size(size),
    _mode(0),
    _bf(size, 0)
  {
  }

  ~BF() {}

  void add_at(const uint64_t p) {
    _bf[p] = 1;
  }

  // Function to add a k-mer to the BF
  void add_kmer(const kmer_t kmer) {
    if (_mode == 0) {
      uint64_t hash = _get_hash(kmer);
      _bf[hash % _size] = 1;
    }
  }

  // Function to test if a k-mer is in the BF
  bool test_kmer(const uint64_t &kmer) const {
    uint64_t hash = _get_hash(kmer);
    return _bf[hash % _size];
  }

  void add_to_kmer_1(vector<uint64_t> &kmers) {
    if (_mode != 1)
      return;

    uint64_t prev = _size;
    for (const auto bf_idx: kmers) {
      if (bf_idx == prev) continue;
      prev = bf_idx;
      const int kmer_rank = _brank(bf_idx);
      ++_set_index[kmer_rank];
    }
  }

  void add_to_kmer_2(vector<uint64_t> &kmers, const uint16_t input_idx) {
    if (_mode != 2)
      return;

    uint64_t prev = _size;
    for (const auto bf_idx: kmers) {
      if (bf_idx == prev) continue;
      prev = bf_idx;
      int kmer_rank = _brank(bf_idx);
      int start_pos = (kmer_rank >= 1) ? _select_bv(kmer_rank) + 1 : 0;
      _index_kmer[start_pos + --_set_index[kmer_rank]] = input_idx;
    }
  }

  // Function that returns the indexes of a given k-mer
  pair<index_kmer_t::const_iterator, index_kmer_t::const_iterator> get_index(const kmer_t &kmer) const {
    int start_pos = 0;
    int end_pos = -1; // in this way, if the kmer is not in the bf, the returned IDView has no next

    #ifndef NDEBUG
    if (_mode != 3)
      return make_pair(_index_kmer.end(), _index_kmer.end());
    #endif

    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      size_t rank_searched = _brank(bf_idx + 1);
      if (rank_searched > 1) { // idxs of the first kmer
        start_pos = _select_bv(rank_searched - 1) + 1;
      }
      end_pos = _select_bv(rank_searched);
    }
    return make_pair(_index_kmer.cbegin()+start_pos, _index_kmer.cbegin()+end_pos);
  }

  /**
   * Method to switch between modes:
   *  - 0: add k-mer to BF
   *  - 1: add indexes to kmers
   *  - 2: get indexes of k-mer
   * Note: it's not possible to go back to a previous mode
   **/
  bool switch_mode(const int new_mode) {
    if(_mode == 0 and new_mode == 1) {
      /**
       * Here we initialize the vector that will contain, for each
       * kmer, the set of idx associated to it. The idxs of the i-th
       * kmer in the Bloom filter, will be at the i-th position in
       * this vector. This vector is filled with the 'add_to_kmer'
       * function.
       **/
      _mode = new_mode;
      util::init_support(_brank, &_bf);
      size_t num_kmer = _brank(_bf.size());
      _set_index.resize(num_kmer, 0);

      return true;
    } else if(_mode == 1 and new_mode == 2) {
      _mode = new_mode;

      // We compute how many idxs we have to store
      int tot_idx = 0;
      for (const int set : _set_index) {
        tot_idx += set;
      }

      /**
       * We build a bit vector that stores the "sizes" of the sets
       * associated to each kmer. The bv has a 1 at the end of each
       * set.  Example: [{1,2}, {1,3,4}, {2}] -> 010011 Underlying
       * data structure for the int_vector that store the
       * concatenation of all the sets.
       **/
      _bv = bit_vector(tot_idx, 0);
      int pos = -1;
      for (const int set: _set_index) {
        pos += set;
        _bv[pos] = 1;
      }
      util::init_support(_select_bv, &_bv);

      /**
       * We merge the idxs associated to each kmer into a single
       * int_vector. This vector is the concatenation of the sets
       * associated to each kmer.
       **/
      _index_kmer.resize(tot_idx);

      return true;
    } else if(_mode == 2 and new_mode == 3) {
      _mode = new_mode;
      set_index_t().swap(_set_index);
      return true;
    } else {
      return false;
    }
  }

private:
  BF() = delete;
  const BF &operator=(const BF &) = delete;
  const BF &operator=(const BF &&) = delete;

  const size_t _size;
  int _mode;
  bit_vector_t _bf;
  rank_t _brank;
  bit_vector_t _bv;
  set_index_t _set_index;
  index_kmer_t _index_kmer;
  select_t _select_bv;
};

#endif
