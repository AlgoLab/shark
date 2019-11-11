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
#include <array>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <string>

#include <sys/mman.h>

#include "kmer_utils.hpp"

using namespace std;
using namespace sdsl;

class KmerBuilder;
class BloomfilterFiller;

#ifndef SHARK_HUGEPAGESIZE
#define SHARK_HUGEPAGESIZE (2 * 1024 * 1024)
#endif

class BF {
  friend class KmerBuilder;
  friend class BloomfilterFiller;

public:

  typedef uint64_t kmer_t;
  typedef uint64_t hash_t;
  typedef bit_vector bit_vector_t;
  typedef bit_vector_t::rank_1_type rank_t;
  typedef vector<int> index_t;
  typedef vector<index_t> set_index_t;
  typedef int_vector<16> index_kmer_t;
  typedef bit_vector_t::select_1_type select_t;

  BF(const size_t size) :
    _size(size),
    _mode(0),
    _bf(size, 0)
  {
    char* const sptr = reinterpret_cast<char*>(_bf.data());
    const size_t soffset = SHARK_HUGEPAGESIZE - (reinterpret_cast<size_t>(sptr) % SHARK_HUGEPAGESIZE);
    char* const eptr = sptr + (((size + 63) >> 6) << 3);
    const size_t eoffset = (reinterpret_cast<size_t>(eptr) % SHARK_HUGEPAGESIZE);
    madvise(sptr + soffset, (eptr - sptr) - eoffset - soffset, MADV_HUGEPAGE);
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

  // Function to add index to a k-mer
  bool add_to_kmer(const uint64_t &kmer, const int &input_idx) {
    if (_mode != 1)
      return false;

    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      int kmer_rank = _brank(bf_idx);
      _set_index[kmer_rank].push_back(input_idx);
      return true;
    }
    return false;
  }

  // Function to add multiple indexes to a k-mer
  bool multiple_add_to_kmer(const uint64_t &kmer, const vector<int> &idxs) {
    // FIXME: can't we just use a single insert (range insertion)? Revert if it's wrong
    if (_mode != 1)
      return false;
    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      int kmer_rank = _brank(bf_idx);
      _set_index[kmer_rank].insert(idxs.begin(), idxs.begin(), idxs.end());
    }
    return true;
  }

  // Function that returns the indexes of a given k-mer
  pair<index_kmer_t::const_iterator, index_kmer_t::const_iterator> get_index(const kmer_t &kmer) const {
    int start_pos = 0;
    int end_pos = -1; // in this way, if the kmer is not in the bf, the returned IDView has no next

    #ifndef NDEBUG
    if (_mode != 2)
      return make_pair(_index_kmer.end(), _index_kmer.end());
    #endif

    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      size_t rank_searched = _brank(bf_idx + 1);
      if (rank_searched == 1) { // idxs of the first kmer
        start_pos = 0;
      } else {
        start_pos = _select_bv(rank_searched - 1) + 1;
      }
      end_pos = _select_bv(rank_searched);
      // FIXME if a k-mer has no indexes the function returns a vector with
      // only one 0
      // FIXME how to handle this situation? is the main that has to manage
      // it?
      // See also comment in switch_mode regarding dummy index
    }
    return make_pair(_index_kmer.begin()+start_pos, _index_kmer.begin()+end_pos);
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
      util::init_support(_brank,&_bf);
      size_t num_kmer = _brank(_bf.size());
      if (num_kmer != 0)
        _set_index.resize(num_kmer, vector<int>());
      return true;
    } else if(_mode == 1 and new_mode == 2) {
      _mode = new_mode;

      // We compute how many idxs we have to store
      int tot_idx = 0;
      for (auto &set : _set_index) {
        sort(set.begin(), set.end());
        set.erase(unique(set.begin(), set.end()), set.end());
        tot_idx += set.size();
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
      for (auto &set : _set_index) {
        pos += set.size();
        _bv[pos] = 1;
      }
      util::init_support(_select_bv,&_bv);

      /**
       * We merge the idxs associated to each kmer into a single
       * int_vector. This vector is the concatenation of the sets
       * associated to each kmer.
       **/
      //int_vector<16> tmp_index_kmer(tot_idx); // uncompressed and temporary
      _index_kmer = int_vector<16>(tot_idx);
      int idx_position = 0;
      for (const auto &set : _set_index) {
        // FIXME: should we check if the set is empty? Maybe saving a
        // dummy index (0)? If so, we cannot use 0 as an index for
        // kmers. Maybe -1 is better. Moreover, in this way, the main
        // must manage this (it decides the idx). Anyway, I (LD) think
        // this can never happen in our context.
        // if ( set.size() != 0) {
        for (const int &set_element : set) {
          //tmp_index_kmer[idx_position] = set_element;
          _index_kmer[idx_position] = set_element;
          ++idx_position;
        }
        // } else { tmp_index_kmer[idx_position] = 0; idx_position++; }
      }

      // _index_kmer = dac_vector<>(tmp_index_kmer);;
      /**
         On a small input (80 genes), dac compression seems the best one:
           Vec, Size, MB
           Orig, 7645238, 29.1643
           VLC, 7645238, 8.76068 //vlc_vector
           DAC, 7645238, 7.54858 // dav_vector
           ENC_eliasdelta, 7645238, 24.4415 // enc_vector<coder::elias_delta>
           ENC_eliasgamma, 7645238, 35.0523 // enc_vector<coder::elias_gamma>
           ENC_fib, 7645238, 26.5783 // enc_vector<coder::fibonacci>
           BitCompress, 1672395, 6.37969 // util::bit_compress()
       **/

      // FIXME: is this the best way to release _set_index?
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

  size_t _size;
  int _mode;
  bit_vector_t _bf;
  rank_t _brank;
  bit_vector_t _bv;
  set_index_t _set_index;
  index_kmer_t _index_kmer;
  select_t _select_bv;
};

#endif
