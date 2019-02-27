#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include "MurmurHash3.hpp"
#include <algorithm>
#include <array>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/util.hpp>
#include <string>

using namespace std;
using namespace sdsl;

class BF;
class KmerBuilder;
class BloomfilterFiller;

class IDView {
private:
  int _b, _e;
  int _p;
  BF *_bf;

public:
  IDView() : _b(), _e(), _p(), _bf(nullptr) {}
  IDView(const size_t &b, const size_t &e, BF *bf)
      : _b(b), _e(e), _p(b), _bf(bf) {}
  IDView &operator=(const IDView &rhs);
  ~IDView(){};

  bool has_next() const {return _p <= _e; }
  int_vector<16>::value_type* get_next();
  size_t size() const { return _e - _b + 1; }
  void clear();
};

class BF {
  friend class IDView;
  friend class KmerBuilder;
  friend class BloomfilterFiller;

public:
  uint64_t _get_hash(const uint64_t &kmer) const {
    const uint64_t *key = &kmer;
    array<uint64_t, 2> hashes;
    MurmurHash3_x64_128(key, sizeof(uint64_t), 0,
                        reinterpret_cast<void *>(&hashes));
    return hashes[0];
  }

public:
  BF(const size_t size) : _mode(0), _bf(size, 0) {
    _size = size;
  };
  ~BF() {}

  void add_at(const uint64_t p) {
    _bf[p] = 1;
  }

  // Function to add a k-mer to the BF
  void add_kmer(const uint64_t &kmer) {
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
  IDView get_index(const uint64_t &kmer) {
    int start_pos = 0;
    int end_pos = -1; // in this way, if the kmer is not in the bf, the returned IDView has no next

    if (_mode != 2)
      return IDView(start_pos, end_pos, nullptr);

    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      size_t rank_searched = _brank(bf_idx + 1);
      if (rank_searched == 1) { // idxs of the first kmer
        start_pos = 0;
        end_pos = _select_bv(rank_searched);
      } else {
        start_pos = _select_bv(rank_searched - 1) + 1;
        end_pos = _select_bv(rank_searched);
      }
      // FIXME if a k-mer has no indexes the function returns a vector with
      // only one 0
      // FIXME how to handle this situation? is the main that has to manage
      // it?
      // See also comment in switch_mode regarding dummy index
    }
    return IDView(start_pos, end_pos, this);
  }

  /**
   * Method to switch between modes:
   *  - 0: add k-mer to BF
   *  - 1: add indexes to kmers
   *  - 2: get indexes of k-mer
   * Note: it's not possible to go back to a previous mode
   **/
  bool switch_mode(const int &new_mode) {
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
      vector<vector<int>>().swap(_set_index);
      return true;
    } else {
      return false;
    }
  }

private:
  BF() {}
  const BF &operator=(const BF &other) { return *this; }
  const BF &operator=(const BF &&other) { return *this; }

  size_t _size;
  int _mode;
  bit_vector _bf;
  rank_support_v<1> _brank;
  bit_vector _bv;
  vector<vector<int>> _set_index;
  int_vector<16> _index_kmer;
  //dac_vector<> _index_kmer;
  select_support_mcl<1> _select_bv;
};

#endif
