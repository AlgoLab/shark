#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include "KMC/kmc_api/kmc_file.h"
#include "MurmurHash3.hpp"
#include <algorithm>
#include <array>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/select_support.hpp>
#include <string>

using namespace std;
using namespace sdsl;

static const char RCN[128] = {
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 0
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 10
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 20
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 30
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 40
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 50
    0,   0,   0, 0,   0,   'T', 0,   'G', 0,   0,   // 60
    0,   'C', 0, 0,   0,   0,   0,   0,   'N', 0,   // 70
    0,   0,   0, 0,   'A', 0,   0,   0,   0,   0,   // 80
    0,   0,   0, 0,   0,   0,   0,   'T', 0,   'G', // 90
    0,   0,   0, 'G', 0,   0,   0,   0,   0,   0,   // 100
    'N', 0,   0, 0,   0,   0,   'A', 0,   0,   0,   // 110
    0,   0,   0, 0,   0,   0,   0,   0              // 120
};

class BF {

private:
  static const char _rcc(const char &c) { return RCN[c]; }

  static const string _rc(const string &kmer) {
    string rc(kmer);
    transform(rc.begin(), rc.end(), rc.begin(), _rcc);
    reverse(rc.begin(), rc.end());
    return rc;
  }

  const string _minrc(const string &kmer) const { return min(kmer, _rc(kmer)); }

  uint64_t _get_hash(const string &kmer) const {
    string k = _minrc(kmer);
    array<uint64_t, 2> hashes;
    MurmurHash3_x64_128(k.c_str(), k.size(), 0,
                        reinterpret_cast<void *>(&hashes));
    return hashes[0];
  }

public:
  BF(const size_t size) : _mode(0), _check_mode(false), _bf(size, 0) {
    _size = size;
  };
  ~BF() {}

  // function to add a k-mer to BF and initialize vectors
  void add_kmer(const string &kmer) {
    //_check_mode = false --> _mode = 0, user can add k-mers to BF
    if (!_check_mode) {
      uint64_t hash = _get_hash(kmer);
      _bf[hash % _size] = 1;
    }
  }

  // function to test the presence of a k-mer in the BF
  bool test_kmer(const string &kmer) const {
    uint64_t hash = _get_hash(kmer);
    return _bf[hash % _size];
  }

  // switch between modes: 0 = add k-mer to BF, 1 = add indexes of k-mer, 2 =
  // get indexes of k-mer
  // it's not possible to go back to a previous mode
  bool switch_mode(int user_input) {
    int tot_idx = 0;
    _mode = user_input;

    // if check_mode is false the only mode used is _mode = 0
    if (!_check_mode) {
      if (_mode == 1) {
        size_t num_kmer;
        _check_mode = true;
        _brank = rank_support_v<1>(&_bf);
        // FIXME: +1 ? What happens if there are 0 kmers in the bf?
        // we need to add 1 because _set_index[0] will always be empty, because
        // index 0 doesn't match any rank value
        num_kmer = _brank.rank(_bf.size()) + 1;
        // test

        if (num_kmer != 1)
          _set_index.resize(num_kmer, vector<int>());
        return true;
      } else
        return false;

    } else if (_check_mode) {

      // if check_mode is true _mode = 1 has been used

      if (_mode == 2) {

        // sorting indexes for each k-mer
        for (auto set : _set_index)
          if (set.size())
            sort(begin(set), end(set), less<int>());

        // remove duplicates
        for (int i = 0; i < _set_index.size(); i++)
          _set_index[i].erase(
              unique(_set_index[i].begin(), _set_index[i].end()),
              _set_index[i].end());

        for (auto set : _set_index)
          tot_idx += set.size();

        _bv = bit_vector(tot_idx, 0);

        // storage indexes in the int_vector
        _index_kmer = int_vector<64>(tot_idx);
        int idx_position = 0;

        // i starts from 1 because the firse element in _set_index will always
        // be an empty vector, because no k-mer has rank = 0
        for (int i = 1; i < _set_index.size(); i++) {

          // _set_index[i].size != 0 because a vector in _set_index with size 0
          // is an empty vector, and there are no indexes to add to _index_kmer
          if (_set_index[i].size() != 0) {
            for (int &set_element : _set_index[i]) {
              _index_kmer[idx_position] = set_element;
              idx_position++;
            }
          } else {
            // when a vector is empty we assign a 0 in _index_kmer (there is a
            // k-mer in BF that has no indexes)
            _index_kmer[idx_position] = 0;
            idx_position++;
          }
        }

        // FIXME: compression of index k-mer
        // util::bit_compress(_index_kmer);

        int pos = -1;
        // set _bv elements to 1 at the end of each index range
        for (int i = 1; i < _set_index.size(); i++) {
          pos += _set_index[i].size();
          _bv[pos] = 1;
        }
        _select_bv = select_support_mcl<1>(&_bv);

        // TODO: release _set_index here
        // is this the best way to do this?
        _set_index.clear();
        return true;
      } else
        return false;
    } else
      return false;
  }

  // add index of a given k-mer

  bool add_to_kmer(const string &kmer, int input_idx) {
    if (_mode != 1)
      return false;

    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;

    if (_bf[bf_idx]) {
      // bf_idx + 1 because _brank start from position 0
      int kmer_rank = _brank(bf_idx + 1);
      // storage in _set_index, the index, in position equals to k-mer's rank on
      // bloom filter
      _set_index[kmer_rank].push_back(input_idx);

      return true;
    }
    return false;
  }

  // add a vector of indexes of a k-mer

  bool multiple_add_to_kmer(const string &kmer, vector<int> idx_vector) {
    for (int i = 0; i < idx_vector.size(); i++)
      if (!add_to_kmer(kmer, idx_vector[i]))
        return false;
    return true;
  }

  // function that returns the indexes of a given k-mer
  // FIXME: to improve this function you could avoid copying back a vector<int>
  //        and provide an iterator over _index_kmer instead.

  vector<int> get_index(const string &kmer) {
    if (_mode != 2)
      return {};

    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    size_t rank_searched = _brank(bf_idx + 1);
    vector<int> index_res;
    int start_pos = 0;
    int end_pos = 0;
    int bv_idx;

    // select on _bv
    if (_bf[bf_idx]) {

      // for a single value _select_bv must be done in one step, because with
      // the steps below it fails for _select_bv(0)
      if (rank_searched == 1) {
        bv_idx = _select_bv(1);
        index_res.push_back(_index_kmer[bv_idx]);
        return index_res;
      } else {
        start_pos = _select_bv(rank_searched - 1) + 1;
        end_pos = _select_bv(rank_searched);
      }

      // output vector
      // FIXME if a k-mer has no indexes the function returns a vector with only
      // one 0
      // FIXME how to handle this situation? is the main that has to manage it?
      for (int i = start_pos; i <= end_pos; i++) {
        index_res.push_back(_index_kmer[i]);
      }

      return index_res;
    }
    return {};
    // FIXME: return an empty vector? (return vector<int>())
    // throw an exception?
  }

private:
  BF() {}
  const BF &operator=(const BF &other) { return *this; }
  const BF &operator=(const BF &&other) { return *this; }

  size_t _size;
  int _mode;
  bool _check_mode;
  bit_vector _bf;
  rank_support_v<1> _brank;
  bit_vector _bv;
  vector<vector<int>> _set_index;
  int_vector<64> _index_kmer;
  select_support_mcl<1> _select_bv;
};

#endif
