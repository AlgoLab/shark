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
  BF(const size_t size) : _mode(0), check_mode(false), _bf(size, 0) { _size = size; };
  ~BF() {}

  // function to add a k-mer to BF and initialize vectors
  void add_kmer(const string &kmer) {
    uint64_t hash = _get_hash(kmer);
    _bf[hash % _size] = 1;
  }

  // function to test the presence of a k-mer in the BF
  bool test_kmer(const string &kmer) const {
    uint64_t hash = _get_hash(kmer);
    return _bf[hash % _size];
  }

  // switch between modes: 0 = add k-mer to BF, 1 = add indexes of k-mer, 2 =
  // get indexes of k-mer
  void switch_mode(int user_input /* FIXME: i? */) {
    // FIXME: It's possible to switch back from _mode=1 to _mode=0
    int tot_idx = 0;

    int pos = -1; // FIXME: not used if i=1
    _mode = user_input;

    //if check_mode is false the only mode used is _mode = 0
    if (!check_mode){
	    if (_mode == 1) {
	      check_mode = true;
	      _brank = rank_support_v<1>(&_bf);
	      // FIXME: +1 ? What happens if there are 0 kmers in the bf?
	      num_kmer = _brank.rank(_bf.size()) + 1;
	      set_index.resize(num_kmer, vector<int>());
	    }
    }
    
    //if check_mode is true _mode = 1 has been used
    if (check_mode){
	    if (_mode == 2) {
	      
	      for(auto set : set_index) 
	      	tot_idx += set.size();
	      bv = bit_vector(tot_idx, 0);

	      // sorting indexes for each k-mer
	      for(auto set : set_index) 
		if (set.size())
		  std::sort(begin(set), end(set), std::less<int>());

	      // storage indexes in the int_vector
	      // FIXME: i is already defined in this function (it's a parameter, change name)
	      //        We can also use a c++11 range-for here
	      index_kmer = int_vector<64>(tot_idx);
	      int tmp = 0; // FIXME: ugly variable name
	      for (int i = 0; i < set_index.size(); i++) {
		// (?) BUG: why set_index[i][j] != 0?
		for (int j = 0; j < set_index[i].size() && set_index[i][j] != 0; j++) {
		  index_kmer[tmp] = set_index[i][j];
		  tmp++;
		}
	      }

	      // compression of index k-mer
	      util::bit_compress(index_kmer);

	      int pos = -1;
	      // set bv elements to 1 at the end of each index range
	      for(auto set : set_index){
		pos += set.size();
		bv[pos] = 1;
	      }
	      // TODO: shouldn't you build a rank support data structure over bv here?
	      //       what about the select support data structure?
	    }
	    // TODO: release set_index here
    }
  }

  // add index of a given k-mer

  bool add_to_kmer(const string &kmer, int input_idx) {
    if (_mode != 1)
      return false;

    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;

    if (_bf[bf_idx]) {
      int kmer_rank = _brank(bf_idx + 1);
      // storage in set_index, the index, in position equals to k-mer's rank on
      // bloom filter
      set_index[kmer_rank].push_back(input_idx);
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
  //        and provide an iterator over index_kmer instead.
  vector<int> get_index(const string &kmer) {
    if (_mode != 2)
      return vector<int>();

    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    size_t rank_searched = _brank(bf_idx + 1);
    vector<int> index_res;
    int start_pos = 0;
    int end_pos = 0;
    int idx_bv; // FIXME: we use bf_idx and idx_bv.  Be consistent

    // select on bv
    if (_bf[bf_idx]) {
      // FIXME: We build the select support data structure each time we perform a query!
      //        Highly inefficient
      _select_bv = select_support_mcl<1>(&bv);
      // TODO: describe why you check if rank_searched == 1
      if (rank_searched == 1) {
        idx_bv = _select_bv(1);
      } else {
        start_pos = _select_bv(rank_searched - 1) + 1;
        end_pos = _select_bv(rank_searched);
      }
      // output vector
      for (int i = start_pos; i <= end_pos; i++)
        index_res.push_back(index_kmer[i]);
      return index_res;
    }
    // FIXME: return missing?
  }

private:
  BF() {}
  const BF &operator=(const BF &other) { return *this; }
  const BF &operator=(const BF &&other) { return *this; }

  // TODO: private attributes names should be consistent (either all start with _ or none does)
  size_t num_kmer; // TODO: this is a local variable not an attribute of the class
  size_t _size;
  int _mode;
  bool check_mode;
  bit_vector _bf;
  rank_support_v<1> _brank;
  rank_support_v<1> _rank; // TODO: unused attribute (?)
  vector<int> index_for_bv; // TODO: unused attribute
  bit_vector bv;
  vector<vector<int>> set_index;
  int_vector<64> index_kmer;
  select_support_mcl<1> _select_bv;
};

#endif
