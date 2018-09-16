#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include <algorithm>
#include <array>
#include <sdsl/bit_vectors.hpp>
#include <string>

#include "MurmurHash3.hpp"
#include "KMC/kmc_api/kmc_file.h"

using namespace std;
using namespace sdsl;
typedef bit_vector::size_type size_type;



static const char RCN[128] = {
        0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   //  0
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
    BF(const size_t size) :  _bf(size, 0) { _size = size; } vector<vector<int> > set_index;  
    ~BF() {}



    void add_kmer(const string &kmer) {
        uint64_t hash = _get_hash(kmer);
        _bf[hash % _size] = 1;
    }

    bool test_kmer(const string &kmer) const {
        uint64_t hash = _get_hash(kmer);
        return _bf[hash % _size];
    }
    
    bool get_rank_kmer(const string &kmer) const {
    	if(!_mode)
    		return false;
        uint64_t hash = _get_hash(kmer);
        size_t bf_idx = hash % _size;
        if(_bf[bf_idx])
        	return _brank(bf_idx);
    }
    
    int get_idx (const string &kmer){
        uint64_t hash = _get_hash(kmer);
        size_t bf_idx = hash % _size;
    	return bf_idx;
    }
    
   void switch_mode() {
	_mode = true;
	_brank = rank_support_v<1>(&_bf);
	bv = bit_vector(100,0);

    }
  
    bool add_to_kmer(const string &kmer, ifstream& input) {

    	if (!_mode)
            return false;

        uint64_t hash = _get_hash(kmer);
	int pos;
        size_t bf_idx = hash % _size;
	int i = 0;
	int key = bf_idx;
	string line;
	int count = 0;
	int newrank = 0;

	if (_bf[bf_idx]) {
	        size_t cnts_idx = _brank(bf_idx);
	        int num_idx=0;
		if (input.is_open()) {
	        	while (getline(input, line)) {
		       		int a = std::stoi(line);
	               		index.push_back(a);
				num_idx++;
				}
		}
		input.close();
		_rank = rank_support_v<1>(&bv);
		newrank=_rank(bv.size());

		if (newrank != 0){
			for(int i=0; count != newrank && i < bv.size() ;i++ ){
				pos = i;
				if(bv[i] == 1)
					count++;
			}
		}else
			pos = 0;

		pos = pos+num_idx;
	
		if(bv[pos]!=1)
			bv[pos] = 1;

        	return true;
	}
    }

    vector<int> get_index(const string &kmer){
    	uint64_t hash = _get_hash(kmer);
        size_t bf_idx = hash % _size;
        size_t rank_searched = _brank(bf_idx+1);
        vector<int> index_res;
 	int count =0;
	int pos;
	int i;
	int end_pos;
	
        if (_bf[bf_idx]){
		if(rank_searched ==1){
			pos = 0;
			for(i=0; bv[i] ==0 ;i++ );
			end_pos =i;

		}else{
			for(i=0; count != rank_searched-1 && i < bv.size() ;i++ ){
				pos = i;
				if(bv[i] == 1)
					count++;
			}
			while(bv[i]!=1)
				i++;
		        end_pos=i;
		}
		for(int i=pos; i <end_pos; i++)
			index_res.push_back(index[i]);

                return index_res;
        }
    }  

private:
    BF() {}
    const BF &operator=(const BF &other) { return *this; }
    const BF &operator=(const BF &&other) { return *this; }

    bool _mode; // false = write, true = read
    size_t _size;
    bit_vector _bf;
    rank_support_v<1> _brank;
    rank_support_v<1> _rank;
    int_vector<8> _counts;
    vector<int> index;
    bit_vector bv;


};

#endif
