#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include <algorithm>
#include <array>
#include <sdsl/bit_vectors.hpp>
#include <string>
#include <sdsl/int_vector.hpp>

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
    BF(const size_t size) :  _bf(size, 0) { _size = size; }  ; 
    ~BF() {}



    void add_kmer(const string &kmer) {
        uint64_t hash = _get_hash(kmer);
        _bf[hash % _size] = 1;
	bv = bit_vector(100,0);
    }

    bool test_kmer(const string &kmer) const {
        uint64_t hash = _get_hash(kmer);
        return _bf[hash % _size];
    }
    
    //to be removed
    bool get_rank_kmer(const string &kmer) const {
    	
        uint64_t hash = _get_hash(kmer);
        size_t bf_idx = hash % _size;
        if(_bf[bf_idx])
        	return _brank(bf_idx);
    }
    
    //to be removed
    int get_idx (const string &kmer){
        uint64_t hash = _get_hash(kmer);
        size_t bf_idx = hash % _size;
    	return bf_idx;
    }
    
    //function to add indexes of a known kmer
    bool add_to_kmer(const string &kmer, ifstream& input) {
	
	_brank = rank_support_v<1>(&_bf);
    	uint64_t hash = _get_hash(kmer);
	int pos;
        size_t bf_idx = hash % _size;
	int i = 0;
	int key = bf_idx;
	string line;
	int count = 0;
	int newrank = 0;
	int num_idx;


	if (_bf[bf_idx]) {
	        size_t cnts_idx = _brank(bf_idx);
		if(set_index[key].size() == 0)
		        int num_idx=0;
		else
			int num_idx = set_index[key].size();
		if (input.is_open()) {
	        	while (getline(input, line)) {
		       		int a = std::stoi(line);
				set_index[key][num_idx];	               		
				//index.push_back(a);
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

		for(int j = 0; j<set_index.size()-1;j++){
			set_index[j].insert( set_index[j].end(), set_index[j+1].begin(), set_index[j+1].end() );
		}

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
			//index_res.push_back(index[i]);

                return index_res;
        }
    }  

private:
    BF() {}
    const BF &operator=(const BF &other) { return *this; }
    const BF &operator=(const BF &&other) { return *this; }


    size_t _size;
    uint8_t size_set = 10000;
    bit_vector _bf;
    rank_support_v<1> _brank;
    rank_support_v<1> _rank;
    vector<int> index;
    bit_vector bv;
    vector< vector<int> > set_index;
    int_vector<> index_kmer();

};

#endif
