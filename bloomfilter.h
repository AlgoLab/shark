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
    BF(const size_t size) :  _bf(size, 0) { _size = size; } ; 
    ~BF() {}


    //function to add a k-mer to BF and initialize vectors
    
    void add_kmer(const string &kmer) {
        uint64_t hash = _get_hash(kmer);
        _bf[hash % _size] = 1;
	bv = bit_vector(100,0);
	set_index.resize(100, vector<int>(50));
	index_kmer = int_vector<64> (100);
	index_for_bv = vector<int>(100);
	
    }

    //function to test the presence of a k-mer in the BF
    
    bool test_kmer(const string &kmer) const {
        uint64_t hash = _get_hash(kmer);
        return _bf[hash % _size];
    }
      
    
    //function to add indexes of a given kmer
    
    bool add_to_kmer(const string &kmer, ifstream& input) {
	
	_brank = rank_support_v<1>(&_bf);
    	uint64_t hash = _get_hash(kmer);
        size_t bf_idx = hash % _size;
	int i = 0;
	string line;
	int num_idx;
	
	if (_bf[bf_idx]) {
		
	        int kmer_rank = _brank(bf_idx+1);
		//num_idx is the position where indexes start to be stored
		if(set_index[kmer_rank][0] == 0){
		       num_idx=0;
		}else{
			while(set_index[kmer_rank][i] != 0)
				i++;
		}		
		num_idx = i;
		//reading of the file with the indexes and storage in set_index, in a position according to the k-mer rank
		if (input.is_open()) {
	        	while (getline(input, line)) {
		       		int a = std::stoi(line);
				set_index[kmer_rank][num_idx]=a;	               		
				num_idx++;
				}
		}
		input.close();

		//vector that contains the number of indexes of each k-mer
		index_for_bv[kmer_rank] = num_idx;

        	return true;
	}
    }

    
    

    //function that returns the indexes of a given k-mer
    
    vector<int> get_index(const string &kmer){
    	
	uint64_t hash = _get_hash(kmer);
        size_t bf_idx = hash % _size;
        size_t rank_searched = _brank(bf_idx+1);
        vector<int> index_res;
 	int count =0;
	int pos=0;
	int end_pos=0;
	int k=0;
	int tmp=0;
	int pos_tmp = 0;
	int h = 0;


	
   	//construction of the int_vector<> that contains all indexes
	for(int i=0; i< set_index.size() ;i++){
		k = index_for_bv[i];
		if( k != 0)
			//sorting of the indexes of each k-mer, before the storage in the int_vector
			std::sort(begin(set_index[i]), begin(set_index[i])+k, std::less<int>());
	}


	//storage of indexes in the int_vector
	for(int i=0; i< set_index.size();i++){
		for(int j =0; j < set_index[i].size() && set_index[i][j] != 0;j++){
			index_kmer[tmp]=set_index[i][j];
			tmp++;
		}	
	}



	for(int i = 1; i < index_for_bv.size(); i++){
		for(int j = pos_tmp; j < bv.size() && j < pos_tmp + index_for_bv[i];j++){
			bv[j]=0;
			h = j;
		}
		bv[h] = 1;
		pos_tmp = h+1;
	}
			
	//search in the bit vector for the start position and end position in index_kmer of a sequence of indexes
	if (_bf[bf_idx]){
		if(rank_searched == 1){
			pos = 0;
			for(k=0; bv[k] == 0 ;k++ );
			end_pos = k;
		}else{
			for(k=0; count != rank_searched-1 && k < bv.size() ;k++ ){
				pos = k;
				if(bv[k] == 1)
					count++;
			}
			pos = k;
			while(bv[k]!=1)
				k++;
		        end_pos=k;
		}
		
		//output vector
		for(int i=pos; i <= end_pos; i++)
			index_res.push_back(index_kmer[i]);

                return index_res;
        }
    }  

private:
    BF() {}
    const BF &operator=(const BF &other) { return *this; }
    const BF &operator=(const BF &&other) { return *this; }


    size_t _size;

    bit_vector _bf;
    rank_support_v<1> _brank;
    rank_support_v<1> _rank;
    vector<int> index_for_bv;
    bit_vector bv;
    vector< vector<int> > set_index;
    int_vector<64> index_kmer;

};

#endif
