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
    BF(const size_t size) :  _bf(size, 0) { _size = size; } ;  
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
        cout << bf_idx<<endl;
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
	_rank = rank_support_v<1>(&bv);
        cout << "A bit vector of length "<< bv.size() << " was created." << endl;
    }
  

    bool add_to_kmer(const string &kmer, ifstream& input) {

    	if (!_mode)
            return false;

        uint64_t hash = _get_hash(kmer);
	int pos;
        size_t bf_idx = hash % _size;
	int i = 0;
	int key = bf_idx;
	cout << key << " sono la chiave del kmer" << kmer << endl;
	string line;
	int cont = 0;
	int newrank = 0;

	cout << "sono la dim di bv " << bv.size() << endl;
	cout << "sono bv " << endl;

	/*for (i =0; i< bv.size(); i++){
		cout << "stampo bv" << endl;
		cout << bv[i] << " "; 
    	}
*/
        if (_bf[bf_idx]) {
        cout << bf_idx << endl;
        size_t cnts_idx = _brank(bf_idx);
	cout << " stampo il rank del kmer considerato" << endl;
        cout << cnts_idx << endl;
        int num_idx=0;

	if (input.is_open()) {

        	while (getline(input, line)) {
			cout << " sto leggendo il file " << line << endl; 
	       		int a = std::stoi(line);
               		index.push_back(a);
			num_idx++;
			}
	}
	input.close();
	cout << "numero di indici letti " << num_idx << endl;
	cout << "chiudo il file " << endl;
        cout << bv.size() << " size bv " << endl;
	
	_rank = rank_support_v<1>(&bv);
	newrank=_rank(bv.size());

	cout << "il rank è " << newrank << endl; 

	if (newrank != 0){
		for(int i=0; cont != newrank && i < bv.size() ;i++ ){
//			cout << "sono i " << i << endl;
//			cout << "sono nel for che cerca la posizione giusta in bv" << endl;
//			cout << "sono rank " << newrank << endl;
//			cout << "sono cont " << cont << endl;
			pos = i;
			if(bv[i] == 1)
				cont++;
		}
//		cout << "sono fuori dal for" << endl;
//		cout << "sono i " << i << endl;
//		cout << "sono pos " << pos << endl;
		
	}else
		pos = 0;

	cout << pos << " sono la posizione cercata" << endl;
	pos = pos+num_idx;
	cout << "sono le dim di bv " << bv.size() << endl;
	cout << pos << " sono la posizione finale degli indici" << endl;

	if(bv[pos]!=1)
		bv[pos] = 1;
	cout << "ho modificato bv " << endl;
	cout << index.size() << " sono la dim di index" << endl;

/*	for(int j =0; j < index.size();j++){
		cout << " stampo l'array di indici" << endl;
		cout << index[j] << " ";
	}
*/
/*	for (i =0; i< bv.size(); i++){
		cout << "stampo bv" << endl;
		cout << bv[i] << " "; 
    	}
*/
        return true;
	}
    }

    vector<int> get_index(const string &kmer){
    	cout << " sono in get new" << endl;
    	uint64_t hash = _get_hash(kmer);
        size_t bf_idx = hash % _size;
        size_t rank_searched = _brank(bf_idx+1);
        vector<int> index_res;
 	int cont =0;
	int pos;
	int i;
	int pos_fin;
	cout << "sono il kmer che voglio cercare " << kmer << endl;
	/*    for(i=0; i< bv.size(); i++)
		cout << bv[i] << endl;
	    for(int j =0; j < index.size();j++){
		cout << " stampo l'array di indici" << endl;
		cout << index[j] << " ";
	    }
	  */
	cout << "sono il rank del kmer che si sta cercando " << rank_searched << endl;
        if (_bf[bf_idx]){
		if(rank_searched ==1){
			pos = 0;
			for(i=0; bv[i] ==0 ;i++ );
			pos_fin =i;
			cout << pos_fin << " sono pos fin di rank 0" << endl;
		}else{
			for(i=0; cont != rank_searched-1 && i < bv.size() ;i++ ){
				cout << "sono i " << i << endl;
				cout << "sono nel for che cerca la posizione giusta in bv" << endl;
				cout << "sono rank " << rank_searched-1 << endl;
				cout << "sono cont " << cont << endl;
				pos = i;
				if(bv[i] == 1)
					cont++;
			}
			while(bv[i]!=1)
				i++;
		        pos_fin=i;
		}
		cout << "sono pos " << pos << endl;
		cout << "sono pos fin " << pos_fin << endl;
		for(int i=pos; i <pos_fin; i++)
			index_res.push_back(index[i]);
		for(int i =0; i<index_res.size();i++)
			cout << index_res[i] << endl;
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
    map<int, bit_vector> myMap;
    vector<int> index;
    bit_vector bv;


};

#endif
