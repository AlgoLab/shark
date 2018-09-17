#include <iostream>
#include <fstream>
#include "sdsl/int_vector.hpp"
//#include <SDL/begin_code.h>
#include "sdsl/util.hpp"
#include "sdsl/int_vector.hpp" // for the bit_vector class

#include <vector>       // std::vector
#include <algorithm>
#include <array>
#include <functional>  // for std::greater
#include "sdsl/rank_support.hpp" // for rank data structures

#include "bloomfilter.h"
#include "hts_log.h"
#include "KMC/kmc_api/kmc_file.h"
#include "kseq.h"


using namespace std;
using namespace sdsl;
const size_t sizebloom = 100;
BF bloom(sizebloom);

typedef bit_vector::size_type size_type;
int const bit_size = 64;

template <typename T1, typename T2>
struct less_second {
    typedef pair<T1, T2> type;
    bool operator ()(type const& a, type const& b) const {
        return a.second < b.second;
    }
};

   vector<int> search(const string &kmer, BF &bloomfilter){
	std::vector<int> null;
	if(bloomfilter.test_kmer(kmer)){
		return bloomfilter.get_index(kmer);
	}else
		return null;
  }


/*****************************************
 * Main
 *****************************************/

int main(int argc, char *argv[]) {

	BF bloom(sizebloom);
	vector<int> output;
	vector<int> index_vector;
	ifstream kmer_file("input/kmer.txt");
	ifstream input;
	map<int,string> mappa;
	char filename[64];

	
	//aggiungo al bloom filter i kmer presi da un file contenente i kmer che voglio cercare
	std::string line;
	int i=0;
	while (std::getline(kmer_file, line)){
		//aggiungo il kmer scelto al bloom filter

		bloom.add_kmer(line);
	}

	kmer_file.close();
	kmer_file.open("input/kmer.txt");

	while (std::getline(kmer_file, line)){
		int rank = bloom.get_rank_kmer(line);
		int index = bloom.get_idx(line);
		index_vector.push_back(index);
		mappa.insert(pair<int,string>(index,line));
	}

	//ordino il vettore che contiene gli indici
	std::sort(index_vector.begin(), index_vector.end(), std::less<int>());
	
	//aggiungo tutti i kmer al bloom filter
	for(int j = 0; j < index_vector.size(); j++){
		string kmer_ref = mappa.at(index_vector[j]);
		input.open("input/input"+kmer_ref+".txt");
		bloom.add_to_kmer(kmer_ref,input);	
		
	}
	

	cout << "************************* \n faccio search \n *****************" << endl;
	output = search("CTT", bloom);


	//stampo il vettore restituito dal search
	if(!output.empty()){

		for(int j = 0; j < output.size(); j++){
			cout << output[j] << " " ;
		}
	}else
		cout << "il kmer non è presente nel bloom filter" << endl;

	cout << "*************************" << endl;
	
}

