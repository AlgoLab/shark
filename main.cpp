#include <iostream>
#include <fstream>
#include <sdsl/int_vector.hpp>
#include <SDL/begin_code.h>
#include <sdsl/util.hpp>
#include <sdsl/int_vector.hpp> // for the bit_vector class

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
		cout << "sono nell'if di search" << endl;
		return bloomfilter.get_index(kmer);
	}else
		return null;
  }


/*****************************************
 * Main
 *****************************************/

int main(int argc, char *argv[]) {

	cout << "inizializzo bloom \n";
	BF bloom(sizebloom);
	vector<int> output;
	vector<int> index;
	ifstream kmer_file("kmer.txt");
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
	bloom.switch_mode();
	kmer_file.close();
	kmer_file.open("kmer.txt");

	while (std::getline(kmer_file, line)){
		int rank = bloom.get_rank_kmer(line);
		int indice = bloom.get_idx(line);
		index.push_back(indice);
		mappa.insert(pair<int,string>(indice,line));
	}

	//ordino il vettore che contiene gli indici
	std::sort(index.begin(), index.end(), std::less<int>());
	
	//aggiungo tutti i kmer al bloom filter
	cout << "sto per entrare nel for" << endl;
	for(int j = 0; j < index.size(); j++){
		cout << "sono nel for" << endl;
		cout << mappa.at(index[j]) << endl;
		string kmer_ref = mappa.at(index[j]);
		cout << "faccio increment" << endl;

		input.open("input"+kmer_ref+".txt");
		bloom.add_to_kmer(kmer_ref,input);	
		
	}
	
	cout << "************************* \n faccio search \n *****************" << endl;
	output = search("CGG", bloom);

	//stampo il vettore restituito dal search
	if(!output.empty()){

		for(int j = 0; j < output.size(); j++){
			cout << output[j] << " " ;
		}
	}else
		cout << "il kmer non è presente nel bloom filter" << endl;
	
}

