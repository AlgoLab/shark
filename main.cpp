#include <iostream>
#include <fstream>
#include "sdsl/int_vector.hpp"
#include "sdsl/util.hpp"
#include "sdsl/int_vector.hpp" // for the bit_vector class
#include <vector>       // std::vector
#include <algorithm>
#include "bloomfilter.h"

using namespace std;
using namespace sdsl;
const size_t sizebloom = 100;
BF bloom(sizebloom);

typedef bit_vector::size_type size_type;


template <typename T1, typename T2>
struct less_second {
    typedef pair<T1, T2> type;
    bool operator ()(type const& a, type const& b) const {
        return a.second < b.second;
    }
};

//function search, returns indexes in a vector
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
	vector<string> kmer_vector;
	ifstream kmer_file("input/kmer.txt");
	ifstream input;
	ifstream input2("input/input2.txt");
	std::string line;
	
	//add k-mer to bloom filter
	while (std::getline(kmer_file, line)){
		bloom.add_kmer(line);
		kmer_vector.push_back(line);
	}
	kmer_file.close();
		
	//add k-mer indexes to data structure
	for(int j = 0; j < kmer_vector.size(); j++){
		string kmer_ref = kmer_vector[j];
		input.open("input/input"+kmer_ref+".txt");
		bloom.add_to_kmer(kmer_ref,input);
	}
	
	//testing with some k-mers
	bloom.add_to_kmer("ATC", input2);
	cout <<  "searching... " << endl;
	output = search("ATC", bloom);

	//print the indexes vector
	if(!output.empty()){
		for(int j = 0; j < output.size(); j++)
			cout << output[j] << " " ;

	}else
		cout << "il kmer non è presente nel bloom filter" << endl;
	
}

