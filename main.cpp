#include "bloomfilter.h"
#include "kseq.h"
#include "sdsl/int_vector.hpp"
#include "sdsl/int_vector.hpp" // for the bit_vector class
#include "sdsl/util.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector> // std::vector
#include <zlib.h>

const size_t sizebloom = 1000000;
BF bloom(sizebloom);

// function search, returns indexes in a vector
vector<int> search(const string &kmer, BF &bloomfilter) {
  return bloomfilter.get_index(kmer);
}

// STEP 1: declare the type of file handler and the read() function

KSEQ_INIT(gzFile, gzread)
/*****************************************
 * Main
 *****************************************/

int main(int argc, char *argv[]) {
  gzFile fp;
  kseq_t *seq;
  int l;
  map<string, int> MapID;
  int mapped_ID = 0;
  string key_map;
  string input_seq;
  const int n = 5; // Assumed positive number smaller than str.size()
  const int n1 = n - 1;
  vector<string> transcript_kmers;
  BF bloom(sizebloom);

  fp = gzopen("example/chrY_mod.fa", "r"); // STEP 2: open the file handler
  seq = kseq_init(fp);                   // STEP 3: initialize seq

  ofstream kmer_idx_FILE;
  kmer_idx_FILE.open("example/output.txt");
  while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence

    // seq name is the key map for the transcript, and has an assigned int
    key_map = seq->name.s;

    MapID[key_map] = mapped_ID;

    input_seq = seq->seq.s;
    transcript_kmers.resize(input_seq.size() - n1);
    transform(input_seq.cbegin(), input_seq.cend() - n1,
              transcript_kmers.begin(),
              [n](const auto &i) { return string(&i, n); });

    // add all elements of result to BF
    for (auto &kmer : transcript_kmers) {
      bloom.add_kmer(kmer);
      kmer_idx_FILE << kmer << " " << mapped_ID << " " << endl;
    }

    mapped_ID++;
  }

  kmer_idx_FILE.close();

  kseq_destroy(seq);
  gzclose(fp);
  l = 0;
  fp = gzopen("example/chrY_mod.fa", "r");
  seq = kseq_init(fp);

  ifstream FILE_input;
  FILE_input.open("example/output.txt");
  int idx;
  string kmer;
  bloom.switch_mode(1);
  while (FILE_input >> kmer >> idx) {

    bloom.add_to_kmer(kmer, idx);
  }
  FILE_input.close();

  bloom.switch_mode(2);

  vector<int> result;
  result = bloom.get_index("TTTTT");

  for (int &i : result)
    cout << i << " ";
  cout << endl;

  printf("return value: %d\n", l);
  kseq_destroy(seq); // STEP 5: destroy seq
  gzclose(fp);       // STEP 6: close the file handler

  
  return 0;
}
