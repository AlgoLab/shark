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
  map<int, string> legend_ID;
  int mapped_ID = 0;
  string key_map;
  string input_seq;
  const int n = 5; // Assumed positive number smaller than str.size()
  const int n1 = n - 1;
  vector<string> transcript_kmers;
  BF bloom(sizebloom);
  string name_transcript;

  fp = gzopen("example/chrY_mod.fa", "r"); // STEP 2: open the file handler
  seq = kseq_init(fp);                     // STEP 3: initialize seq

  // open and read the .fa
  ofstream kmer_idx_FILE;
  kmer_idx_FILE.open("example/output.txt");
  while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence of transcript

    // seq name is the key map for the transcript, and has an assigned int
    name_transcript = seq->name.s;
    legend_ID[mapped_ID] = name_transcript;

    // split each sequence of a transcritp in k-mer with k=n
    input_seq = seq->seq.s;
    transcript_kmers.resize(input_seq.size() - n1);
    transform(input_seq.cbegin(), input_seq.cend() - n1,
              transcript_kmers.begin(),
              [n](const auto &i) { return string(&i, n); });

    // add all k-mers to BF
    for (auto &kmer : transcript_kmers) {
      bloom.add_kmer(kmer);
      kmer_idx_FILE << kmer << " " << mapped_ID << " " << endl;
    }

    mapped_ID++;
  }

  kmer_idx_FILE.close();
  printf("return value: %d\n", l);
  kseq_destroy(seq); // STEP 5: destroy seq
  gzclose(fp);       // STEP 6: close the file handler

  ifstream FILE_input;
  FILE_input.open("example/output.txt");
  int idx;
  string kmer;
  bloom.switch_mode(1);
  // read a file that contains a k-mer and the relative id of the transcript
  while (FILE_input >> kmer >> idx) {
    // add each k-mer and its id to BF
    bloom.add_to_kmer(kmer, idx);
  }
  FILE_input.close();

  // test
  bloom.switch_mode(2);
  vector<int> result;
  result = bloom.get_index("TTTTT");

  // open .fq file that contains the reads
  fp = gzopen("example/chrY_reads.fq", "r"); // STEP 2: open the file handler
  seq = kseq_init(fp);                       // STEP 3: initialize seq
  string read_header;
  string read_seq;
  vector<string> read_kmers;
  map<int, int> classification_id;

  map<string, int> read_mapped;

  while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence of reads
    read_header = seq->name.s;
    read_seq = seq->seq.s;

    input_seq = seq->seq.s;
    transcript_kmers.resize(input_seq.size() - n1);
    transform(input_seq.cbegin(), input_seq.cend() - n1,
              transcript_kmers.begin(),
              [n](const auto &i) { return string(&i, n); });

    read_kmers.resize(read_seq.size() - n1);
    transform(read_seq.cbegin(), read_seq.cend() - n1, read_kmers.begin(),
              [n](const auto &i) { return string(&i, n); });

    // search all kmers of each read in the BF and store relative indexes

    vector<int> id_kmer;

    for (auto &kmer : read_kmers) {
      id_kmer = bloom.get_index(kmer);
      for (auto &id : id_kmer)
        classification_id[id] += 1;
    }
    int index_found;
    int max = 0;
    for (auto it_class = classification_id.cbegin();
         it_class != classification_id.cend(); ++it_class) {

      if (it_class->second > max) {

        max = it_class->second;
        index_found = it_class->first;
      }
    }

    read_mapped[read_seq] = index_found;
    classification_id.clear();
  }

  ofstream final_id;
  final_id.open("example/final_id.txt");
  string header_transcript;
  for (auto iterator_read = read_mapped.cbegin();
       iterator_read != read_mapped.cend(); ++iterator_read) {

    auto iterator_indexes = legend_ID.find(iterator_read->second);

    if (iterator_indexes != legend_ID.end()) {
      header_transcript = iterator_indexes->second;
    }
    final_id << iterator_read->first << " " << header_transcript << endl;
  }
  final_id.close();

  printf("return value: %d\n", l);
  kseq_destroy(seq); // STEP 5: destroy seq
  gzclose(fp);       // STEP 6: close the file handler

  return 0;
}
