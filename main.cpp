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
  gzFile transcript_file;
  kseq_t *seq;
  int file_line;
  map<int, string> legend_ID;
  int mapped_ID = 0;

  const int kmer_length = 60;

  vector<string> transcript_kmers_vec;
  BF bloom(sizebloom);
  string name_transcript;

  transcript_file =
      gzopen("example/chrY_mod.fa", "r"); // STEP 2: open the file handler
  seq = kseq_init(transcript_file);       // STEP 3: initialize seq

  // open and read the .fa
  ofstream kmer_idx_FILE;
  kmer_idx_FILE.open("example/kmer_idx.txt");
  while ((file_line = kseq_read(seq)) >=
         0) { // STEP 4: read sequence of transcript
    string input_seq;
    // seq name is the key map for the transcript, and has an assigned int
    name_transcript = seq->name.s;

    legend_ID[mapped_ID] = name_transcript;

    // split each sequence of a transcritp in k-mer with k=n
    input_seq = seq->seq.s;
    transcript_kmers_vec.resize(input_seq.size() - (kmer_length - 1));
    transform(input_seq.cbegin(), input_seq.cend() - (kmer_length - 1),
              transcript_kmers_vec.begin(),
              [kmer_length](const auto &i) { return string(&i, kmer_length); });

    // add all k-mers to BF

    // FIXME: devo creare un file perchè non è possibile switchare tra mod 0 e
    // mod 1 per aggiungere volta per volta prima kmer e poi indici associati
    for (auto &kmer : transcript_kmers_vec) {
      bloom.add_kmer(kmer);
      kmer_idx_FILE << kmer << " " << mapped_ID << " " << endl;
    }

    mapped_ID++;
  }

  kmer_idx_FILE.close();
  printf("return value: %d\n", file_line);
  kseq_destroy(seq);        // STEP 5: destroy seq
  gzclose(transcript_file); // STEP 6: close the file handler

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

  bloom.switch_mode(2);

  gzFile read_file;
  string read_seq;
  vector<string> read_kmers_vec;
  map<int, int> classification_id;
  ofstream final_id;
  vector<int> id_kmer;
  final_id.open("example/final_id.txt");

  // open .fq file that contains the reads
  read_file =
      gzopen("example/chrY_reads.fq", "r"); // STEP 2: open the file handler
  seq = kseq_init(read_file);               // STEP 3: initialize seq

  while ((file_line = kseq_read(seq)) >= 0) { // STEP 4: read sequence of reads

    read_seq = seq->seq.s;

    read_kmers_vec.resize(read_seq.size() - (kmer_length - 1));
    transform(read_seq.cbegin(), read_seq.cend() - (kmer_length - 1),
              read_kmers_vec.begin(),
              [kmer_length](const auto &i) { return string(&i, kmer_length); });

    // search all kmers of each read in the BF and store relative indexes

    for (auto &kmer : read_kmers_vec) {
      id_kmer = bloom.get_index(kmer);

      for (auto &id : id_kmer) {

        classification_id[id]++;
      }
    }

    int max = 0;

    // save in file all headers of transcripts probably associated to the read
    for (auto it_class = classification_id.cbegin();
         it_class != classification_id.cend(); it_class++) {
      // it_class->second is the number of times that the index
      // (it_class->first) has been found
      if (it_class->second >= max) {

        max = it_class->second;
      }
    }
    // search and store in a file all elements of classification with
    // max(classification[i])
    for (auto it_class = classification_id.cbegin();
         it_class != classification_id.cend(); it_class++) {
      if (it_class->second == max) {
        // legend_ID[it_class->first] is the name of the transcript, mapped with
        // index it_class->first
        final_id << read_seq << " " << legend_ID[it_class->first] << endl;
      }
    }

    classification_id.clear();
    id_kmer.clear();
    read_kmers_vec.clear();
  }

  final_id.close();

  printf("return value: %d\n", file_line);
  kseq_destroy(seq);  // STEP 5: destroy seq
  gzclose(read_file); // STEP 6: close the file handler

  return 0;
}
