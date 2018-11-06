#include "bloomfilter.h"
#include "kseq.h"
#include "sdsl/int_vector.hpp"
#include "sdsl/int_vector.hpp" // for the bit_vector class
#include "sdsl/util.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
//#include <seqan/align.h>
#include <stdio.h>
#include <string>
#include <vector> // std::vector
#include <zlib.h>

static size_t sizebloom = (10e8);

// function search, returns indexes in a vector
IDView search(const string &kmer, BF &bloomfilter) {
  return bloomfilter.get_index(kmer);
}
// using namespace seqan;

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

/*****************************************
 * Main
 *****************************************/

int main(int argc, char *argv[]) {
  gzFile transcript_file;
  kseq_t *seq;
  long file_line;
  map<long, string> legend_ID;
  long mapped_ID = 0;
  const int kmer_length = 60;
  vector<string> transcript_kmers_vec;
  BF bloom(sizebloom);
  string name_transcript;
  string transcript_name = "";
  string read_name = "";
  // typedef String<char> TSequence;             // sequence type

  time_t start = time(0);

  if (argc > 1) {
    transcript_name = argv[1];
    read_name = argv[2];
  } else {
    cout << "Error in input" << endl;
    return 1;
  }

  time_t read_time = time(0);

  transcript_file =
      gzopen(transcript_name.c_str(), "r"); // STEP 2: open the file handler
  seq = kseq_init(transcript_file);         // STEP 3: initialize seq
  long size = 0;
  string input_seq;
  // open and read the .fa
  while ((file_line = kseq_read(seq)) >=
         0) { // STEP 4: read sequence of transcript

    // seq name is the key map for the transcript, and has an assigned int
    name_transcript = seq->name.s;
    input_seq = seq->seq.s;
    //  cout << name_transcript << endl;

    if (input_seq.size() >= kmer_length) {
      legend_ID[mapped_ID] = name_transcript;
      // split each sequence of a transcritp in k-mer with k=n

      transcript_kmers_vec.resize(input_seq.size() - (kmer_length - 1));
      transform(input_seq.cbegin(), input_seq.cend() - (kmer_length - 1),
                transcript_kmers_vec.begin(), [kmer_length](const auto &i) {
                  return string(&i, kmer_length);
                });

      // add all k-mers to BF

      for (auto &kmer : transcript_kmers_vec) {
        bloom.add_kmer(kmer);
      }
      size = size + transcript_kmers_vec.size();
      //    cout <<size <<endl;
      mapped_ID++;
      //    cout << "map " << mapped_ID << endl;

      transcript_kmers_vec.clear();
    }
    input_seq.clear();
  }

  printf("return value: %ld\n", file_line);
  kseq_destroy(seq);        // STEP 5: destroy seq
  gzclose(transcript_file); // STEP 6: close the file handler

  cout << "Transcript file processed" << endl;

  int seconds = 0;
  int minutes = 0;
  time_t end_read_time = time(0);
  int time_reading = difftime(end_read_time, read_time) * 1000.0;

  seconds = (time_reading / 1000) % 60;
  minutes = (int)((time_reading / (1000 * 60)) % 60);

  cout << "Time to process transcripts:  " << minutes << ":" << seconds << endl;

  bloom.switch_mode(1);

  time_t add_time = time(0);

  transcript_file = gzopen(transcript_name.c_str(), "r");
  seq = kseq_init(transcript_file);
  long idx = 0;
  // open and read the .fa, every time a kmer is found the relative index is
  // added to BF
  while ((file_line = kseq_read(seq)) >= 0) {
    name_transcript = seq->name.s;
    input_seq = seq->seq.s;

    if (input_seq.size() >= kmer_length) {

      // split each sequence of a transcritp in k-mer with k=n

      transcript_kmers_vec.resize(input_seq.size() - (kmer_length - 1));
      transform(input_seq.cbegin(), input_seq.cend() - (kmer_length - 1),
                transcript_kmers_vec.begin(), [kmer_length](const auto &i) {
                  return string(&i, kmer_length);
                });

      // add for each k-mer its id to BF
      for (auto &kmer : transcript_kmers_vec)
        bloom.add_to_kmer(kmer, idx);

      idx++;
    }
  }

  printf("return value: %ld\n", file_line);
  kseq_destroy(seq);        // STEP 5: destroy seq
  gzclose(transcript_file); // STEP 6: close the file handler

  cout << "Transcript indexes added to Bloom filter" << endl;
  time_t end_add_time = time(0);
  int time_adding = difftime(end_add_time, add_time) * 1000.0;

  seconds = (int)((time_adding / 1000) % 60);
  minutes = (int)((time_adding / (1000 * 60)) % 60);

  cout << "Time to add indexes:  " << minutes << ":" << seconds << endl;

  bloom.switch_mode(2);

  time_t association_time = time(0);
  gzFile read_file;
  string read_seq;
  vector<string> read_kmers_vec;
  map<long, long> classification_id;
  ofstream final_id("id_results.fa");
  // vector<int> id_kmer;
  IDView id_kmer;

  // FILE * pFile;
  // pFile = fopen ("myfile.txt", "wb");

  // open .fq file that contains the reads
  read_file = gzopen(read_name.c_str(), "r"); // STEP 2: open the file handler
  seq = kseq_init(read_file);                 // STEP 3: initialize seq

  while ((file_line = kseq_read(seq)) >= 0) { // STEP 4: read sequence of reads

    read_seq = seq->seq.s;

    if (read_seq.size() >= kmer_length) {

      read_kmers_vec.resize(read_seq.size() - (kmer_length - 1));
      transform(read_seq.cbegin(), read_seq.cend() - (kmer_length - 1),
                read_kmers_vec.begin(), [kmer_length](const auto &i) {
                  return string(&i, kmer_length);
                });

      // search all kmers of each read in the BF and store relative indexes

      for (auto &kmer : read_kmers_vec) {
        id_kmer = bloom.get_index(kmer);
        //		int* buffer = &id_kmer[0];
        //		fwrite (buffer, sizeof(int), sizeof(buffer), pFile);
        while (id_kmer.has_next())
          classification_id[id_kmer.get_next()]++;
      }
    }

    long max = std::max(std::max_element(begin(classification_id), end(classification_id),
                                         [](const std::pair<long, long> &a, const std::pair<long, long> &b) {
                                           return a.second < b.second;
                                         })->second,
                        (long)3);
    // save in file all headers of transcripts probably associated to the read
    // search and store in a file all elements of classification with
    // max(classification[i])
    for (auto it_class = classification_id.cbegin();
         it_class != classification_id.cend(); it_class++) {
      if (it_class->second == max) {
        // legend_ID[it_class->first] is the name of the transcript, mapped
        // with index it_class->first
        final_id << ">" << legend_ID[it_class->first] << endl;
        final_id << read_seq << endl;
      }
    }

    classification_id.clear();
    id_kmer.clear();
    read_kmers_vec.clear();
  }

  final_id.close();
  // fclose(pFile);

  printf("return value: %ld\n", file_line);
  kseq_destroy(seq);  // STEP 5: destroy seq
  gzclose(read_file); // STEP 6: close the file handler

  cout << "Association done." << endl;
  time_t end_assoc_time = time(0);
  int time_assoc = difftime(end_assoc_time, association_time) * 1000.0;

  seconds = (int)((time_assoc / 1000) % 60);
  minutes = (int)((time_assoc / (1000 * 60)) % 60);

  cout << "Time to associate reads to transcript:  " << minutes << ":"
       << seconds << endl;

  time_t end = time(0);
  int time_proc = difftime(end, start) * 1000.0;

  seconds = (int)((time_proc / 1000) % 60);
  minutes = (int)((time_proc / (1000 * 60)) % 60);

  cout << "Time of entire process:  " << minutes << ":" << seconds << endl;

  cout << "Reading transcripts: " << (int)(time_reading / time_proc) * 100
       << endl;
  cout << "Adding indexes to BF: " << (int)(time_adding / time_proc) * 100
       << endl;
  cout << "Associating reads to transcript: "
       << (int)(time_assoc / time_proc) * 100 << endl;

  return 0;
}
