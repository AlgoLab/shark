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


static size_t sizebloom = (10e11);


// function search, returns indexes in a vector
vector<int> search(const string &kmer, BF &bloomfilter) {
  return bloomfilter.get_index(kmer);
}
//using namespace seqan;

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
  string transcript_name = "";
  string read_name = "";
  //typedef String<char> TSequence;             // sequence type

  if (argc > 1) {
    transcript_name = argv[1];
    read_name = argv[2];
  } else {
    cout << "Error in input" << endl;
    return 1;
  }

  transcript_file =
      gzopen(transcript_name.c_str(), "r"); // STEP 2: open the file handler
  seq = kseq_init(transcript_file);         // STEP 3: initialize seq

  // open and read the .fa
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
    }

    mapped_ID++;
  }

  printf("return value: %d\n", file_line);
  kseq_destroy(seq);        // STEP 5: destroy seq
  gzclose(transcript_file); // STEP 6: close the file handler

cout << "Transcript file processed" << endl;

  bloom.switch_mode(1);

  transcript_file = gzopen(transcript_name.c_str(), "r");
  seq = kseq_init(transcript_file);
  int idx = 0;
  // open and read the .fa, every time a kmer is found the relative index is
  // added to BF
  while ((file_line = kseq_read(seq)) >= 0) {
    string input_seq;

    // split each sequence of a transcritp in k-mer with k=n
    input_seq = seq->seq.s;
    transcript_kmers_vec.resize(input_seq.size() - (kmer_length - 1));
    transform(input_seq.cbegin(), input_seq.cend() - (kmer_length - 1),
              transcript_kmers_vec.begin(),
              [kmer_length](const auto &i) { return string(&i, kmer_length); });

    // add for each k-mer its id to BF
    for (auto &kmer : transcript_kmers_vec)
      bloom.add_to_kmer(kmer, idx);

    idx++;
  }

cout << "Transcript indexes added to Bloom filter" << endl;

  bloom.switch_mode(2);

  gzFile read_file;
  string read_seq;
  vector<string> read_kmers_vec;
  map<int, int> classification_id;
  ofstream final_id;
  vector<int> id_kmer;
  final_id.open("id_results.fa");

  // open .fq file that contains the reads
  read_file = gzopen(read_name.c_str(), "r"); // STEP 2: open the file handler
  seq = kseq_init(read_file);                 // STEP 3: initialize seq

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
        final_id << ">" << legend_ID[it_class->first] << endl;
        final_id << read_seq << endl;
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

cout << "Association done." << endl;

/****************************************
*ALIGNMENT
****************************************/
/*
  // read fasta created as output in the preious step
  gzFile test_file;
  test_file =
      gzopen("final_id.fa", "r"); // STEP 2: open the file handler
  seq = kseq_init(test_file);             // STEP 3: initialize seq
  // seq1--> read
  // seq2--> transcript
  ofstream al_file;
  al_file.open("alignment.fa");
  while ((file_line = kseq_read(seq)) >= 0) { // read final_id
    kseq_t *seq_tran;
    string transcript_assigned;
    transcript_assigned = seq->name.s;
    TSequence seq1 = seq->seq.s;
    gzFile transcript_global;
    string name_tr;
 	TSequence seq2;
    transcript_global = gzopen("example/chrY_mod.fa", "r");
    seq_tran = kseq_init(transcript_global);
    while ((file_line = kseq_read(seq_tran)) >= 0) { // read final_id
      name_tr = seq_tran->name.s;
      if (name_tr == transcript_assigned) {
        seq2 = seq_tran->seq.s;
        Align<String<char>> ali;
        resize(rows(ali), 2);
        assignSource(row(ali, 0), seq1);
        assignSource(row(ali, 1), seq2);
        al_file << ">" << name_tr << endl;
        al_file << ";Score = "
                << localAlignment(ali, Score<int>(3, -3, -2, -2), DynamicGaps())
                << ";" << endl;
        al_file << seq1 << endl;
      }
    }
    kseq_destroy(seq_tran);     // STEP 5: destroy seq
    gzclose(transcript_global); // STEP 6: close the file handler
  }
  al_file.close();
  printf("return value: %d\n", file_line);
  kseq_destroy(seq);  // STEP 5: destroy seq
  gzclose(test_file); // STEP 6: close the file handler
  //read fasta with reads and scores of their alignment with transcript
  gzFile align_fasta;
  align_fasta =
      gzopen("example/alignment.fa", "r"); 
  seq = kseq_init(align_fasta);  
  while ((file_line = kseq_read(seq)) >= 0) { 
  		cout << seq->name.s << endl;
  		cout << seq->comment.l << endl;
  		cout << "seq: " << seq->seq.s << endl;
  }
  printf("return value: %d\n", file_line);
  kseq_destroy(seq);  
  gzclose(align_fasta); 
*/


  return 0;
}