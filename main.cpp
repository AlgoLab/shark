#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

#include <zlib.h>

#include "kseq.h"
#include "sdsl/int_vector.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/util.hpp"

#include "argument_parser.hpp"
#include "bloomfilter.h"

using namespace std;

KSEQ_INIT(gzFile, gzread)

auto start_t = chrono::high_resolution_clock::now();

void pelapsed(const string &s = "") {
  auto now_t = chrono::high_resolution_clock::now();
  cerr << "[read-filter/" << s << "] Time elapsed "
       << chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()/1000
       << endl;
}

void analyze_read(const kseq_t *seq, BF &bloom, map<int, string> legend_ID, const uint &k, const uint &c) {
  map<int, int> classification_id;
  string read_seq = seq->seq.s;

  if (read_seq.size() >= k) {
    string kmer (read_seq, 0, k);
    transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
    IDView id_kmer = bloom.get_index(kmer);
    while (id_kmer.has_next())
      ++classification_id[id_kmer.get_next()];

    for (uint p = k; p < read_seq.size(); ++p) {
      char c = toupper(read_seq[p]);
      kmer.erase(0, 1);
      kmer += c;
      id_kmer = bloom.get_index(kmer);
      while (id_kmer.has_next())
        ++classification_id[id_kmer.get_next()];
    }
  }

  int max = std::max(max_element(begin(classification_id),
                                 end(classification_id),
                                 [](const pair<int, int> &a,
                                    const pair<int, int> &b) {
                                   return a.second < b.second;
                                 })
                     ->second,
                     (int)3);
  // save in file all headers of transcripts probably associated to the read
  // search and store in a file all elements of classification with
  // max(classification[i])
  for (auto it_class = classification_id.cbegin();
       it_class != classification_id.cend(); ++it_class) {
    if (it_class->second == max && it_class->second >= opt::c) {
      // legendid[it_class->first] is the name of the transcript, mapped
      /// with index it_class->first
      cout << seq->name.s << "\t" << legend_ID.at(it_class->first) << endl;
    }
  }
}

/*****************************************
 * Main
 *****************************************/
int main(int argc, char *argv[]) {
  parse_arguments(argc, argv);

  /*** 0. Check input files and initialize variables **************************/
  // Transcripts
  gzFile transcript_file = gzopen(opt::fasta_path.c_str(), "r");
  kseq_t *seq = kseq_init(transcript_file);

  // Sample 1
  gzFile read1_file = gzopen(opt::sample1_path.c_str(), "r");
  seq = kseq_init(read1_file);

  // Sample 2
  gzFile read2_file;
  if(opt::paired_flag) {
    read2_file = gzopen(opt::sample2_path.c_str(), "r");
    seq = kseq_init(read2_file);
  }

  BF bloom(opt::bf_size);
  map<int, string> legend_ID;
  int file_line;

  if(opt::verbose) {
    cerr << "Transcripts: " << opt::fasta_path << endl;
    cerr << "Sample 1: " << opt::sample1_path << endl;
    if(opt::paired_flag)
      cerr << "Sample 2: " << opt::sample2_path << endl;
    cerr << "K-mer length: " << opt::k << endl;
    cerr << "Threshold value: " << opt::c << endl << endl;
  }

  /****************************************************************************/

  /*** 1. First iteration over transcripts ************************************/
  // open and read the .fa
  seq = kseq_init(transcript_file);
  int mapped_ID = 0;
  while ((file_line = kseq_read(seq)) >= 0) {
    // seq name is the key map for the transcript, and has an assigned int
    string input_name = seq->name.s;
    string input_seq = seq->seq.s;
    vector<string> transcript_kmers_vec;

    if (input_seq.size() >= opt::k) {
      legend_ID[mapped_ID] = input_name;

      // Build kmers and add them to bf
      string kmer (input_seq, 0, opt::k);
      transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
      bloom.add_kmer(kmer);
      for (uint p = opt::k; p < input_seq.size(); ++p) {
        char c = toupper(input_seq[p]);
        kmer.erase(0, 1);
        kmer += c;
        bloom.add_kmer(kmer);
      }

      ++mapped_ID;
    }
  }

  kseq_destroy(seq);
  gzclose(transcript_file);

  pelapsed("Transcript file processed (" + to_string(mapped_ID) + ")");

  bloom.switch_mode(1);

  /****************************************************************************/

  /*** 2. Second iteration over transcripts ***********************************/
  transcript_file = gzopen(opt::fasta_path.c_str(), "r");
  seq = kseq_init(transcript_file);
  int idx = 0;
  // open and read the .fa, every time a kmer is found the relative index is
  // added to BF
  while ((file_line = kseq_read(seq)) >= 0) {
    string input_seq = seq->seq.s;
    vector<string> transcript_kmers_vec;

    if (input_seq.size() >= opt::k) {
      // Build kmers and store their idx in the bf
      string kmer (input_seq, 0, opt::k);
      transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
      bloom.add_to_kmer(kmer, idx);
      for (uint p = opt::k; p < input_seq.size(); ++p) {
        char c = toupper(input_seq[p]);
        kmer.erase(0, 1);
        kmer += c;
        bloom.add_to_kmer(kmer, idx);
      }

      ++idx;
    }
  }

  kseq_destroy(seq);
  gzclose(transcript_file);

  pelapsed("BF created from transcripts");

  bloom.switch_mode(2);

  /****************************************************************************/

  /*** 3. Iteration over first sample *****************************************/
  // open file that contains the reads
  seq = kseq_init(read1_file);
  while ((file_line = kseq_read(seq)) >= 0) {
    analyze_read(seq, bloom, legend_ID, opt::k, opt::c);
  }

  kseq_destroy(seq);
  gzclose(read1_file);

  pelapsed("First sample completed");

  /****************************************************************************/

  /*** 4. Iteration over second sample ****************************************/
  if(opt::paired_flag){
    read2_file = gzopen(opt::sample2_path.c_str(), "r");

    // open file that contains the reads
    seq = kseq_init(read2_file);
    while ((file_line = kseq_read(seq)) >= 0) {
      analyze_read(seq, bloom, legend_ID, opt::k, opt::c);
    }

    kseq_destroy(seq);
    gzclose(read2_file);

    pelapsed("Second sample completed");
  }

  /****************************************************************************/

  pelapsed("Association done");
  return 0;
}
