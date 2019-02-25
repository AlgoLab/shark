#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_set>

#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#include "sdsl/int_vector.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/util.hpp"

#include "argument_parser.hpp"
#include "bloomfilter.h"
#include "BloomfilterFiller.hpp"
#include "KmerBuilder.hpp"
#include "FastaSplitter.hpp"
#include "ReadAnalyzer.hpp"
#include "ReadOutput.hpp"

using namespace std;

auto start_t = chrono::high_resolution_clock::now();

void pelapsed(const string &s = "") {
  auto now_t = chrono::high_resolution_clock::now();
  cerr << "[read-filter/" << s << "] Time elapsed "
       << chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()/1000
       << endl;
}

void analyze_read(const kseq_t *seq, BF &bloom, const vector<string> &legend_ID, const uint &k, const int &c) {
  map<int, int> classification_id;
  string read_seq = seq->seq.s;

  if (read_seq.size() >= k) {
    string kmer (read_seq, 0, k);
    transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
    IDView id_kmer = bloom.get_index(kmer);
    while (id_kmer.has_next()) {
      int x = id_kmer.get_next();
      ++classification_id[x];
    }

    for (uint p = k; p < read_seq.size(); ++p) {
      char c = toupper(read_seq[p]);
      kmer.erase(0, 1);
      kmer += c;
      id_kmer = bloom.get_index(kmer);
      while (id_kmer.has_next()) {
        int x = id_kmer.get_next();
        ++classification_id[x];
      }
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

  if(max >= c) {
    unordered_set<string> output_names;
    for (auto it_class = classification_id.cbegin(); it_class != classification_id.cend(); ++it_class) {
      if (it_class->second == max) {
        string associated_name = legend_ID[it_class->first];
        if(output_names.find(associated_name) == output_names.end()) {
          // legendid[it_class->first] is the name of the transcript, mapped
          /// with index it_class->first
          output_names.insert(associated_name);
          cout << seq->name.s << " " << associated_name << " " << read_seq << endl;
        }
      }
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
  gzFile ref_file = gzopen(opt::fasta_path.c_str(), "r");
  kseq_t *seq = kseq_init(ref_file);

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
  vector<string> legend_ID;
  int file_line;

  if(opt::verbose) {
    cerr << "Reference texts: " << opt::fasta_path << endl;
    cerr << "Sample 1: " << opt::sample1_path << endl;
    if(opt::paired_flag)
      cerr << "Sample 2: " << opt::sample2_path << endl;
    cerr << "K-mer length: " << opt::k << endl;
    cerr << "Threshold value: " << opt::c << endl << endl;
  }

  /****************************************************************************/

  /*** 1. First iteration over transcripts ************************************/
  // tbb::filter_t<void, vector<pair<string, string>>*>
  //   tr(tbb::filter::serial_in_order, FastaSplitter(opt::fasta_path, 100));
  // tbb::filter_t<vector<pair<string, string>>*, vector<uint64_t>*>
  //   kb(tbb::filter::parallel, KmerBuilder(opt::k));
  // tbb::filter_t<vector<uint64_t>*, void>
  //   bff(tbb::filter::serial_out_of_order, BloomfilterFiller(&bloom));

  // tbb::filter_t<void, void> pipeline = tr & kb & bff;
  // tbb::parallel_pipeline(opt::nThreads, pipeline);

  seq = kseq_init(ref_file);
  while ((file_line = kseq_read(seq)) >= 0) {
    string input_seq = seq->seq.s;

    if (input_seq.size() >= opt::k) {
      // Build kmers and add them to bf
      string kmer (input_seq, 0, opt::k);
      // transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
      bloom.add_kmer(kmer);
      for (uint p = opt::k; p < input_seq.size(); ++p) {
        char c = input_seq[p]; //toupper(input_seq[p]);
        kmer.erase(0, 1);
        kmer += c;
        bloom.add_kmer(kmer);
      }
    }
  }

  kseq_destroy(seq);
  gzclose(ref_file);

  pelapsed("Transcript file processed");

  bloom.switch_mode(1);

  /****************************************************************************/

  /*** 2. Second iteration over transcripts ***********************************/
  ref_file = gzopen(opt::fasta_path.c_str(), "r");
  seq = kseq_init(ref_file);
  int nidx = 0;
  // open and read the .fa, every time a kmer is found the relative index is
  // added to BF
  while ((file_line = kseq_read(seq)) >= 0) {
    string input_name = seq->name.s;
    string input_seq = seq->seq.s;

    if (input_seq.size() >= opt::k) {
      // Build kmers and store their nidx in the bf
      string kmer (input_seq, 0, opt::k);
      // transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
      bloom.add_to_kmer(kmer, nidx);
      for (uint p = opt::k; p < input_seq.size(); ++p) {
        char c = input_seq[p]; //toupper(input_seq[p]);
        kmer.erase(0, 1);
        kmer += c;
        bloom.add_to_kmer(kmer, nidx);
      }
    }

    legend_ID.push_back(input_name);
    ++nidx;
  }

  kseq_destroy(seq);
  gzclose(ref_file);

  pelapsed("BF created from transcripts (" + to_string(nidx+1) + " genes)");

  bloom.switch_mode(2);

  /****************************************************************************/

  /*** 3. Iteration over first sample *****************************************/
  gzFile f = gzopen(opt::sample1_path.c_str(), "r");
  kseq_t *sseq = kseq_init(f);
  tbb::filter_t<void, vector<pair<string, string>>*>
    sr(tbb::filter::serial_in_order, FastaSplitter(sseq, 10000));
  tbb::filter_t<vector<pair<string, string>>*, vector<array<string, 3>>*>
    ra(tbb::filter::parallel, ReadAnalyzer(&bloom, legend_ID, opt::k, opt::c));
  tbb::filter_t<vector<array<string, 3>>*, void>
    so(tbb::filter::serial_out_of_order, ReadOutput());

  tbb::filter_t<void, void> pipeline_reads = sr & ra & so;
  tbb::parallel_pipeline(opt::nThreads, pipeline_reads);

  // seq = kseq_init(read1_file);
  // while ((file_line = kseq_read(seq)) >= 0) {
  //   analyze_read(seq, bloom, legend_ID, opt::k, opt::c);
  // }

  kseq_destroy(sseq);
  gzclose(f);

  pelapsed("First sample completed");

  /****************************************************************************/

  /*** 4. Iteration over second sample ****************************************/
  if(opt::paired_flag){
    read2_file = gzopen(opt::sample2_path.c_str(), "r");

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
