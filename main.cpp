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
#include "sdsl/util.hpp"

#include "argument_parser.hpp"
#include "bloomfilter.h"
#include "BloomfilterFiller.hpp"
#include "KmerBuilder.hpp"
#include "FastaSplitter.hpp"
#include "ReadAnalyzer.hpp"
#include "ReadOutput.hpp"
#include "kmer_utils.hpp"

using namespace std;

auto start_t = chrono::high_resolution_clock::now();

void pelapsed(const string &s = "") {
  auto now_t = chrono::high_resolution_clock::now();
  cerr << "[shark/" << s << "] Time elapsed "
       << chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()/1000
       << endl;
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
  kseq_destroy(seq);
  gzclose(ref_file);

  // Sample 1
  gzFile read1_file = gzopen(opt::sample1_path.c_str(), "r");
  seq = kseq_init(read1_file);
  kseq_destroy(seq);
  gzclose(read1_file);

  // Sample 2
  gzFile read2_file;
  if(opt::paired_flag) {
    read2_file = gzopen(opt::sample2_path.c_str(), "r");
    seq = kseq_init(read2_file);
    kseq_destroy(seq);
    gzclose(read2_file);
  }

  BF bloom(opt::bf_size);
  vector<string> legend_ID;
  int seq_len;

  if(opt::verbose) {
    cerr << "Reference texts: " << opt::fasta_path << endl;
    cerr << "Sample 1: " << opt::sample1_path << endl;
    if(opt::paired_flag)
      cerr << "Sample 2: " << opt::sample2_path << endl;
    cerr << "K-mer length: " << opt::k << endl;
    cerr << "Threshold value: " << opt::c << endl;
    cerr << "Only single associations: " << (opt::single ? "Yes" : "No") << endl;
    cerr << "Minimum base quality: " << static_cast<int>(opt::min_quality) << endl;
    cerr << endl;
  }

  /****************************************************************************/

  /*** 1. First iteration over transcripts ************************************/
  {
    ref_file = gzopen(opt::fasta_path.c_str(), "r");
    kseq_t *refseq = kseq_init(ref_file);
    tbb::filter_t<void, vector<pair<string, string>>*>
      tr(tbb::filter::serial_in_order, FastaSplitter(refseq, 100));
    tbb::filter_t<vector<pair<string, string>>*, vector<uint64_t>*>
      kb(tbb::filter::parallel, KmerBuilder(opt::k));
    tbb::filter_t<vector<uint64_t>*, void>
      bff(tbb::filter::serial_out_of_order, BloomfilterFiller(&bloom));

    tbb::filter_t<void, void> pipeline = tr & kb & bff;
    tbb::parallel_pipeline(opt::nThreads, pipeline);

    kseq_destroy(refseq);
    gzclose(ref_file);
  }

  pelapsed("Transcript file processed");

  bloom.switch_mode(1);

  pelapsed("First switch performed");
  /****************************************************************************/
                                                                        \
  /*** 2. Second iteration over transcripts ***********************************/
  ref_file = gzopen(opt::fasta_path.c_str(), "r");
  seq = kseq_init(ref_file);
  int nidx = 0;
  // open and read the .fa, every time a kmer is found the relative index is
  // added to BF
  while ((seq_len = kseq_read(seq)) >= 0) {
    string input_name = seq->name.s;
    legend_ID.push_back(input_name);

    if ((uint)seq_len >= opt::k) {
      int _p = 0;
      uint64_t kmer = build_kmer(seq->seq.s, _p, opt::k);
      if(kmer == (uint64_t)-1) continue;
      uint64_t rckmer = revcompl(kmer, opt::k);
      bloom.add_to_kmer(min(kmer, rckmer), nidx);
      for (int p = _p; p < seq_len; ++p) {
        uint8_t new_char = to_int[seq->seq.s[p]];
        if(new_char == 0) { // Found a char different from A, C, G, T
          ++p; // we skip this character then we build a new kmer
          kmer = build_kmer(seq->seq.s, p, opt::k);
          if(kmer == (uint64_t)-1) break;
          rckmer = revcompl(kmer, opt::k);
          --p; // p must point to the ending position of the kmer, it will be incremented by the for
        } else {
          --new_char; // A is 1 but it should be 0
          kmer = lsappend(kmer, new_char, opt::k);
          rckmer = rsprepend(rckmer, reverse_char(new_char), opt::k);
        }
        bloom.add_to_kmer(min(kmer, rckmer), nidx);
      }
    }
    ++nidx;
  }
  kseq_destroy(seq);
  gzclose(ref_file);

  pelapsed("BF created from transcripts (" + to_string(nidx) + " genes)");

  bloom.switch_mode(2);
  pelapsed("Second switch performed");

  /****************************************************************************/

  /*** 3. Iteration over first sample *****************************************/
  {
    gzFile f = gzopen(opt::sample1_path.c_str(), "r");
    kseq_t *sseq = kseq_init(f);
    tbb::filter_t<void, vector<pair<string, string>>*>
      sr(tbb::filter::serial_in_order, FastaSplitter(sseq, 50000, opt::min_quality));
    tbb::filter_t<vector<pair<string, string>>*, vector<array<string, 4>>*>
      ra(tbb::filter::parallel, ReadAnalyzer(&bloom, legend_ID, opt::k, opt::c, opt::single));
    tbb::filter_t<vector<array<string, 4>>*, void>
      so(tbb::filter::serial_out_of_order, ReadOutput());

    tbb::filter_t<void, void> pipeline_reads = sr & ra & so;
    tbb::parallel_pipeline(opt::nThreads, pipeline_reads);

    kseq_destroy(sseq);
    gzclose(f);
  }
  pelapsed("First sample completed");

  /****************************************************************************/

  /*** 4. Iteration over second sample ****************************************/
  if(opt::paired_flag){
    read2_file = gzopen(opt::sample2_path.c_str(), "r");

    kseq_t *sseq = kseq_init(read2_file);
    tbb::filter_t<void, vector<pair<string, string>>*>
      sr2(tbb::filter::serial_in_order, FastaSplitter(sseq, 50000, opt::min_quality));
    tbb::filter_t<vector<pair<string, string>>*, vector<array<string, 4>>*>
      ra2(tbb::filter::parallel, ReadAnalyzer(&bloom, legend_ID, opt::k, opt::c, opt::single));
    tbb::filter_t<vector<array<string, 4>>*, void>
      so2(tbb::filter::serial_out_of_order, ReadOutput());
    tbb::filter_t<void, void> pipeline_reads2 = sr2 & ra2 & so2;
    tbb::parallel_pipeline(opt::nThreads, pipeline_reads2);

    kseq_destroy(sseq);
    gzclose(read2_file);

    pelapsed("Second sample completed");
  }

  /****************************************************************************/

  pelapsed("Association done");
  return 0;
}
