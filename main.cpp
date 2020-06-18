/**
 * shark - Mapping-free filtering of useless RNA-Seq reads
 * Copyright (C) 2019 Tamara Ceccato, Luca Denti, Yuri Pirola, Marco Previtali
 *
 * This file is part of shark.
 *
 * shark is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * shark is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with shark; see the file LICENSE. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "common.hpp"
#include "argument_parser.hpp"
#include "bloomfilter.h"
#include "BloomfilterFiller.hpp"
#include "KmerBuilder.hpp"
#include "FastaSplitter.hpp"
#include "FastqSplitter.hpp"
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
  gzFile read2_file = nullptr;
  if(opt::paired_flag) {
    read2_file = gzopen(opt::sample2_path.c_str(), "r");
    seq = kseq_init(read2_file);
    kseq_destroy(seq);
    gzclose(read2_file);
  }

  BF bloom(opt::bf_size);
  vector<string> legend_ID;
  legend_ID.reserve(100);

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
    tbb::filter_t<void, vector<string>*>
      tr(tbb::filter::serial_in_order, FastaSplitter(refseq, 100, &legend_ID));
    tbb::filter_t<vector<string>*, vector<vector<uint64_t>>*>
      kb(tbb::filter::parallel, KmerBuilder(opt::k, opt::bf_size));
    tbb::filter_t<vector<vector<uint64_t>>*, void>
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
  /*** 2a. Second iteration over transcripts **********************************/
  {
    ref_file = gzopen(opt::fasta_path.c_str(), "r");
    kseq_t *refseq = kseq_init(ref_file);
    tbb::filter_t<void, vector<string>*>
      tr(tbb::filter::serial_in_order, FastaSplitter(refseq, 100));
    tbb::filter_t<vector<string>*, vector<vector<uint64_t>>*>
      kb(tbb::filter::parallel, KmerBuilder(opt::k, opt::bf_size));

    tbb::filter_t<void, void> pipeline =
      tr &
      kb &
      tbb::make_filter<vector<vector<uint64_t>>*, void>(
                                   tbb::filter::serial_in_order,
                                   [&](vector<vector<uint64_t>>* kmerss) {
                                     if (kmerss == nullptr) return;
                                     for (auto& kmers: *kmerss) {
                                       bloom.add_to_kmer_1(kmers);
                                     }
                                     delete kmerss;
                                   });
    tbb::parallel_pipeline(opt::nThreads, pipeline);

    kseq_destroy(refseq);
    gzclose(ref_file);
  }

  pelapsed("Set sizes computed");
  bloom.switch_mode(2);
  pelapsed("Second switch performed");
  /****************************************************************************/
  /*** 2b. Third iteration over transcripts ***********************************/
  {
    ref_file = gzopen(opt::fasta_path.c_str(), "r");
    kseq_t *refseq = kseq_init(ref_file);
    int nidx = 0;
    tbb::filter_t<void, vector<string>*>
      tr(tbb::filter::serial_in_order, FastaSplitter(refseq, 100));
    tbb::filter_t<vector<string>*, vector<vector<uint64_t>>*>
      kb(tbb::filter::parallel, KmerBuilder(opt::k, opt::bf_size));

    tbb::filter_t<void, void> pipeline =
      tr &
      kb &
      tbb::make_filter<vector<vector<uint64_t>>*, void>(
                                   tbb::filter::serial_in_order,
                                   [&](vector<vector<uint64_t>>* kmerss) {
                                     if (kmerss == nullptr) return;
                                     for (auto& kmers: *kmerss) {
                                       bloom.add_to_kmer_2(kmers, nidx);
                                       ++nidx;
                                     }
                                     delete kmerss;
                                   });
    tbb::parallel_pipeline(opt::nThreads, pipeline);

    kseq_destroy(refseq);
    gzclose(ref_file);
    pelapsed("BF created from transcripts (" + to_string(nidx) + " genes)");
  }

  bloom.switch_mode(3);
  pelapsed("Third switch performed");

  /****************************************************************************/

  /*** 3. Iteration over the sample *****************************************/
  {
    kseq_t *sseq1 = nullptr, *sseq2 = nullptr;
    FILE *out1 = nullptr, *out2 = nullptr;
    read1_file = gzopen(opt::sample1_path.c_str(), "r");
    sseq1 = kseq_init(read1_file);
    if (opt::out1_path != "") {
      out1 = fopen(opt::out1_path.c_str(), "w");
    }
    if(opt::paired_flag) {
      read2_file = gzopen(opt::sample2_path.c_str(), "r");
      sseq2 = kseq_init(read2_file);
      if (opt::out2_path != "") {
        out2 = fopen(opt::out2_path.c_str(), "w");
      }
    }

    tbb::filter_t<void, FastqSplitter::output_t*>
      sr(tbb::filter::serial_in_order, FastqSplitter(sseq1, sseq2, 50000, opt::min_quality, out1 != nullptr));
    tbb::filter_t<FastqSplitter::output_t*, ReadAnalyzer::output_t*>
      ra(tbb::filter::parallel, ReadAnalyzer(&bloom, legend_ID, opt::k, opt::c, opt::single));
    tbb::filter_t<ReadAnalyzer::output_t*, void>
      so(tbb::filter::serial_in_order, ReadOutput(out1, out2));

    tbb::filter_t<void, void> pipeline_reads = sr & ra & so;
    tbb::parallel_pipeline(opt::nThreads, pipeline_reads);

    kseq_destroy(sseq1);
    gzclose(read1_file);
    if(opt::paired_flag) {
      kseq_destroy(sseq2);
      gzclose(read2_file);
    }
    if (out1 != nullptr) fclose(out1);
    if (out2 != nullptr) fclose(out2);
  }
  pelapsed("Sample completed");

  /****************************************************************************/

  pelapsed("Association done");
  return 0;
}
