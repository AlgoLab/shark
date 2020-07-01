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
#include <thread>

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


void reference_1st_pass(FastaSplitter& fs, KmerBuilder& kb, BloomfilterFiller& bff) {
  while (true) {
    vector<pair<string, string>>* r_fs = fs();
    if (r_fs == nullptr) return;
    vector<uint64_t>* r_kb = kb(r_fs);
    bff(r_kb);
  }
}

void read_analysis(FastqSplitter& fs, ReadAnalyzer& ra, ReadOutput& ro) {
  FastqSplitter::output_t reads;
  ReadAnalyzer::output_t associations;
  while (true) {
    fs(reads);
    if (reads.empty()) return;
    ra(reads, associations);
    ro(associations);
    reads.clear();
    associations.clear();
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

    FastaSplitter fs(refseq, 100, &legend_ID);
    KmerBuilder kb(opt::k);
    BloomfilterFiller bff(&bloom);

    std::vector<std::thread> threads;
    while (static_cast<int>(threads.size()) < opt::nThreads)
      threads.emplace_back(reference_1st_pass, std::ref(fs), std::ref(kb), std::ref(bff));
    for (auto& t: threads)
      t.join();

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
  vector<uint64_t> kmers;
  while ((seq_len = kseq_read(seq)) >= 0) {
    kmers.clear();

    if ((uint)seq_len >= opt::k) {
      int _p = 0;
      uint64_t kmer = build_kmer(seq->seq.s, _p, opt::k);
      if(kmer == (uint64_t)-1) continue;
      uint64_t rckmer = revcompl(kmer, opt::k);
      kmers.push_back(min(kmer, rckmer));
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
        kmers.push_back(min(kmer, rckmer));
      }
      bloom.add_to_kmer(kmers, nidx);
    }
    ++nidx;
  }
  kseq_destroy(seq);
  gzclose(ref_file);

  pelapsed("BF created from transcripts (" + to_string(nidx) + " genes)");

  bloom.switch_mode(2);
  pelapsed("Second switch performed");

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

    FastqSplitter fs(sseq1, sseq2, 50000, opt::min_quality, out1 != nullptr);
    ReadAnalyzer ra(&bloom, legend_ID, opt::k, opt::c, opt::single);
    ReadOutput ro(out1, out2);

    std::vector<std::thread> threads;
    while (static_cast<int>(threads.size()) < opt::nThreads)
      threads.emplace_back(read_analysis, std::ref(fs), std::ref(ra), std::ref(ro));
    for (auto& t: threads)
      t.join();

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
