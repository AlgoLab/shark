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
#include <mutex>
#include <condition_variable>

#include <zlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "common.hpp"
#include "argument_parser.hpp"
#include "bloomfilter.h"
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


class BFFiller1st {
public:
  BFFiller1st(BF& _bf) : bf(_bf) {}

  void operator()(const vector<uint64_t>& positions) {
    std::lock_guard<std::mutex> lock(mtx);
    for(const auto p : positions) {
      bf.add_at(p);
    }
  }

private:
  BF& bf;
  std::mutex mtx;

};

template <size_t N>
void reference_1st_pass(FastaSplitter<N>& fs, KmerBuilder& kb, BFFiller1st& bff1) {
  vector<string> references;
  references.reserve(N);
  vector<vector<kmer_t>> kmer_poss;
  kmer_poss.reserve(N);
  vector<kmer_t> kmer_pos;
  size_t base_idx;
  while (true) {
    fs(references, base_idx);
    if (references.empty()) return;
    kb(references, kmer_poss);
    references.clear();
    for (auto& v: kmer_poss) {
      const size_t prev_size = kmer_pos.size();
      kmer_pos.resize(prev_size + v.size());
      std::copy(v.begin(), v.end(), kmer_pos.begin() + prev_size);
      v.clear();
    }
    kmer_poss.clear();
    bff1(kmer_pos);
    kmer_pos.clear();
  }
}

class BFFiller2nd {
public:
  BFFiller2nd(BF& _bf) : bf(_bf), next_idx(0) {}

  void operator()(std::vector<std::vector<uint64_t>>& kmer_poss, size_t idx) {
    {
      std::unique_lock<std::mutex> lock(mtx);
      cv.wait(lock, [&]{return idx == next_idx;});
      for(auto& kmer_pos : kmer_poss) {
        bf.add_to_kmer(kmer_pos, idx);
        kmer_pos.clear();
        ++idx;
      }
      next_idx = idx;
    }
    cv.notify_all();
  }

private:
  BF& bf;
  std::mutex mtx;
  size_t next_idx;
  std::condition_variable cv;
};


template <size_t N>
void reference_2nd_pass(FastaSplitter<N>& fs, KmerBuilder& kb, BFFiller2nd& bff2) {
  vector<string> references;
  references.reserve(N);
  vector<vector<kmer_t>> kmer_poss;
  kmer_poss.reserve(N);
  size_t base_idx;
  while (true) {
    fs(references, base_idx);
    if (references.empty()) return;
    kb(references, kmer_poss);
    references.clear();
    bff2(kmer_poss, base_idx);
    kmer_poss.clear();
  }
}

template <size_t N>
void read_analysis(FastqSplitter<N>& fs, ReadAnalyzer& ra, ReadOutput& ro) {
  typename FastqSplitter<N>::output_t reads;
  ReadAnalyzer::output_t associations;
  while (true) {
    fs(reads);
    if (reads.empty()) return;
    ra(reads, associations);
    reads.clear();
    ro(associations);
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
  constexpr size_t N = 128;
  {
    ref_file = gzopen(opt::fasta_path.c_str(), "r");
    kseq_t *refseq = kseq_init(ref_file);

    FastaSplitter<N> fs(refseq, &legend_ID);
    KmerBuilder kb(opt::k);
    BFFiller1st bff1(bloom);

    std::vector<std::thread> threads;
    while (static_cast<int>(threads.size()) < opt::nThreads)
      threads.emplace_back(reference_1st_pass<N>, std::ref(fs), std::ref(kb), std::ref(bff1));
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
  {
    ref_file = gzopen(opt::fasta_path.c_str(), "r");
    kseq_t *refseq = kseq_init(ref_file);

    FastaSplitter<N> fs(refseq);
    KmerBuilder kb(opt::k);
    BFFiller2nd bff2(bloom);

    std::vector<std::thread> threads;
    const size_t nThreads = std::min(static_cast<size_t>(opt::nThreads), (legend_ID.size() + N - 1) / N);
    while (threads.size() < nThreads)
      threads.emplace_back(reference_2nd_pass<N>, std::ref(fs), std::ref(kb), std::ref(bff2));
    for (auto& t: threads)
      t.join();

    kseq_destroy(refseq);
    gzclose(ref_file);
  }

  pelapsed("BF created from transcripts");

  bloom.switch_mode(2);
  pelapsed("Second switch performed");

  /****************************************************************************/

  /*** 3. Iteration over the sample *****************************************/
  constexpr size_t NS = 1 << 16;
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

    FastqSplitter<NS> fs(sseq1, sseq2, opt::min_quality, out1 != nullptr);
    ReadAnalyzer ra(bloom, legend_ID, opt::k, opt::c, opt::single);
    ReadOutput ro(out1, out2);

    std::vector<std::thread> threads;
    while (static_cast<int>(threads.size()) < opt::nThreads)
      threads.emplace_back(read_analysis<NS>, std::ref(fs), std::ref(ra), std::ref(ro));
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
