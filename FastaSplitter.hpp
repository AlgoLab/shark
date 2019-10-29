#ifndef FASTA_SPLITTER_HPP
#define FASTA_SPLITTER_HPP

#include "kseq.h"
#include <zlib.h>
#include <string>
#include <vector>
#include <memory>

#include <tbb/pipeline.h>

using namespace std;

class FastaSplitter {
public:
  FastaSplitter(kseq_t * const _seq1, kseq_t * const _seq2, const int _maxnum, const char _min_quality = 0)
    : seq1(_seq1), seq2(_seq2), maxnum(_maxnum), min_quality(_min_quality)
  { }

  ~FastaSplitter() {
  }

  vector<pair<string, string>>* operator()(tbb::flow_control &fc) const {
    vector<pair<string, string>>* const fasta = new vector<pair<string, string>>();
    fasta->reserve(maxnum);
    int seq_len1, seq_len2;
    if (min_quality == 0) {
      if (seq2 == nullptr) {
        while(fasta->size() < maxnum && (seq_len1 = kseq_read(seq1)) >= 0) {
          fasta->emplace_back(seq1->name.s, seq1->seq.s);
        }
      } else {
        while(fasta->size() < maxnum && (seq_len1 = kseq_read(seq1)) >= 0 && (seq_len2 = kseq_read(seq2)) >= 0) {
          fasta->emplace_back(seq1->name.s, string(seq1->seq.s) + string(seq2->seq.s));
        }
      }
    } else {
      const char mq = min_quality + 33;
      if (seq2 == nullptr) {
        while(fasta->size() < maxnum && (seq_len1 = kseq_read(seq1)) >= 0) {
          if (!seq1->qual.l) {
            fasta->emplace_back(seq1->name.s, seq1->seq.s);
          } else {
            fasta->emplace_back(seq1->name.s, mask_seq(seq1->seq.s, seq1->qual.s, seq1->qual.l, mq));
          }
        }
      } else {
        while(fasta->size() < maxnum && (seq_len1 = kseq_read(seq1)) >= 0 && (seq_len2 = kseq_read(seq2)) >= 0) {
          if (!seq1->qual.l || !seq2->qual.l) {
            fasta->emplace_back(seq1->name.s, string(seq1->seq.s) + string(seq2->seq.s));
          } else {
            fasta->emplace_back(seq1->name.s,
                                mask_seq(string(seq1->seq.s) + string(seq2->seq.s),
                                         string(seq1->qual.s) + string(seq2->qual.s),
                                         mq
                                         )
                                );
          }
        }
      }
    }
    if(fasta->size() > 0) return fasta;
    fc.stop();
    delete fasta;
    return NULL;
  }
private:
  kseq_t * const seq1;
  kseq_t * const seq2;
  const size_t maxnum;
  const char min_quality;

  static string mask_seq(string seq, const char* const qual, const size_t l, const char min_quality) {
    for (size_t i = 0; i < l; ++i) {
      if (qual[i] < min_quality) seq[i] = seq[i] - 64;
    }
    return seq;
  }

  static string mask_seq(string seq, const string& qual, const char min_quality) {
    return mask_seq(seq, qual.c_str(), qual.length(), min_quality);
  }
};

#endif
