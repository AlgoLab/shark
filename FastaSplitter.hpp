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
  FastaSplitter(kseq_t * const _seq, const int _maxnum, const char _min_quality = 0)
    : seq(_seq), maxnum(_maxnum), min_quality(_min_quality)
  { }

  ~FastaSplitter() {
  }

  vector<pair<string, string>>* operator()(tbb::flow_control &fc) const {
    vector<pair<string, string>>* const fasta = new vector<pair<string, string>>();
    fasta->reserve(maxnum);
    int seq_len;
    if (min_quality == 0) {
      while(fasta->size() < maxnum && (seq_len = kseq_read(seq)) >= 0) {
        fasta->emplace_back(seq->name.s, seq->seq.s);
      }
    } else {
      const char mq = min_quality + 33;
      while(fasta->size() < maxnum && (seq_len = kseq_read(seq)) >= 0) {
        if (!seq->qual.l) {
          fasta->emplace_back(seq->name.s, seq->seq.s);
        } else {
          fasta->emplace_back(seq->name.s, mask_seq(seq->seq.s, seq->qual.s, seq->qual.l, mq));
        }
      }
    }
    if(fasta->size() > 0) return fasta;
    fc.stop();
    delete fasta;
    return NULL;
  }
private:
  kseq_t * const seq;
  const size_t maxnum;
  const char min_quality;

  static string mask_seq(string seq, const char* const qual, const size_t l, const char min_quality) {
    for (size_t i = 0; i < l; ++i) {
      if (qual[i] < min_quality) seq[i] = seq[i] - 64;
    }
    return seq;
  }
};

#endif
