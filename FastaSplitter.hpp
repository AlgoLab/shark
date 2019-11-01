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
  FastaSplitter(kseq_t * const _seq, const int _maxnum)
    : seq(_seq), maxnum(_maxnum)
  { }

  ~FastaSplitter() {
  }

  vector<pair<string, string>>* operator()(tbb::flow_control &fc) const {
    vector<pair<string, string>>* const fasta = new vector<pair<string, string>>();
    fasta->reserve(maxnum);
    int seq_len;
    while(fasta->size() < maxnum && (seq_len = kseq_read(seq)) >= 0) {
      fasta->emplace_back(seq->name.s, seq->seq.s);
    }
    if(fasta->size() > 0) return fasta;
    fc.stop();
    delete fasta;
    return NULL;
  }

private:
  kseq_t * const seq;
  const size_t maxnum;

};

#endif
