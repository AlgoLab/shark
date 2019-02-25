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
  FastaSplitter(kseq_t * _seq, const int _maxnum) : maxnum(_maxnum) {
    seq = _seq;
  }

  ~FastaSplitter() {
  }

  vector<pair<string, string>>* operator()(tbb::flow_control &fc) const {
    vector<pair<string, string>>* fasta = new vector<pair<string, string>>();
    int file_line;
    // kseq_t *seq;
    while(fasta->size() < maxnum && (file_line = kseq_read(seq)) >= 0) {
      fasta->push_back(make_pair(string(seq->name.s), string(seq->seq.s)));
    }
    if(fasta->size() > 0) return fasta;
    fc.stop();
    return NULL;
  }
private:
  mutable kseq_t *seq;
  size_t maxnum;
};

#endif
