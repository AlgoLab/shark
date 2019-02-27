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
    int seq_len;
    while(fasta->size() < maxnum && (seq_len = kseq_read(seq)) >= 0) {
      fasta->push_back(make_pair(std::move(string(seq->name.s)),
				 std::move(string(seq->seq.s))));
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
