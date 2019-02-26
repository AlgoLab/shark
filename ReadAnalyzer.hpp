#ifndef READANALYZER_HPP
#define READANALYZER_HPP

#include "bloomfilter.h"
#include "kmer_utils.hpp"
#include <vector>
#include <array>

using namespace std;

class ReadAnalyzer {
private:
  BF *bf;
  vector<string> legend_ID;
  uint k;
  int c;

public:
  ReadAnalyzer(BF *_bf, vector<string> &_legend_ID, uint _k, int _c) :
    bf(_bf), legend_ID(_legend_ID), k(_k), c(_c) {}

  vector<array<string, 3>>* operator()(vector<pair<string, string>> *reads) const {
    vector<array<string, 3>> *associations = new vector<array<string, 3>>();
    for(const auto & p : *reads) {
      map<int, int> classification_id;
      string read_name = p.first;
      string read_seq = p.second;
      if(read_seq.size() >= k) {
        int _pos = 0;
        uint64_t kmer = build_kmer(read_seq, &_pos, k);
        if(kmer == (uint64_t)-1) continue;
        uint64_t rckmer = revcompl(kmer, k);
        IDView id_kmer = bf->get_index(min(kmer, rckmer));
        while (id_kmer.has_next())
          ++classification_id[id_kmer.get_next()];

        for (int pos = _pos; pos < (int)read_seq.size(); ++pos) {
          uint8_t new_char = to_int[read_seq[pos]];
          if(new_char == 0) { // Found a char different from A, C, G, T
            ++pos; // we skip this character then we build a new kmer
            kmer = build_kmer(read_seq, &pos, k);
            if(kmer == (uint64_t)-1) break;
            rckmer = revcompl(kmer, k);
            --pos; // p must point to the ending position of the kmer, it will be incremented by the for
          } else {
            --new_char; // A is 1 but it should be 0
            kmer = lsappend(kmer, new_char, k);
            rckmer = rsprepend(rckmer, reverse_char(new_char), k);
          }
          id_kmer = bf->get_index(min(kmer, rckmer));
          while (id_kmer.has_next())
            ++classification_id[id_kmer.get_next()];
        }
      }

      int max = 0;
      set<int> genes_idx;
      for(map<int,int>::iterator it=classification_id.begin(); it!=classification_id.end(); ++it) {
        if(it->second == max) {
          genes_idx.insert(it->first);
        } else if(it->second > max) {
          genes_idx.clear();
          max = it->second;
          genes_idx.insert(it->first);
        }
      }

      if(max >= c) {
        for(const auto &idx : genes_idx) {
          array<string, 3> elem;
          elem[0] = read_name;
          elem[1] = idx;
          elem[2] = read_seq;
          associations->push_back(elem);
        }
      }
    }
    delete reads;
    if(associations->size())
      return associations;
    else
      return NULL;
  }
};

#endif
