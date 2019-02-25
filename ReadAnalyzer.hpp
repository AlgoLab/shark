#ifndef READANALYZER_HPP
#define READANALYZER_HPP

#include "bloomfilter.h"
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
        string kmer(read_seq, 0, k);
        transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
        IDView id_kmer = bf->get_index(kmer);
        while (id_kmer.has_next()) {
          int x = id_kmer.get_next();
          ++classification_id[x];
        }
        for (uint p = k; p < read_seq.size(); ++p) {
          char c = toupper(read_seq[p]);
          kmer.erase(0, 1);
          kmer += c;
          id_kmer = bf->get_index(kmer);
          while (id_kmer.has_next()) {
            int x = id_kmer.get_next();
            ++classification_id[x];
          }
        }
      }
      int max = std::max(max_element(begin(classification_id),
                                     end(classification_id),
                                     [](const pair<int, int> &a,
                                        const pair<int, int> &b) {
                                       return a.second < b.second;
                                     })
                         ->second,
                         (int)3);
      if(max >= c) {
        unordered_set<string> output_names;
        for (auto it_class = classification_id.cbegin();
             it_class != classification_id.cend(); ++it_class) {
          if (it_class->second == max) {
            string associated_name = legend_ID[it_class->first];
            if(output_names.find(associated_name) == output_names.end()) {
              output_names.insert(associated_name);
              array<string, 3> elem;
              elem[0] = read_name;
              elem[1] = associated_name;
              elem[2] = read_seq;
              associations->push_back(elem);
            }
          }
        }
      }
    }
    if(associations->size())
      return associations;
    else
      return NULL;
  }
};
      
#endif
