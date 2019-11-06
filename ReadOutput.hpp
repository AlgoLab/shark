#ifndef READOUTPUT__HPP
#define READOUTPUT__HPP

#include <iostream>
#include <vector>
#include "common.hpp"

class ReadOutput {
public:
  ReadOutput(FILE* const _out1 = nullptr, FILE* const _out2 = nullptr)
    : out1(_out1), out2(_out2)
  { }

  void operator()(std::vector<assoc_t> *associations) const {
    if(associations) {
      string previd = "";
      for(const auto & a : *associations) {
        const sharseq_t& s1 = a.second.first;
        const sharseq_t& s2 = a.second.second;
        printf("%s %s\n", s1.id.c_str(), a.first.c_str());
        if (out1 != nullptr && previd != s1.id)
          fprintf(out1, "@%s\n%s\n+\n%s\n", s1.id.c_str(), s1.seq.c_str(), s1.qual.c_str());
        if (out2 != nullptr && previd != s1.id)
          fprintf(out2, "@%s\n%s\n+\n%s\n", s2.id.c_str(), s2.seq.c_str(), s2.qual.c_str());
        previd = std::move(s1.id);
      }
      delete associations;
    }
  }

private:
  FILE* const out1;
  FILE* const out2;
};

#endif
