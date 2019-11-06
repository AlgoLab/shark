#ifndef SHARK_COMMON_HPP
#define SHARK_COMMON_HPP

#include <string>
#include <utility>

using std::string;

struct sharseq_t {
  string id, seq, qual;
};

typedef std::pair<string, std::pair<sharseq_t, sharseq_t>> elem_t;

typedef std::pair<string, std::pair<sharseq_t, sharseq_t>> assoc_t;

#endif
