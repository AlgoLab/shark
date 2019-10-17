#ifndef _KMER_UTILS_HPP
#define _KMER_UTILS_HPP

using namespace std;

static const uint8_t to_int[128] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 10
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 20
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 30
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 40
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 50
                                    0, 0, 0, 0, 0, 1, 0, 2, 0, 0, // 60
                                    0, 3, 0, 0, 0, 0, 0, 0, 0, 0, // 70
                                    0, 0, 0, 0, 4, 0, 0, 0, 0, 0, // 80
                                    0, 0, 0, 0, 0, 0, 0, 1, 0, 2, // 90
                                    0, 0, 0, 3, 0, 0, 0, 0, 0, 0, // 100
                                    0, 0, 0, 0, 0, 0, 4, 0, 0, 0, // 110
                                    0, 0, 0, 0, 0, 0, 0, 0};      // 120

inline uint8_t reverse_char(const uint8_t c) {
  return ((~c) & 3);
}

uint64_t revcompl(uint64_t kmer, const uint8_t k) {
  uint64_t rckmer = 0;
  kmer = ~kmer;
  for(uint i = 0; i < k; ++i) {
    rckmer = (rckmer << 2) | (kmer & 3);
    kmer >>= 2;
  }
  return rckmer;
}

int64_t build_kmer(const string &seq, int &p, const uint8_t k) {
  for(int _p = p; _p < (int)seq.size() && _p < p+k; ++_p) {
    if(to_int[seq[_p]] == 0) p = _p + 1;
  }
  if(p+k > (int)seq.size()) {
    p = seq.size();
    return -1;
  }
  // we can build the k-mer starting at p
  uint64_t kmer = 0;
  for(int end = p + k; p < end; ++p) {
    kmer = (kmer << 2) | (to_int[seq[p]] - 1);
  }
  return kmer;
}

inline uint64_t lsappend(const uint64_t kmer, const uint64_t c, const uint64_t k) { // left shift and append
  return ((kmer << 2) | c) & ((1UL << 2*k)-1);
}

inline uint64_t rsprepend(const uint64_t kmer, const uint64_t c, const uint64_t k) { // right shift and prepend
  return (kmer >> 2) | (c << (2*k - 2));
}

#endif
