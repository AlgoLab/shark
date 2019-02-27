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

inline uint8_t reverse_char(const uint8_t &c) {
  return ((~c) & 3);
}

uint64_t revcompl(const uint64_t &kmer, const uint8_t &k) {
  uint64_t rckmer = 0;
  for(uint8_t i = 0; i<=2*k - 2; ++(++i)) {
    rckmer = (rckmer << 2) | (~(kmer >> i) & 3);
  }
  return rckmer;
}

int64_t build_kmer(const string &seq, int *p, const uint8_t &k) {
  bool found_kmer = false;
  uint64_t kmer = 0;
  uint curr_klen = 0;
  while(!found_kmer && *p<(int)seq.size()) {
    uint8_t new_char = to_int[seq[*p]];
    if(new_char == 0) {
      // Found a char different from A, C, G, T
      kmer = 0;
      curr_klen = 0;
    } else {
      kmer = (kmer << 2) | (new_char-1); // -1 since N is 0 in the table, A is 1 but it should be 0
      ++curr_klen;
    }
    ++(*p);
    found_kmer = curr_klen == k;
  }
  if(!found_kmer)
    return -1;
  else
    return kmer;
}

inline uint64_t lsappend(const uint64_t &kmer, const uint8_t &c, const uint8_t &k) { // left shift and append
  return ((kmer << 2) | c) & (((uint64_t)1 << 2*k)-1);
}

inline uint64_t rsprepend(const uint64_t &kmer, const uint8_t &c, const uint8_t &k) { // right shift and prepend
  return (kmer >> 2) | ((uint64_t)c << (2*k - 2));
}

#endif
