#include "bloomfilter.h"

IDView &IDView::operator=(const IDView &rhs) {
  if (this == &rhs)
    return *this;
  _b = rhs._b;
  _e = rhs._e;
  _p = rhs._p;
  _bf = rhs._bf;
  return *this;
}

int_vector<64>::value_type IDView::get_next() { return _bf->_index_kmer[_p++]; }

void IDView::clear() {
  _b = 0;
  _e = 0;
  _p = 0;
  _bf = nullptr;
}
