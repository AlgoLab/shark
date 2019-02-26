#ifndef READOUTPUT__HPP
#define READOUTPUT__HPP

#include <iostream>
#include <vector>
#include <array>

class ReadOutput {
 public:
  ReadOutput() { }

  void operator()(vector<array<string, 3>> *associations) const {
    if(associations) {
      for(const auto & a : *associations) {
        cout << a[0] << " " << a[1] << " " << a[2] << '\n';
      }
      delete associations;
    }
  }
  
};

#endif
