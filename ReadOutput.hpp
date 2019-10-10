#ifndef READOUTPUT__HPP
#define READOUTPUT__HPP

#include <iostream>
#include <vector>
#include <array>

class ReadOutput {
 public:
  ReadOutput() { }

  void operator()(vector<array<string, 4>> *associations) const {
    if(associations) {
      for(const auto & a : *associations) {
        printf("%s %s %s %s\n", a[0].c_str(), a[1].c_str(), a[2].c_str(), a[3].c_str());
        //cout << a[0] << " " << a[1] << " " << a[2] << '\n';
      }
      delete associations;
    }
  }
  
};

#endif
