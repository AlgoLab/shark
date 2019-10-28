#ifndef READOUTPUT__HPP
#define READOUTPUT__HPP

#include <iostream>
#include <vector>
#include <array>

class ReadOutput {
 private:
  static void mask_seq(string& s) {
    for (size_t i = 0; i < s.length(); ++i) {
      if (s[i] < 64) {
        s[i] += 64;
        if (65 <= s[i] && s[i] <= 90) { // upper case
          s[i] += 32;
        }
      } else {
        if (97 <= s[i] && s[i] <= 122) { // lower case
          s[i] -= 32;
        }
      }
    }
  }
 public:
  ReadOutput() { }

  void operator()(vector<array<string, 4>> *associations) const {
    if(associations) {
      for(auto & a : *associations) {
        mask_seq(a[2]);
        printf("%s %s %s %s\n", a[0].c_str(), a[1].c_str(), a[2].c_str(), a[3].c_str());
        //cout << a[0] << " " << a[1] << " " << a[2] << '\n';
      }
      delete associations;
    }
  }
};

#endif
