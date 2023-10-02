#ifndef GUARD_MapTString_h
#define GUARD_MapTString_h

/*------------------------------------------------------------------------------
  MapTString
 Created : 2023-07-19  okugawa
 Main class of MapTString program.
------------------------------------------------------------------------------*/
#include <TString.h>
#include <unordered_map>

namespace std {

  template <>
  struct hash<TString> {
    std::size_t operator()(const TString & k) const {
      return std::hash<std::string>()(k.Data());
    };
  };
  // A hash function used to hash a pair of any kind
  struct hash_pair {
      // template <class TString, class TString>
      size_t operator()(const pair<TString, TString>& p) const
      {
          auto hash1 = hash<TString>{}(p.first);
          auto hash2 = hash<TString>{}(p.second);
  
          if (hash1 != hash2) {
              return hash1 ^ hash2;             
          }
          
          // If hash1 == hash2, their XOR is zero.
            return hash1;
      }
  };
}

#endif