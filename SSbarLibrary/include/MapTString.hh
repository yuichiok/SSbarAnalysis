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
  
}

#endif