#ifndef GUARD_AnalysisConfig_h
#define GUARD_AnalysisConfig_h

#include "MapTString.hh"
#include "ConfigReader.hh"

#include <iostream>
#include <string>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <TString.h>
#include <TRegexp.h>
#include <utility>
#include <vector>

// class AnalysisConfig : public ConfigReader
class AnalysisConfig
{
  public:
    AnalysisConfig();
   ~AnalysisConfig(){}

    ConfigReader cr;
    void SetConfig(TString fnc = "etc/SSbarAnalysisConfig_default.ini");

  // Config variables
  // Gen cuts
    // int gen_quark;
    std::vector<int> gen_quarks;

  // PFO cuts
    int   PFO_TPCHits_min;
    float PFO_p_min;
    float PFO_p_max;
    float PFO_offset_max;
  
  // LPFO cuts
    float LPFO_p_min;
    float LPFO_p_max;

  protected:
    template <typename T>
    void getListFromString(std::string& str, std::vector<T>& list)
    { // Simple function that extracts numbers from a string.
      // Feeds string into a stringstream and, while there is still something to
      //   read out, ouputs the entry from the stream into an variable and
      //   pushes the variable into the output vector.
        list.clear();
        std::stringstream strm(str);
        while(true) {
            T n;
            strm >> n;
            if(!strm) break;
            list.push_back(n);
        }
    }

};



#endif
