#ifndef GUARD_AnalysisConfig_h
#define GUARD_AnalysisConfig_h

#include <iostream>
#include <string>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <TString.h>
#include <TRegexp.h>
#include <map>
#include <utility>
#include <vector>
#include "ConfigReader.hh"

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
    int gen_quark;

  // PFO cuts
    int   PFO_TPCHits_max;
    float PFO_p_min;
    float PFO_p_max;
    float PFO_offset_max;
  
  // LPFO cuts
    float LPFO_p_min;
    float LPFO_p_max;

};



#endif
