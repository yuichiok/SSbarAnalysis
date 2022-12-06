#ifndef GUARD_ConfigReader_h
#define GUARD_ConfigReader_h

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <TString.h>

using std::cout; using std::endl;

class ConfigReader {
  public:
    ConfigReader(){}
    virtual ~ConfigReader(){}
    
    void SetConfigReader(TString fnc)
    {
      TString fn_config;
      boost::property_tree::ptree pt;

      fn_config = fnc;
      // Open and read in config file
        boost::property_tree::ini_parser::read_ini(fn_config.Data(), pt);
    }

  private:

};



#endif
