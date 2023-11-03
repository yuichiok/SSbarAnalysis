/*------------------------------------------------------------------------------
   main.cpp
------------------------------------------------------------------------------*/
#define _GLIBCXX_USE_CXX11_ABI 0

// HEADERS
#include <iostream>                 // stdlib
#include <string>
#include <TFile.h>                  // ROOT class headers
#include <TString.h>
#include <TTree.h>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "timestamp.hh"
#include "NtupleProcessor.hh"

int main(int argc, char* argv[])
// int main()
{

  using std::cout; using std::endl;
  namespace po = boost::program_options;

  // C++ version check
    cout << "Current C++ version: ";
    if      (__cplusplus == 201703L) cout << "C++17\n";
    else if (__cplusplus == 201402L) cout << "C++14\n";
    else if (__cplusplus == 201103L) cout << "C++11\n";
    else if (__cplusplus == 199711L) cout << "C++98\n";
    else cout << "pre-standard C++\n";

  // Input check
    try
    {
      if ( argv[1] ){
        cout << "Processing file: " << argv[1] << endl;
      }else{
        throw 0;
      }

      // TString anconfig = "etc/SSbarAnalysisConfig_default.ini";
      TString anconfig = "etc/CCbarAnalysisConfig_default.ini";

    // Record the time main starts processing.
      std::string ts_mainBegin  = timeStamp();
      std::string fts_mainBegin = fileTimeStamp();

    // BEGIN OUTPUT
      cout << "\n\n"
              "=============================================\n"
              "===============SSbarProcessor================\n"
              "  Processing Begun: " << ts_mainBegin << "\n"
              "\n";

      NtupleProcessor ntplproc(argv[1],anconfig,"",-1);

    // CLOSING OUTPUT.
      std::string ts_mainEnd = timeStamp();

      cout << "\n"
              "   Completion time: " << ts_mainEnd <<      "\n"
              "==========SSbarProcessor - FINISHED==========\n"
              "=============================================\n" << endl;

      return 0;

    }
    catch(int exc)
    {
      std::cerr << "[Error] Input file not given.\n";
    }

}