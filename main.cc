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
#include "SSbarLibrary/timestamp.hh"
#include "NtupleProcessor/include/NtupleProcessor.hh"

using namespace std;

// int main(int argc, char* argv[])
int main()
{
  // Record the time main starts processing.
  TString ts_mainBegin  = timeStamp();
  TString fts_mainBegin = fileTimeStamp();

  // BEGIN OUTPUT
  cout << "\n\n"
          "============================================\n"
          "==========SSbarHistogramExtractor===========\n"
          "  Processing Begun: " << ts_mainBegin << "\n"
          "\n";

  NtupleProcessor ntplproc("",-1);

  // CLOSING OUTPUT.
    TString ts_mainEnd = timeStamp();

    cout << "\n"
            "   Completion time: " << ts_mainEnd <<      "\n"
            "=====SSbarHistogramExtractor - FINISHED=====\n"
            "============================================\n" << endl;

    return 0;
}