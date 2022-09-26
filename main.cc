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
#include "timestamp.hh"
#include "NtupleProcessor.hh"

using namespace std;

// int main(int argc, char* argv[])
int main()
{

  if (__cplusplus == 201703L) std::cout << "C++17\n";
  else if (__cplusplus == 201402L) std::cout << "C++14\n";
  else if (__cplusplus == 201103L) std::cout << "C++11\n";
  else if (__cplusplus == 199711L) std::cout << "C++98\n";
  else std::cout << "pre-standard C++\n";


  // Record the time main starts processing.
  string ts_mainBegin  = timeStamp();
  string fts_mainBegin = fileTimeStamp();

  // BEGIN OUTPUT
  cout << "\n\n"
          "=============================================\n"
          "===============SSbarProcessor================\n"
          "  Processing Begun: " << ts_mainBegin << "\n"
          "\n";

  // NtupleProcessor ntplproc(argv[1],-1);
  NtupleProcessor ntplproc("data/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n002.d_dstm_15162_000.root",-1);

  // CLOSING OUTPUT.
    string ts_mainEnd = timeStamp();

    cout << "\n"
            "   Completion time: " << ts_mainEnd <<      "\n"
            "==========SSbarProcessor - FINISHED==========\n"
            "=============================================\n" << endl;

    return 0;
}