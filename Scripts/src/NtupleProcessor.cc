
/*------------------------------------------------------------------------------
NtupleProcessor.cpp
 Created : 2022-09-02  okugawa
------------------------------------------------------------------------------*/

#include <iostream>
#include <TString.h>
#include <TFile.h> 
#include "NtupleProcessor.hh"

using std::cout;   using std::endl;
using std::string;

NtupleProcessor::NtupleProcessor(TString o, int me)
: eAnalyzer(o), tIter(eAnalyzer), options(o), maxEvents(me)
{

  // PARAM output
    cout << "  [NtupleProcessor]\n"
            "    Options:    " << options   << "\n"
            "    MaxEntries: " << me   << "\n"
    << endl;

  ntupleFile = TFile::Open(options);
  if(!ntupleFile) cout << "NtupleProcessor: ERROR: Unable to open file " << options << endl;
  TTree *ntuple   = (TTree*) ntupleFile->Get("Stats");
  if(!ntuple) cout << "NtupleProcessor: ERROR: Unable to open ttree in " << options << endl;
  ntuple->Process(&tIter, "");

}