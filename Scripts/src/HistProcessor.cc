
/*------------------------------------------------------------------------------
HistProcessor.cpp
 Created : 2022-09-02  okugawa
------------------------------------------------------------------------------*/

#include <iostream>
#include <TString.h>
#include <TFile.h> 
#include "HistProcessor.hh"

using std::cout;   using std::endl;
using std::string;

HistProcessor::HistProcessor(TString o, int me)
: eAnalyzer(o), tIter(eAnalyzer), options(o), maxEvents(me)
{

  // PARAM output
    cout << "  [HistProcessor]\n"
            "    Options:    " << options   << "\n"
            "    MaxEntries: " << me   << "\n"
    << endl;

  ntupleFile = TFile::Open(options);
  if(!ntupleFile) cout << "HistProcessor: ERROR: Unable to open file " << options << endl;
  TTree *ntuple   = (TTree*) ntupleFile->Get("event");
  if(!ntuple) cout << "HistProcessor: ERROR: Unable to open ttree in " << options << endl;
  ntuple->Process(&tIter, "");

}