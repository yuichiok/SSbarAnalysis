
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

  // TString dummy_label = "rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n002.d_dstm_15162_000";
  // TString filename = "data/" + dummy_label + ".root";


  ntupleFile = TFile::Open(o);
  if(!ntupleFile) cout << "NtupleProcessor: ERROR: Unable to open file " << o << endl;
  TTree *ntuple   = (TTree*) ntupleFile->Get("Stats");
  if(!ntuple) cout << "NtupleProcessor: ERROR: Unable to open ttree in " << o << endl;
  ntuple->Process(&tIter, "");



  // Takes input options and processes appropriately
  //   Options that can be specified:
  //     - Operating points
  //     - Differential variables (default: JetPT)
  //     - Dataset ntuple file to run on (default: all)
  //     - Desired number of entry to run on per file and entries to skip.
  //     - Modifications to the analysis (systematic analyses, etc)
  // For each dataset that needs to be run over...
  //   Opens the appropriate file and tree
  //   Creates a TreeIterator, EventHandler, and the desired HistogramMakers
  //   Runs the TreeIterator.

}