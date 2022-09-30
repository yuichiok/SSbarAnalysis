
/*------------------------------------------------------------------------------
EventAnalyzer.cpp
 Created : 2022-09-05  okugawa
------------------------------------------------------------------------------*/

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include <TBranch.h>
#include <TLeaf.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>

#include "EventAnalyzer.hh"

using std::cout;   using std::endl;
typedef unsigned int Index;

ClassImp(TEvent)
ClassImp(MC_QQbar)
ClassImp(TreeVariables)
ClassImp(LPFO_Info)

EventAnalyzer::EventAnalyzer(TString o)
: options(o)
{
    _fs.SetNames(o.Data());
    patEventsAnalyzed = 0;
    entriesInNtuple   = 0;
}

Bool_t EventAnalyzer::InitReadTree(TTree* tree)
{
  // Maps TTree to class' variables.
  // TO DO: Implement check for correct mapping, return result?
  //    - Set up exception handling for negative result.

    entriesInNtuple = tree->GetEntries();

  // Set branch addresses and branch pointers
    if (!tree) return false;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

  // Read Tree
    // fChain->SetBranchAddress("Event",&_eve);
    // fChain->SetBranchAddress("MC",&_mc);
    // fChain->SetBranchAddress("Stats_LPFO",&_stats_lpfo);
    fChain->SetBranchAddress("Data_LPFO",&_data_lpfo);

    Notify();

    return true;

}

void EventAnalyzer::Analyze(Long64_t entry)
{

  if(_data_lpfo->lpfo_config) cout << "HERE!!!!" << endl;

}

Bool_t EventAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}