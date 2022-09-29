
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

ClassImp(MC_QQbar)
ClassImp(TreeVariables)
ClassImp(LPFO_Info)

EventAnalyzer::EventAnalyzer(TString o)
: options(o)
{
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
    // TreeReader reader;
    // reader.InitializeMCReadTree(fChain, _mc, _branch);
    // reader.InitializeJetReadTree(fChain, _jet, _branch);
    // reader.InitializePFOReadTree(fChain, _pfo, _branch);

    Notify();

    return true;

}

void EventAnalyzer::InitWriteTree()
{

  // Initialize Write Tree
    _hfilename = TString(_fs.GetOutName_withPath());
    _hfile = new TFile( _hfilename, "RECREATE", _hfilename ) ;

    _hTree     = new TTree( "event", "tree" );
    _hTree->Branch("Event", &_eve);
    _hTree->Branch("MC", &_mc);
    _hTree->Branch("Stats_LPFO", &_stats_lpfo);
    _hTree->Branch("Data_LPFO", &_data_lpfo);

}

void EventAnalyzer::WriteFile()
{
  // Write Tree
    _hfile->Write();
    _hfile->Close();

}

void EventAnalyzer::Analyze(Long64_t entry)
{


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