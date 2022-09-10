
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
#include "../include/EventAnalyzer.hh"
#include "../include/TreeReader.hh"

using std::cout;   using std::endl;

EventAnalyzer::EventAnalyzer(TString o)
: options(o)
{
  // TEST output
    cout << "    NtupleProcessor: Created." << endl;

    patEventsAnalyzed = 0;
    entriesInNtuple   = 0;

}

bool EventAnalyzer::MapTree(TTree* tree)
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

    TreeReader reader;
    reader.InitializeMCReadTree(fChain, _mc, _branch);
    reader.InitializeJetReadTree(fChain, _jet, _branch);
    reader.InitializePFOReadTree(fChain, _pfo, _branch);

    Notify();

    return true;

}

void EventAnalyzer::Analyze()
{

  // PFO Analysis
		for (int ipfo = 0; ipfo < _pfo.pfo_n; ipfo++)
		{
      cout << _pfo.pfo_n << endl;



    }




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

bool EventAnalyzer::Select()
{ // Evaluates the class' list of event selection criteria

  /*
  Must initialize
    - Float_t MINP_CUT (= 20 GeV)
    - TString PROCESS  (= "SS")
    - TString FILE_OUT (?)
  */

  // Options
    MCParticlePair PROCESS  = ss;

    bool gpp = GenPairPicker( _mc.mc_quark_pdg[0], PROCESS );

    return gpp;

}

bool EventAnalyzer::GenPairPicker( Float_t mc_particle, MCParticlePair pair )
{
    Float_t abs_mc_particle = fabs(mc_particle);

    bool isGoodPair = (abs_mc_particle == pair) ? true : false;

    return isGoodPair;
}