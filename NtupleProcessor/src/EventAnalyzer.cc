
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
#include "../include/EventAnalyzer.hh"
#include "../include/TreeReader.hh"
#include "../include/PFOTools.hh"

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

void EventAnalyzer::Analyze(Long64_t entry)
{

  // PFO Analysis
    cout << "event " << entry << endl;
    PFOTools pfot( &_pfo );
    if ( !pfot.ValidPFO() ) return;

    for (int ijet=0; ijet < 2; ijet++ )
    {
      LPFO[ijet] = pfot.GetSortedJet(ijet).at(0);
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

  vector<bool> boolNest;

  // Options

  // QQbar check
    MCParticlePair PROCESS  = ss;
    boolNest.push_back( GenPairPicker( _mc.mc_quark_pdg[0], PROCESS ) );

  return false;


}

bool EventAnalyzer::GenPairPicker ( Float_t mc_particle, MCParticlePair pair )
{
    Float_t abs_mc_particle = fabs(mc_particle);

    bool isGoodPair = (abs_mc_particle == pair) ? true : false;

    return isGoodPair;
}

bool EventAnalyzer::ISRPicker ( Float_t Kvcut = 25)
{
  using namespace ROOT::Math;

	if (_jet.jet_E[0] < 0.5 || _jet.jet_E[1] < 0.5)
		return false;

  // UNDER CONSTRUCTION
  // USE TLorentzVector

  XYZTVector jet0(_jet.jet_px[0],_jet.jet_py[0],_jet.jet_pz[0],_jet.jet_E[0]);
  XYZTVector jet1(_jet.jet_px[1],_jet.jet_py[1],_jet.jet_pz[1],_jet.jet_E[1]);

	Double_t ssmass = (jet0 + jet1).M();

	// TVector3 v1(jet_px[0], jet_py[0], jet_pz[0]);
	// TVector3 v2(jet_px[1], jet_py[1], jet_pz[1]);
	// VecOP vop;
	// float acol = vop.GetSinacol(v1, v2);

	// double jet0_p = sqrt(pow(jet_px[0], 2) + pow(jet_py[0], 2) + pow(jet_pz[0], 2));
	// double jet1_p = sqrt(pow(jet_px[1], 2) + pow(jet_py[1], 2) + pow(jet_pz[1], 2));

	// float costheta_j0;
	// VecOP p_j0(jet_px[0], jet_py[0], jet_pz[0]);
	// costheta_j0 = fabs(p_j0.GetCostheta());

	// float costheta_j1;
	// VecOP p_j1(jet_px[1], jet_py[1], jet_pz[1]);
	// costheta_j1 = fabs(p_j1.GetCostheta());

	// float Kv = 250. * acol / (acol + sqrt(1 - costheta_j0 * costheta_j0) + sqrt(1 - costheta_j1 * costheta_j1));
	// float K1 = jet0_p * acol / sqrt(1 - costheta_j1 * costheta_j1);
	// float K2 = jet1_p * acol / sqrt(1 - costheta_j0 * costheta_j0);

	// if (type == 1)
	// 	if (Kv < Kvcut)
	// 		return true;
	// if (type == 2)
	// 	if (Kv < Kvcut && bbmass > 130)
	// 		return true;

	return false;

}