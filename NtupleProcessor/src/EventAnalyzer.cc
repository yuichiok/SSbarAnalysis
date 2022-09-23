
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
#include "TreeReader.hh"
#include "PFOTools.hh"
#include "VectorTools.hh"

using std::cout;   using std::endl;
typedef unsigned int Index;

template<typename T>
void pop_front(std::vector<T>& vec)
{
    assert(!vec.empty());
    vec.front() = std::move(vec.back());
    vec.pop_back();
}

EventAnalyzer::EventAnalyzer(TString o)
: options(o)
{
    patEventsAnalyzed = 0;
    entriesInNtuple   = 0;
}

Bool_t EventAnalyzer::MapTree(TTree* tree)
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
    PFOTools pfot( &_pfo );
    if ( !pfot.ValidPFO() ) return;

    vector<PFO_Info> sorted_jets[2] = {pfot.GetSortedJet(0), pfot.GetSortedJet(1)};

    for (int ijet=0; ijet < 2; ijet++ )
    {
      LPFO[ijet]  = sorted_jets[ijet].at(0);
      SPFOs[ijet] = sorted_jets[ijet];
      // SPFOs[ijet].erase(SPFOs[ijet].begin());
      pop_front(SPFOs[ijet]); // faster algorithm wise?
    }

  // Selection
    Select( kLPFO );
  


}

Bool_t EventAnalyzer::Select(Selector sel)
{ // Evaluates the class' list of event selection criteria

  /*
  Must initialize
    - Float_t MINP_CUT (= 20 GeV)
    - TString PROCESS  (= "SS")
    - TString FILE_OUT (?)
  */

  vector<Bool_t> CutTrigger;

  // Options

  switch (sel) {

    case kMC:
      // QQbar checks
        CutTrigger.push_back( GenPairPicker( _mc.mc_quark_pdg[0], kSS ) );
        CutTrigger.push_back( ISRPicker( 35 ) );
        break;

    case kLPFO:
      // LPFO checks
        CutTrigger.push_back( is_charge_config( kOpposite ) );    // Charge opposite check
        CutTrigger.push_back( PFO_Quality_checks() );             // Double Tagger
        break;

    default:
      break;

  }

  for (auto trigger : CutTrigger ){
    if (!trigger) { return false; }
  }

  return true;


}

Bool_t EventAnalyzer::GenPairPicker ( Float_t mc_particle, MCParticlePair pair )
{
    Float_t abs_mc_particle = fabs(mc_particle);

    Bool_t isGoodPair = (abs_mc_particle == pair) ? true : false;

    return isGoodPair;
}

Bool_t EventAnalyzer::ISRPicker ( Float_t Kvcut = 25)
{
  using namespace ROOT::Math;

	if (_jet.jet_E[0] < 0.5 || _jet.jet_E[1] < 0.5)
		return false;

  VectorTools    jet_vec[2];
  Float_t     jet_abscos[2] = {0};
  for (int ijet=0; ijet < 2; ijet++)
  {
    jet_vec[ijet].SetCoordinates(_jet.jet_px[ijet],_jet.jet_py[ijet],_jet.jet_pz[ijet],_jet.jet_E[ijet]);
    jet_abscos[ijet] = fabs( std::cos( jet_vec[ijet].v3().Theta() ) );
  }

	Float_t ssmass = (jet_vec[0].v4() + jet_vec[1].v4()).M();
  Float_t acol   = VectorTools::GetSinacol( jet_vec[0].v3(), jet_vec[1].v3() );

	Float_t Kv = 250. * acol / (acol + sqrt(1 - jet_abscos[0] * jet_abscos[0]) + sqrt(1 - jet_abscos[1] * jet_abscos[1]));
  // Float_t K[2] = {0};
  //         K[0] = jet_vec[0].v3().R() * acol / sqrt(1 - jet_abscos[1] * jet_abscos[1]);
  //         K[1] = jet_vec[1].v3().R() * acol / sqrt(1 - jet_abscos[0] * jet_abscos[0]);

  if (Kv < Kvcut && ssmass > 130)
    return true;

	return false;

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

Bool_t EventAnalyzer::is_charge_config( ChargeConfig cc )
{

  Int_t charge_config = LPFO[0].charge * LPFO[1].charge;

  switch (cc)
  {
    case kSame:
      if ( charge_config > 0 ) return true;
      break;
    
    case kOpposite:
      if ( charge_config < 0 ) return true;
      break;

    default:
      return false;
      break;
  }

  return false;

}

Bool_t EventAnalyzer::PFO_Quality_checks()
{
  vector<Bool_t> CutTrigger;

  for (auto iLPFO : LPFO ){
    CutTrigger.push_back( is_momentum( iLPFO, 20.0, 60.0 ) );     // MIN/MAX momentum check
    CutTrigger.push_back( is_tpc_hits( iLPFO, 210 ) );            // Number of TPC hit check
    CutTrigger.push_back( is_offset_small( iLPFO, 1.0 ) );        // Offset distance check
  }
  
  for (auto trigger : CutTrigger ){
    if (!trigger) { return false; }
  }

  return true;

}

Bool_t EventAnalyzer::is_momentum( PFO_Info iPFO, Float_t MINP_CUT, Float_t MAXP_CUT )
{
  return ( MINP_CUT < iPFO.p_mag && iPFO.p_mag < MAXP_CUT);
}

Bool_t EventAnalyzer::is_tpc_hits( PFO_Info iPFO, Int_t MIN_TPC_HITS )
{
  return ( MIN_TPC_HITS < iPFO.tpc_hits );
}

Bool_t EventAnalyzer::is_offset_small( PFO_Info iPFO, Int_t MAX_OFFSET )
{
  return ( iPFO.pv < MAX_OFFSET );
}