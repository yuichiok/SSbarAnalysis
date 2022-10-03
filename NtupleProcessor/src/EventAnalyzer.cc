
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
    TreeReader reader;
    reader.InitializeMCReadTree(fChain, _mc, _branch);
    reader.InitializeJetReadTree(fChain, _jet, _branch);
    reader.InitializePFOReadTree(fChain, _pfo, _branch);

    Notify();

    return true;

}

void EventAnalyzer::CreateFile()
{
  // Write Tree
    _hfilename = TString(_fs.GetOutName_withPath());
    _hfile = new TFile( _hfilename, "RECREATE", _hfilename ) ;

}

void EventAnalyzer::InitWriteTree()
{

  // Initialize Write Tree

    _hTree     = new TTree( "event", "tree" );
    _hTree->Branch("Event", &_eve);
    _hTree->Branch("MC", &_mc);
    _hTree->Branch("Stats_LPFO", &_stats_lpfo);
    _hTree->Branch("Data_LPFO", &_data_lpfo);

}

void EventAnalyzer::InitHists()
{
    _hm.InitializeHists();
}

void EventAnalyzer::WriteFile()
{
  // Write Tree
    _hfile->Write();
    _hfile->Close();

}

void EventAnalyzer::Analyze(Long64_t entry)
{
  Bool_t TreeWrite = 0;

  // PFO Analysis
    PFOTools pfot( &_pfo );
    if ( !pfot.ValidPFO() ) {
      _eve.eve_valid_lpfo = 0;
      if(TreeWrite) _hTree->Fill();
      // return;
    }else{ _eve.eve_valid_lpfo = 1; }
    

  // Fill raw LPFO info
    writer.WriteLPFO_Info(pfot,&_pfo,&_stats_lpfo);
    
  ////////////////
  // Selections //
  ////////////////

    vector<Bool_t> CutTrigger;

  // Valid LPFO
    CutTrigger.push_back(_eve.eve_valid_lpfo);

  // Base Selection (mom, tpc_hit, offset)
    Bool_t LPFO_double_quality    = true;
    for ( auto iLPFO : pfot.LPFO ){
      if( !pfot.PFO_Quality_checks(iLPFO) ){
        LPFO_double_quality = false;
        break;
      }
    }
    CutTrigger.push_back(LPFO_double_quality);

  // SPFO opposite check
    Bool_t is_gluon_K[2] = {0};
    Bool_t is_there_a_gluon_K = false;
    for ( int ijet=0; ijet<2; ijet++){
      for ( auto iSPFO_K : pfot.SPFOs_K[ijet] ){
        Bool_t charge_opposite = iSPFO_K.pfo_charge * pfot.LPFO[ijet].pfo_charge < 0;
        Bool_t momentum_above  = iSPFO_K.p_mag > 10;
        if( charge_opposite && momentum_above ) is_gluon_K[ijet] = true;
      }
    }
    for ( auto ibool : is_gluon_K ){
      if( ibool ) is_there_a_gluon_K = true;
    }
    CutTrigger.push_back(!is_there_a_gluon_K);

  // dEdx dist PDG check
    enum PDGConfig { noKPi, K_K, K_Pi, Pi_Pi };
    Int_t dEdx_pdg_check = -1;
    
    if     (   pfot.isKaon(pfot.LPFO[0]) && pfot.isKaon(pfot.LPFO[1]) )  {  dEdx_pdg_check = K_K;    }
    else if(   pfot.isPion(pfot.LPFO[0]) && pfot.isPion(pfot.LPFO[1]) )  {  dEdx_pdg_check = Pi_Pi;  }
    else if( ( pfot.isKaon(pfot.LPFO[0]) && pfot.isPion(pfot.LPFO[1]) ) ||
             ( pfot.isKaon(pfot.LPFO[1]) && pfot.isPion(pfot.LPFO[0]) ) ){  dEdx_pdg_check = K_Pi;   }
    else{ dEdx_pdg_check = noKPi; }

  // charge config check
    Bool_t charge_check = false;
    switch ( dEdx_pdg_check )
    {
      case K_K:
        charge_check = pfot.is_charge_config(pfot.kOpposite);
        break;
      case K_Pi:
        charge_check = pfot.is_charge_config(pfot.kSame);

    default:
      break;
    }
    CutTrigger.push_back(charge_check);

  // Check all bools
    Bool_t all_checks = true;
    for (auto ibool : CutTrigger){
      if (!ibool) {
        all_checks = false;
        break;
      }
    }

  
  if (all_checks){

    _data_lpfo.lpfo_config = dEdx_pdg_check;
    
    for (int i=0; i<2; i++){
      _data_lpfo.lpfo_p_mag[i]         = pfot.LPFO[i].p_mag;
      _data_lpfo.lpfo_dEdx_dist_pdg[i] = pfot.LPFO[i].dEdx_dist_pdg;
      _data_lpfo.lpfo_cos[i]           = pfot.LPFO[i].cos;
      _data_lpfo.lpfo_qcos[i]          = pfot.LPFO[i].qcos;
    }

  }

  if(TreeWrite) _hTree->Fill();
  







  // Fill Hists can make another class called histogram extractor?

  for ( int imc=0; imc < 2; imc++){
    VectorTools mcqv(_mc.mc_quark_px[imc],_mc.mc_quark_py[imc],_mc.mc_quark_pz[imc],_mc.mc_quark_E[imc]);
    Float_t mc_qq_cos  = std::cos(mcqv.v3().Theta());
    Float_t mc_qq_qcos = (_mc.mc_quark_charge < 0) ? mc_qq_cos : -mc_qq_cos;

    _hm.h[_hm.gen_q_cos]->Fill(mc_qq_cos);
    _hm.h[_hm.gen_q_qcos]->Fill(mc_qq_qcos);

  }

  for ( int imc=0; imc < _mc.mc_stable_n; imc++ ){
    VectorTools mcv(_mc.mc_stable_px[imc],_mc.mc_stable_py[imc],_mc.mc_stable_pz[imc],_mc.mc_stable_E[imc]);
    Float_t mc_stable_cos  = std::cos(mcv.v3().Theta());
    Float_t mc_stable_qcos = (_mc.mc_stable_charge[imc] < 0) ? mc_stable_cos : -mc_stable_cos;

    if(abs(_mc.mc_stable_pdg[imc]) == 321) {
      _hm.h[_hm.gen_K_cos]->Fill(mc_stable_cos);
      _hm.h[_hm.gen_K_qcos]->Fill(mc_stable_qcos);
    }

  }

  for ( int ijet=0; ijet < 2; ijet++ ){

    std::vector<PFO_Info> jet = pfot.GetJet(ijet);

    for (auto ijet : jet ){

      if( pfot.isKaon(ijet) ) _hm.h[_hm.reco_K_cos]->Fill(ijet.cos);

    }

  }



  if(_eve.eve_valid_lpfo){

    for (auto iLPFO : pfot.LPFO){
      if( pfot.isKaon(iLPFO) ) _hm.h[_hm.lpfo_reco_K_mom]->Fill(iLPFO.p_mag);
    }

    for (int i=0; i<2; i++){
      if( abs(_stats_lpfo.lpfo_pdgcheat[i]) == 321 ) _hm.h[_hm.lpfo_gen_K_mom]->Fill(pfot.LPFO[i].p_mag);
    }
  }










  ClearStructs();

}

void EventAnalyzer::ClearStructs()
{
  _eve        = {};
  _stats_lpfo = {};
  _data_lpfo  = {};
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
    CutTrigger.push_back( GenPairPicker( _mc.mc_quark_pdg[0], kSS ) );
    CutTrigger.push_back( ISRPicker( 35 ) );

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