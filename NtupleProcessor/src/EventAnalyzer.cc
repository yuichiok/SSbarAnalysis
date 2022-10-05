
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

    _hInfo     = new TTree( "info", "tree" );
    _hInfo->Branch("sum_jet_E",&sum_jet_E,"sum_jet_E/F");
    _hInfo->Branch("jet_acol",&jet_acol,"jet_acol/F");


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

void EventAnalyzer::AnalyzeGen()
{
    PFOTools mct( &_mc );
    
    PolarAngleGen(mct);

}

void EventAnalyzer::Analyze(Long64_t entry)
{
  Bool_t TreeWrite = 0;

  // MC, PFO Analysis
    PFOTools mct( &_mc );
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
  // Check all does double tagging
  // CheckTrigger = [ Valid_LPFO, Quality, Not_Gluon_K, charge_check ]
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
  




  // Try Stability and Purity Calculation here.
  Gen_Reco_Stats( mct, pfot );



  // Fill Hists can make another class called histogram extractor?
  Bool_t all_K_K = all_checks && (dEdx_pdg_check == K_K);
  PolarAngle(pfot,all_K_K);

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

    VectorTools mcvt[2];
    VectorTools jetvt[2];
    for (int i=0; i<2; i++){
      mcvt[i].SetCoordinates(_mc.mc_quark_px[i],_mc.mc_quark_py[i],_mc.mc_quark_pz[i],_mc.mc_quark_E[i]);
      jetvt[i].SetCoordinates(_jet.jet_px[i],_jet.jet_py[i],_jet.jet_pz[i],_jet.jet_E[i]);
    }

    vector<Bool_t> CutTrigger;

    switch (sel){
      case kQQ:
        CutTrigger.push_back( GenPairPicker( _mc.mc_quark_pdg[0], kSS ) );
        break;
      case kMC:
        CutTrigger.push_back( Cut_ESum( mcvt ) );
        CutTrigger.push_back( Cut_ACol( mcvt ) );
        CutTrigger.push_back( Cut_ISR ( mcvt ) );
        // cout << "Cut MC : " << Cut_ESum( mcvt ) << "," << Cut_ACol( mcvt ) << "," << Cut_ISR ( mcvt ) << "\n";
        break;
      case kReco:
        CutTrigger.push_back( Cut_ESum( jetvt ) );
        CutTrigger.push_back( Cut_ACol( jetvt ) );
        CutTrigger.push_back( Cut_ISR ( jetvt ) );
        // cout << "Cut Jet: " << Cut_ESum( jetvt ) << "," << Cut_ACol( jetvt ) << "," << Cut_ISR ( jetvt ) << "\n";
        break;

      default:
        break;
    }

  // Options

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

Bool_t EventAnalyzer::Cut_ESum ( VectorTools v[2] )
{
  Float_t SumE = v[0].v4().E() + v[1].v4().E();

  return (SumE > 220) ? true : false;
}

Bool_t EventAnalyzer::Cut_ACol ( VectorTools v[2] )
{
  Float_t cosacol = std::cos( VectorTools::GetThetaBetween( v[0].v3(), v[1].v3() ) );

  return (cosacol > 0.95) ? true : false;
}

Bool_t EventAnalyzer::Cut_ISR ( VectorTools v[2] )
{
  using namespace ROOT::Math;

	if (v[0].v4().E() < 0.5 || v[1].v4().E() < 0.5)
		return false;

  Float_t abscos[2] = {0};
  for (int i=0; i < 2; i++)
  {
    abscos[i] = fabs( std::cos( v[i].v3().Theta() ) );
  }

	Float_t ssmass  = (v[0].v4() + v[1].v4()).M();
  Float_t sinacol = std::sin( VectorTools::GetThetaBetween( v[0].v3(), v[1].v3() ) );

	Float_t Kv = 250. * sinacol / (sinacol + sqrt(1 - abscos[0] * abscos[0]) + sqrt(1 - abscos[1] * abscos[1]));
  // Float_t K[2] = {0};
  //         K[0] = jet_vec[0].v3().R() * sinacol / sqrt(1 - abscos[1] * abscos[1]);
  //         K[1] = jet_vec[1].v3().R() * sinacol / sqrt(1 - abscos[0] * abscos[0]);

  if (Kv < 35 && ssmass > 130)
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

void EventAnalyzer::Gen_Reco_Stats( PFOTools mct, PFOTools pfot )
{
  std::vector<PFO_Info> PFO_Collection;
  std::vector<PFO_Info> jet[2] = { pfot.GetJet(0), pfot.GetJet(1) };

  PFO_Collection.reserve( jet[0].size() + jet[1].size() );
  PFO_Collection.insert( PFO_Collection.begin(), jet[0].begin(), jet[0].end() );
  PFO_Collection.insert( PFO_Collection.end(), jet[1].begin(), jet[1].end() );

  std::vector<PFO_Info> PFO_K_Collection;
  for ( auto iPFO : PFO_Collection ){
    if( PFOTools::isKaon(iPFO) ) PFO_K_Collection.push_back(iPFO);
  }

  std::vector<MC_Info> Gen_K_Collection;
  for ( int istable=0; istable < _mc.mc_stable_n; istable++ ){
    if(abs(_mc.mc_stable_pdg[istable]) == 321) {
      Gen_K_Collection.push_back( mct.mc_stable[istable] );
    }
  }

  Int_t N_K_Gen = Gen_K_Collection.size();
  Int_t N_K_PFO = PFO_K_Collection.size();
  Int_t N_corr  = 0;

  Float_t cos_r = 0.02;
  std::vector<PFO_Info> PFO_K_Remain = PFO_K_Collection;

  for ( auto igen : Gen_K_Collection ){
    Float_t min_cos_diff   = 1000.0;
    Int_t   i_min_cos_diff = -1;

    Int_t counter = 0;
    for ( auto iremain : PFO_K_Remain ){
      Float_t cos_diff = igen.cos - iremain.cos;
      if( cos_diff < min_cos_diff ) {
        min_cos_diff = cos_diff;
        i_min_cos_diff = counter;
      }
      counter++;
    }

    if( min_cos_diff < cos_r ) {
      N_corr++;
      PFO_K_Remain.erase( PFO_K_Remain.begin() + i_min_cos_diff );
    }

  }

  Float_t stability = (Float_t) N_corr / (Float_t) N_K_Gen;
  Float_t purity = (Float_t) N_corr / (Float_t) N_K_PFO;


}

void EventAnalyzer::PolarAngleGen(PFOTools mct)
{
  // Gen QQbar
  for ( auto iq : mct.mc_quark ){
    _hm.h[_hm.gen_q_cos]->Fill(iq.cos);
    _hm.h[_hm.gen_q_qcos]->Fill(iq.qcos);
  }

  // Gen K
  for ( int istable=0; istable < _mc.mc_stable_n; istable++ ){
    if(abs(_mc.mc_stable_pdg[istable]) == 321) {
      _hm.h[_hm.gen_K_cos]->Fill(mct.mc_quark[istable].cos);
      _hm.h[_hm.gen_K_qcos]->Fill(mct.mc_quark[istable].qcos);
    }
  }

}

void EventAnalyzer::PolarAngle(PFOTools pfot, Bool_t b_reco)
{
  // Reco K_K
  if(b_reco){

    for ( auto iLPFO : pfot.LPFO ){
      _hm.h[_hm.reco_K_cos]->Fill( iLPFO.cos );
      _hm.h[_hm.reco_K_qcos]->Fill( iLPFO.qcos );
    }
    
  }

}

void EventAnalyzer::Jet_sum_n_acol()
{
  VectorTools jetvt[2];
  for (int i=0; i<2; i++){
    jetvt[i].SetCoordinates(_jet.jet_px[i],_jet.jet_py[i],_jet.jet_pz[i],_jet.jet_E[i]);
  }
  Float_t cosacol = std::cos( VectorTools::GetThetaBetween( jetvt[0].v3(), jetvt[1].v3() ) );

  sum_jet_E = _jet.jet_E[0] + _jet.jet_E[1];
  jet_acol  = cosacol;

  _hm.h[_hm.reco_sum_jetE]->Fill( sum_jet_E );
  _hm.h[_hm.reco_jet_sep]->Fill( cosacol );

  _hInfo->Fill();
}