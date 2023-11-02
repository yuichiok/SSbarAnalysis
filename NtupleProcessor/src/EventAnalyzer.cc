
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

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

EventAnalyzer::EventAnalyzer(TString input, TString fnac, TString o)
: options(o), _config(fnac)
{
    _fs.SetNames(input.Data());
    _anCfg.SetConfig(_config);
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
    reader.InitializeVTXReadTree(fChain,_vtx, _branch);

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
    _hTree     = new TTree( "data", "tree" );
    writer.InitializeDataTree(_hTree,_data);
}

void EventAnalyzer::InitHists()
{
    _hm.InitializeHists();
}

void EventAnalyzer::WriteFile()
{
  // Write Tree
    _hm.WriteLists(_hfile);
    _hTree->Write();
    // _hfile->Write();
    _hfile->Close();

}

void EventAnalyzer::AnalyzeReco(Long64_t entry)
{
  // MC, PFO Analysis
  PFOTools mct( &_mc, _config );
  PFOTools pfot( &_mc, &_pfo, _config );
  
  // cout << "===== " << entry + 1 << " =====" << endl;

    if ( !pfot.ValidPFO() ) {
      _eve.eve_valid_lpfo = 0;
    }else{ _eve.eve_valid_lpfo = 1; }

    // Fill raw LPFO info
    writer.WriteLPFO_Info(pfot,&_pfo,&_stats_lpfo);
    for (int i=0; i<2; i++){
      _data_lpfo.lpfo_p_mag[i]         = pfot.KLPFO[i].p_mag;
      _data_lpfo.lpfo_dEdx_dist_pdg[i] = pfot.KLPFO[i].dEdx_dist_pdg;
      _data_lpfo.lpfo_cos[i]           = pfot.KLPFO[i].cos;
      _data_lpfo.lpfo_qcos[i]          = pfot.KLPFO[i].qcos;
    }

    for (int i=0; i<2; i++){
      _data.LPFO_cos[i]  = pfot.KLPFO[i].cos;
      _data.LPFO_qcos[i] = pfot.KLPFO[i].qcos;
    }

  vector<Bool_t> CutTrigger[3];

  // Valid LPFO
    CutTrigger[kKaon].push_back(_eve.eve_valid_lpfo);
    CutTrigger[kPion].push_back(_eve.eve_valid_lpfo);

  // Base Selection (mom, tpc_hit, offset)
    Bool_t LPFO_double_quality[3]    = {true,true,true};
    for ( int i=0; i<2; i++ ){
      if( !pfot.LPFO_Quality_checks(pfot.KLPFO[i]) )  LPFO_double_quality[kKaon] = false;
      if( !pfot.LPFO_Quality_checks(pfot.PiLPFO[i]) ) LPFO_double_quality[kPion] = false;
    }
    CutTrigger[kKaon].push_back(LPFO_double_quality[kKaon]);
    CutTrigger[kPion].push_back(LPFO_double_quality[kPion]);

  // SPFO opposite check
    Bool_t is_gluon[3][2] = {0};
    Bool_t is_there_a_gluon[3] = {false};
    for ( int ijet=0; ijet<2; ijet++){

      Bool_t check_KSLPFO_PiLPFO = false;

      for ( auto iSPFO_K : pfot.SPFOs_K[ijet] ){
        Bool_t charge_opposite = iSPFO_K.pfo_charge * pfot.KLPFO[ijet].pfo_charge < 0;
        Bool_t charge_opposite_KSLPFO_PiLPFO = iSPFO_K.pfo_charge * pfot.PiLPFO[ijet].pfo_charge < 0;
        Bool_t momentum_above  = iSPFO_K.p_mag > 10;
        if( charge_opposite && momentum_above ) is_gluon[kKaon][ijet] = true;
        if( charge_opposite_KSLPFO_PiLPFO && momentum_above ) check_KSLPFO_PiLPFO = true;
        if( charge_opposite_KSLPFO_PiLPFO && iSPFO_K.p_mag > 15 ) check_KSLPFO_PiLPFO = true;
      }

      for ( auto iSPFO_Pi : pfot.SPFOs_Pi[ijet] ){
        Bool_t charge_opposite = iSPFO_Pi.pfo_charge * pfot.PiLPFO[ijet].pfo_charge < 0;
        Bool_t momentum_above  = iSPFO_Pi.p_mag > 10;
        if( (charge_opposite && momentum_above) || check_KSLPFO_PiLPFO ) is_gluon[kPion][ijet] = true;
        if( charge_opposite && momentum_above ) is_gluon[kPion][ijet] = true;
      }

    }

    for ( int i=0; i<2; i++ ){
      if( is_gluon[kKaon][i] ) is_there_a_gluon[kKaon] = true;
      if( is_gluon[kPion][i] ) is_there_a_gluon[kPion] = true;
    }
    CutTrigger[kKaon].push_back( !is_there_a_gluon[kKaon] );
    CutTrigger[kPion].push_back( !is_there_a_gluon[kPion] );

  // dEdx dist PDG check
    PDGConfig dEdx_pdg_match[3] = {noKPi, noKPi, noKPi};
    
    // if ( pfot.isKaon(pfot.KLPFO[0])  && pfot.isKaon(pfot.KLPFO[1]) && pfot.is_cheatNoOthers(pfot.KLPFO[0]) && pfot.is_cheatNoOthers(pfot.KLPFO[1])  ) { dEdx_pdg_match[kKaon] = K_K; }
    // if ( pfot.isPion(pfot.PiLPFO[0]) && pfot.isPion(pfot.PiLPFO[1]) && pfot.is_cheatNoOthers(pfot.PiLPFO[0]) && pfot.is_cheatNoOthers(pfot.PiLPFO[1]) ) { dEdx_pdg_match[kPion] = Pi_Pi; }

    if ( pfot.isKaon(pfot.KLPFO[0])  && pfot.isKaon(pfot.KLPFO[1])  ) { dEdx_pdg_match[kKaon] = K_K; }
    if ( pfot.isPion(pfot.PiLPFO[0]) && pfot.isPion(pfot.PiLPFO[1]) ) { dEdx_pdg_match[kPion] = Pi_Pi; }


  // charge config check
    Bool_t charge_check[3] = {false};

    if( dEdx_pdg_match[kKaon] == K_K ){
      charge_check[kKaon] = pfot.is_charge_config(pfot.kOpposite,pfot.KLPFO[0].pfo_charge,pfot.KLPFO[1].pfo_charge);
    }
    if( dEdx_pdg_match[kPion] == Pi_Pi ){
      charge_check[kPion] = pfot.is_charge_config(pfot.kOpposite,pfot.PiLPFO[0].pfo_charge,pfot.PiLPFO[1].pfo_charge);
    }

    CutTrigger[kKaon].push_back(charge_check[kKaon]);
    CutTrigger[kPion].push_back(charge_check[kPion]);


  Double_t ctag[2] = {0,0};
  for(Int_t ijet=0; ijet<2; ijet++){
    ctag[ijet] = _jet.jet_ctag[ijet];
  }
  CCbarAnalysis(pfot,CutTrigger,dEdx_pdg_match, ctag, entry + 1);
  
  // _hTree->Fill();
  ClearStructs();

}

void EventAnalyzer::ClearStructs()
{
  _eve        = {};
  _stats_lpfo = {};
  _data_lpfo  = {};

  _data       = {};
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
        CutTrigger.push_back( GenPairPicker( _mc.mc_quark_pdg[0], _anCfg.gen_quarks ) );
        break;
      case kMC:
        CutTrigger.push_back( Cut_ESum( mcvt ) );
        CutTrigger.push_back( Cut_ACol( mcvt ) );
        CutTrigger.push_back( Cut_ISR ( mcvt ) );
        break;
      case kReco:
        CutTrigger.push_back( Cut_ESum( jetvt ) );
        CutTrigger.push_back( Cut_ACol( jetvt ) );
        CutTrigger.push_back( Cut_ISR ( jetvt ) );
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

Bool_t EventAnalyzer::GenPairPicker ( Float_t mc_particle, std::vector<int> input_gen )
{
  for ( auto igen : input_gen ){
    if( fabs(mc_particle) == igen ) return true;
  }
  
  return false;
}

Bool_t EventAnalyzer::Cut_ESum ( VectorTools v[2] )
{
  Float_t SumE = v[0].v4().E() + v[1].v4().E();

  return (SumE > 220);
}

Bool_t EventAnalyzer::Cut_ACol ( VectorTools v[2] )
{
  Float_t cosacol = std::cos( VectorTools::GetThetaBetween( v[0].v3(), v[1].v3() ) );

  return (cosacol > 0.95);
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

  return (Kv < 35 && ssmass > 130);

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

void EventAnalyzer::CCbarAnalysis(PFOTools pfot, vector<Bool_t> cuts[3], PDGConfig double_tag[3], Double_t ctag[2], Double_t event)
{  
  std::vector<PFO_Info> PFO_Collection = pfot.Get_Valid_PFOs();

  Int_t nprong[2][2] = {{0,0},{0,0}};

  for ( long unsigned int i=0; i < PFO_Collection.size(); i++ )
  {
    PFO_Info ipfo = PFO_Collection.at(i);
    Int_t imatch = ipfo.pfo_match;

    if(ipfo.pfo_vtx==0){
      nprong[imatch][0]++;
    }else if(ipfo.pfo_vtx==1){
      nprong[imatch][1]++;
    }
  }

  Bool_t LPFO_checks[3] = {true,true,true};
  
  for (int i=0; i<3; i++ ){

    if ( cuts[i].empty() ) continue;

    for (int j=0; j< cuts[i].size()-1; j++){

      if (!cuts[i].at(j)){
        LPFO_checks[i] = false;
        break;
      }

    }
  }

  Bool_t sign_check[3] = {false, false, false};
  for (int i=0; i<3; i++ ){
    if ( cuts[i].empty() ) continue;
    sign_check[i] = cuts[i].back();
  }

  Double_t maxCtag = (ctag[0] > ctag[1]) ? ctag[0] : ctag[1];
  for(Int_t ijet=0; ijet<2; ijet++) {
    if(LPFO_checks[kKaon]) {
      _hm.h_general[static_cast<int>(HistManager::h_general::cuts)]->Fill(0);
      if(pfot.KLPFO[ijet].p_mag > 15) {
        _hm.h_general[static_cast<int>(HistManager::h_general::cuts)]->Fill(1);
        if (pfot.is_ss() && double_tag[kKaon] == K_K ){
          _hm.h_general[static_cast<int>(HistManager::h_general::cuts)]->Fill(2);
          if(sign_check[kKaon]){
              _hm.h_total[static_cast<int>(HistManager::h_total::d0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_d0), maxCtag);
              _hm.h_total[static_cast<int>(HistManager::h_total::d0_sigma)]->Fill(TMath::Sqrt(pfot.KLPFO[ijet].pfo_d0error), maxCtag);
              _hm.h_total[static_cast<int>(HistManager::h_total::d0_sigma_d0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_d0)/TMath::Sqrt(pfot.KLPFO[ijet].pfo_d0error), maxCtag);
            
              _hm.h_total[static_cast<int>(HistManager::h_total::z0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_z0), maxCtag);
              _hm.h_total[static_cast<int>(HistManager::h_total::z0_sigma)]->Fill(TMath::Sqrt(pfot.KLPFO[ijet].pfo_z0error), maxCtag);
              _hm.h_total[static_cast<int>(HistManager::h_total::z0_sigma_z0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_z0)/TMath::Sqrt(pfot.KLPFO[ijet].pfo_z0error), maxCtag);
            

            _hm.h_general[static_cast<int>(HistManager::h_general::cuts)]->Fill(3);
            _hm.h_general[static_cast<int>(HistManager::h_general::ctag_total)]->Fill(ctag[ijet]);
            if(pfot.KLPFO[ijet].pfo_vtx==-1) {
              _hm.h_garbage[static_cast<int>(HistManager::h_total::d0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_d0), maxCtag);
              _hm.h_garbage[static_cast<int>(HistManager::h_total::d0_sigma)]->Fill(TMath::Sqrt(pfot.KLPFO[ijet].pfo_d0error), maxCtag);
              _hm.h_garbage[static_cast<int>(HistManager::h_total::d0_sigma_d0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_d0)/TMath::Sqrt(pfot.KLPFO[ijet].pfo_d0error), maxCtag);
            
              _hm.h_garbage[static_cast<int>(HistManager::h_garbage::z0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_z0), maxCtag);
              _hm.h_garbage[static_cast<int>(HistManager::h_garbage::z0_sigma)]->Fill(TMath::Sqrt(pfot.KLPFO[ijet].pfo_z0error), maxCtag);
              _hm.h_garbage[static_cast<int>(HistManager::h_garbage::z0_sigma_z0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_z0)/TMath::Sqrt(pfot.KLPFO[ijet].pfo_z0error), maxCtag);
            }
            if(pfot.KLPFO[ijet].pfo_vtx==0 and nprong[pfot.KLPFO[ijet].pfo_match][0]>1) {
              _hm.h_general[static_cast<int>(HistManager::h_general::cuts)]->Fill(4);
              _hm.h_general[static_cast<int>(HistManager::h_general::ctag_primary)]->Fill(ctag[ijet]);

              _hm.h_primary[static_cast<int>(HistManager::h_primary::d0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_d0), maxCtag);
              _hm.h_primary[static_cast<int>(HistManager::h_primary::d0_sigma)]->Fill(TMath::Sqrt(pfot.KLPFO[ijet].pfo_d0error), maxCtag);
              _hm.h_primary[static_cast<int>(HistManager::h_primary::d0_sigma_d0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_d0)/TMath::Sqrt(pfot.KLPFO[ijet].pfo_d0error), maxCtag);
            
              _hm.h_primary[static_cast<int>(HistManager::h_primary::z0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_z0), maxCtag);
              _hm.h_primary[static_cast<int>(HistManager::h_primary::z0_sigma)]->Fill(TMath::Sqrt(pfot.KLPFO[ijet].pfo_z0error), maxCtag);
              _hm.h_primary[static_cast<int>(HistManager::h_primary::z0_sigma_z0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_z0)/TMath::Sqrt(pfot.KLPFO[ijet].pfo_z0error), maxCtag);
            }
            if(pfot.KLPFO[ijet].pfo_vtx==1 and nprong[pfot.KLPFO[ijet].pfo_match][1]>1) {
              _hm.h_general[static_cast<int>(HistManager::h_general::cuts)]->Fill(5);
              _hm.h_general[static_cast<int>(HistManager::h_general::ctag_secondary)]->Fill(ctag[ijet]);

              _hm.h_secondary[static_cast<int>(HistManager::h_secondary::d0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_d0), maxCtag);
              _hm.h_secondary[static_cast<int>(HistManager::h_secondary::d0_sigma)]->Fill(TMath::Sqrt(pfot.KLPFO[ijet].pfo_d0error), maxCtag);
              _hm.h_secondary[static_cast<int>(HistManager::h_secondary::d0_sigma_d0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_d0)/TMath::Sqrt(pfot.KLPFO[ijet].pfo_d0error), maxCtag);
            
              _hm.h_secondary[static_cast<int>(HistManager::h_secondary::z0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_z0), maxCtag);
              _hm.h_secondary[static_cast<int>(HistManager::h_secondary::z0_sigma)]->Fill(TMath::Sqrt(pfot.KLPFO[ijet].pfo_z0error), maxCtag);
              _hm.h_secondary[static_cast<int>(HistManager::h_secondary::z0_sigma_z0)]->Fill(TMath::Abs(pfot.KLPFO[ijet].pfo_z0)/TMath::Sqrt(pfot.KLPFO[ijet].pfo_z0error), maxCtag);
            }
          }
        }
      }
    }
  }
  if ( LPFO_checks[kKaon] && pfot.is_ss() && double_tag[kKaon] == K_K ){
    Int_t ineg = -1;

    if( pfot.KLPFO[0].pfo_charge < 0 ){
      ineg = 0;
    }else{
      ineg = 1;
    }
    _hm.h_cos[static_cast<int>(HistManager::h_cos::cos_theta)]->Fill(pfot.KLPFO[ineg].cos, maxCtag);
    if(pfot.KLPFO[ineg].pfo_vtx==0 and nprong[pfot.KLPFO[ineg].pfo_match][0]>1 and pfot.KLPFO[ineg].pfo_d0 < 0.01 and pfot.KLPFO[1-ineg].pfo_d0 < 0.01) {
      _hm.h_cos[static_cast<int>(HistManager::h_cos::cos_theta_primary)]->Fill(pfot.KLPFO[ineg].cos, maxCtag);
    }
    if(sign_check[kKaon]){
      _hm.h_cos[static_cast<int>(HistManager::h_cos::cos_theta_acc)]->Fill(pfot.KLPFO[ineg].qcos, maxCtag);  
      if(pfot.KLPFO[ineg].pfo_vtx==0 and nprong[pfot.KLPFO[ineg].pfo_match][0]>1 and pfot.KLPFO[ineg].pfo_d0 < 0.01 and pfot.KLPFO[1-ineg].pfo_d0 < 0.01) {
        _hm.h_cos[static_cast<int>(HistManager::h_cos::cos_theta_acc_primary)]->Fill(pfot.KLPFO[ineg].qcos, maxCtag);  
      }
    }
    else{
      _hm.h_cos[static_cast<int>(HistManager::h_cos::cos_theta_rej)]->Fill(pfot.KLPFO[ineg].cos, maxCtag);
      _hm.h_cos[static_cast<int>(HistManager::h_cos::cos_theta_rej)]->Fill(-1*pfot.KLPFO[1-ineg].cos, maxCtag);
      if(pfot.KLPFO[ineg].pfo_vtx==0 and nprong[pfot.KLPFO[ineg].pfo_match][0]>1 and pfot.KLPFO[ineg].pfo_d0 < 0.01 and pfot.KLPFO[1-ineg].pfo_d0 < 0.01) {
        _hm.h_cos[static_cast<int>(HistManager::h_cos::cos_theta_rej_primary)]->Fill(pfot.KLPFO[ineg].cos, maxCtag);
        _hm.h_cos[static_cast<int>(HistManager::h_cos::cos_theta_rej_primary)]->Fill(-1*pfot.KLPFO[1-ineg].cos, maxCtag);
      }
    }
  }
}