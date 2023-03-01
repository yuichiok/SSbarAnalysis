
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

void EventAnalyzer::AnalyzeGen()
{
    PFOTools mct( &_mc, _config );
    PolarAngleGen(mct);
}

void EventAnalyzer::AnalyzeGenReco(PFOTools mct, PFOTools pfot)
{
    Mom_Polar_Gen(mct,pfot);
}

void EventAnalyzer::AnalyzeReco(Long64_t entry)
{
  // MC, PFO Analysis
    PFOTools mct( &_mc, _config );
    PFOTools pfot( &_mc, &_pfo, _config );

    // cout << "evt: " << entry << endl;
    AnalyzeGenReco(mct,pfot);

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

  // jet info
  VectorTools jetvt[2];
  for (int i=0; i<2; i++){
    jetvt[i].SetCoordinates(_jet.jet_px[i],_jet.jet_py[i],_jet.jet_pz[i],_jet.jet_E[i]);
  }
  Float_t jet_cos[2]  = { std::cos( jetvt[0].v3().theta() ), std::cos( jetvt[1].v3().theta() ) };
  _hm.h2_jet[_hm.jet_mult_cos_noISR]->Fill( jet_cos[0], _pfo.pfo_n_j1 );
  _hm.h2_jet[_hm.jet_mult_cos_noISR]->Fill( jet_cos[1], _pfo.pfo_n_j2 );

  ////////////////
  // Selections //
  ////////////////

    enum SelectID { kKaon, kPion, kProton };
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

      for ( auto iSPFO_K : pfot.SPFOs_K[ijet] ){
        Bool_t charge_opposite = iSPFO_K.pfo_charge * pfot.KLPFO[ijet].pfo_charge < 0;
        Bool_t momentum_above  = iSPFO_K.p_mag > 10;
        if( charge_opposite && momentum_above ) is_gluon[kKaon][ijet] = true;
      }

      for ( auto iSPFO_Pi : pfot.SPFOs_Pi[ijet] ){
        Bool_t charge_opposite = iSPFO_Pi.pfo_charge * pfot.PiLPFO[ijet].pfo_charge < 0;
        Bool_t momentum_above  = iSPFO_Pi.p_mag > 10;
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
    PDGConfig dEdx_pdg_match;
    
    if ( pfot.isPion(pfot.PiLPFO[0]) && pfot.isPion(pfot.PiLPFO[1]) ) { dEdx_pdg_match = Pi_Pi; }
    else{ dEdx_pdg_match = noKPi; }

    // if ( pfot.isKaon(pfot.KLPFO[0]) && pfot.isKaon(pfot.KLPFO[1]) ) { dEdx_pdg_match = K_K; }
    // else{ dEdx_pdg_match = noKPi; }

    /*
    if     (   pfot.isKaon(pfot.KLPFO[0]) && pfot.isKaon(pfot.KLPFO[1]) )  {  dEdx_pdg_match = K_K;    }
    else if(   pfot.isPion(pfot.KLPFO[0]) && pfot.isPion(pfot.KLPFO[1]) )  {  dEdx_pdg_match = Pi_Pi;  }
    else if( ( pfot.isKaon(pfot.KLPFO[0]) && pfot.isPion(pfot.KLPFO[1]) ) ||
             ( pfot.isKaon(pfot.KLPFO[1]) && pfot.isPion(pfot.KLPFO[0]) ) ){  dEdx_pdg_match = K_Pi;   }
    else{ dEdx_pdg_match = noKPi; }
    */

  // charge config check
    Bool_t charge_check[3] = {false};
    switch ( dEdx_pdg_match )
    {
      case K_K:
        charge_check[kKaon] = pfot.is_charge_config(pfot.kOpposite,pfot.KLPFO[0].pfo_charge,pfot.KLPFO[1].pfo_charge);
        break;
      case Pi_Pi:
        charge_check[kPion] = pfot.is_charge_config(pfot.kOpposite,pfot.PiLPFO[0].pfo_charge,pfot.PiLPFO[1].pfo_charge);
        break;
    default:
      break;
    }

    CutTrigger[kKaon].push_back(charge_check[kKaon]);
    CutTrigger[kPion].push_back(charge_check[kPion]);

  // Try Stability and Purity Calculation here.
    Int_t   *N_Ks  = Gen_Reco_Stats_Stable( mct, pfot, -1, 1 );
    _data.N_K_Gen  = N_Ks[0];
    _data.N_K_PFO  = N_Ks[1];
    _data.N_K_corr = N_Ks[2];

    Float_t *SPs    = Get_Stable_Purity(N_Ks);
    _data.stability = SPs[0];
    _data.purity    = SPs[1];

    Int_t nbins_cos = _hm.h2[_hm.stable_cos]->GetNbinsX();
    TAxis *xaxis    = _hm.h2[_hm.stable_cos]->GetXaxis();
    for ( int ibin=1; ibin<=nbins_cos; ibin++ ){
      Float_t bin_center = xaxis->GetBinCenter(ibin);
      Float_t bin_width  = xaxis->GetBinWidth(ibin);
      Float_t cos_min    = xaxis->GetBinLowEdge(ibin);
      Float_t cos_max    = cos_min + bin_width;
      Int_t   *dN_Ks     = Gen_Reco_Stats_Stable( mct, pfot, cos_min, cos_max );
      Float_t *dSPs      = Get_Stable_Purity(dN_Ks);

      _hm.h1[_hm.gen_N_K_cos]->Fill( bin_center ,dN_Ks[0]);
      _hm.h1[_hm.reco_N_K_cos]->Fill( bin_center ,dN_Ks[1]);
      _hm.h1[_hm.N_K_corr_cos]->Fill( bin_center ,dN_Ks[2]);

      _hm.h2[_hm.stable_cos]->Fill( bin_center ,dSPs[0]);
      _hm.h2[_hm.purity_cos]->Fill( bin_center ,dSPs[1]);

    }

    Int_t nbins_cos2 = _hm.h1[_hm.gen_N_K_cos2]->GetNbinsX();
    TAxis *xaxis2    = _hm.h1[_hm.gen_N_K_cos2]->GetXaxis();
    for ( int ibin=1; ibin<=nbins_cos2; ibin++ ){
      Float_t bin_center = xaxis->GetBinCenter(ibin);
      Float_t bin_width  = xaxis->GetBinWidth(ibin);
      Float_t cos_min    = xaxis->GetBinLowEdge(ibin);
      Float_t cos_max    = cos_min + bin_width;
      Int_t   *dN_Ks     = Gen_Reco_Stats_Stable( mct, pfot, cos_min, cos_max );
      Float_t *dSPs      = Get_Stable_Purity(dN_Ks);
      _hm.h1[_hm.gen_N_K_cos2]->Fill( bin_center ,dN_Ks[0]);
      _hm.h1[_hm.reco_N_K_cos2]->Fill( bin_center ,dN_Ks[1]);
      _hm.h1[_hm.N_K_corr_cos2]->Fill( bin_center ,dN_Ks[2]);
    }

  // Fill Hists can make another class called histogram extractor?
  // CutTrigger = [ Valid_LPFO, Quality, Not_Gluon_K, charge_check ]
  //                          * Quality = {momentum, tpc hits, offset}

  ProcessDoubleTag(pfot,mct,CutTrigger[kPion],dEdx_pdg_match);
  // ProcessDoubleTag(pfot,mct,CutTrigger[kKaon],dEdx_pdg_match);

  // Fill PFO

  // Initialize counters
  Int_t n_reco_particles[3] = {0};
  Int_t n_gen_particles[3]  = {0};

  TH1F * h_n_reco_particles[3];
  for (int i=0; i < 3; i++){
    h_n_reco_particles[i] = new TH1F(TString::Format("h_n_reco_particles%d",i),TString::Format("h_n_reco_particles%d",i),40,-1,1);
  }

  TH1F * h_n_gen_particles[3];
  for (int i=0; i < 3; i++){
    h_n_gen_particles[i] = new TH1F(TString::Format("h_n_gen_particles%d",i),TString::Format("h_n_gen_particles%d",i),40,-1,1);
  }

  std::vector<PFO_Info> PFO_Collection = pfot.Get_Valid_PFOs();
  _data.n_valid_pfo = PFO_Collection.size();

  for ( long unsigned int i=0; i < PFO_Collection.size(); i++ )
  {
    PFO_Info ipfo = PFO_Collection.at(i);

    Count_Particle(ipfo,321,h_n_reco_particles[0],h_n_gen_particles[0]);
    Count_Particle(ipfo,211,h_n_reco_particles[1],h_n_gen_particles[1]);
    Count_Particle(ipfo,2212,h_n_reco_particles[2],h_n_gen_particles[2]);

    // cheat
    switch ( abs(ipfo.pfo_pdgcheat) ) {
      case 321:
        _hm.h2_dEdx[_hm.gen_K_dEdx_p]->Fill(ipfo.p_mag,ipfo.pfo_dedx);
        _hm.h2_dEdx[_hm.gen_K_KdEdx_dist_cos]->Fill(ipfo.cos,ipfo.pfo_piddedx_k_dedxdist);
        _hm.h2_dEdx[_hm.gen_K_PidEdx_dist_cos]->Fill(ipfo.cos,ipfo.pfo_piddedx_pi_dedxdist);
        break;
      case 211:
        _hm.h2_dEdx[_hm.gen_pi_dEdx_p]->Fill(ipfo.p_mag,ipfo.pfo_dedx);
        _hm.h2_dEdx[_hm.gen_pi_KdEdx_dist_cos]->Fill(ipfo.cos,ipfo.pfo_piddedx_k_dedxdist);
        _hm.h2_dEdx[_hm.gen_pi_PidEdx_dist_cos]->Fill(ipfo.cos,ipfo.pfo_piddedx_pi_dedxdist);
        break;
      case 2212:
        _hm.h2_dEdx[_hm.gen_p_dEdx_p]->Fill(ipfo.p_mag,ipfo.pfo_dedx);
        _hm.h2_dEdx[_hm.gen_p_KdEdx_dist_cos]->Fill(ipfo.cos,ipfo.pfo_piddedx_k_dedxdist);
        _hm.h2_dEdx[_hm.gen_p_PidEdx_dist_cos]->Fill(ipfo.cos,ipfo.pfo_piddedx_pi_dedxdist);
        break;
      default:
        break;
    }

  }


  // Access cheated Kaon information
  if( pfot.PFO_cheat_Ks[0].size() && pfot.PFO_cheat_Ks[1].size() ){

    vector<Bool_t> Cheat_K_CutTrigger;
    // quality check
    Bool_t cheat_K_double_quality    = true;
    for ( auto iLPFO : pfot.cheat_KLPFO ){
      if( !pfot.LPFO_Quality_checks(iLPFO) ){
        cheat_K_double_quality = false;
        break;
      }
    }
    Cheat_K_CutTrigger.push_back(cheat_K_double_quality);

    // SPFO opposite check
    Bool_t is_cheat_gluon_K[2] = {0};
    Bool_t is_there_a_cheat_gluon_K = false;
    for ( int ijet=0; ijet<2; ijet++){
      if( pfot.SPFOs_cheat_K[ijet].size() ){
        for ( auto iSPFO_K : pfot.SPFOs_cheat_K[ijet] ){
          Bool_t charge_opposite = iSPFO_K.pfo_pdgcheat * pfot.cheat_KLPFO[ijet].pfo_pdgcheat < 0;
          Bool_t momentum_above  = iSPFO_K.p_mag > 10;
          if( charge_opposite && momentum_above ) is_cheat_gluon_K[ijet] = true;
        }
      }
    }
    for ( auto ibool : is_cheat_gluon_K ){
      if( ibool ) is_there_a_cheat_gluon_K = true;
    }
    Cheat_K_CutTrigger.push_back(!is_there_a_cheat_gluon_K);

    // charge check
    Bool_t cheat_K_charge_check = pfot.is_charge_config(pfot.kOpposite,pfot.cheat_KLPFO[0].pfo_pdgcheat,pfot.cheat_KLPFO[1].pfo_pdgcheat);
    Cheat_K_CutTrigger.push_back(cheat_K_charge_check);

    Bool_t cheat_K_all_checks = true;
    for (auto ibool : Cheat_K_CutTrigger){
      if (!ibool) {
        cheat_K_all_checks = false;
        break;
      }
    }

    if ( cheat_K_all_checks ){

      Float_t cheat_gen_angle = abs(pfot.cheat_KLPFO[0].cos) * sgn( -_mc.mc_quark_charge[0] ) * mct.mc_quark[0].cos / abs(mct.mc_quark[0].cos);
      _hm.h1[_hm.cheat_K_cos]->Fill(pfot.cheat_KLPFO[0].cos);
      _hm.h1[_hm.cheat_K_qcos]->Fill(cheat_gen_angle);

    }

  }

  // Fill Event by Event hists
  Float_t * particle_ratios_reco = Particle_Ratios( h_n_reco_particles, 0 );
  Float_t * particle_ratios_gen  = Particle_Ratios( h_n_gen_particles, 1 );

  _hm.h1_particle_ratio[_hm.K_rate_reco]->Fill(particle_ratios_reco[0]);
  _hm.h1_particle_ratio[_hm.pi_rate_reco]->Fill(particle_ratios_reco[1]);
  _hm.h1_particle_ratio[_hm.p_rate_reco]->Fill(particle_ratios_reco[2]);

  _hm.h1_particle_ratio[_hm.K_rate_gen]->Fill(particle_ratios_gen[0]);
  _hm.h1_particle_ratio[_hm.pi_rate_gen]->Fill(particle_ratios_gen[1]);
  _hm.h1_particle_ratio[_hm.p_rate_gen]->Fill(particle_ratios_gen[2]);



  _hm.h2[_hm.nK_gen_reco]->Fill(n_reco_particles[0],n_gen_particles[0]);
  _hm.h2[_hm.npi_gen_reco]->Fill(n_reco_particles[1],n_gen_particles[1]);
  _hm.h2[_hm.np_gen_reco]->Fill(n_reco_particles[2],n_gen_particles[2]);

  if(_eve.eve_valid_lpfo){

    for (auto iLPFO : pfot.KLPFO){
      if( pfot.isKaon(iLPFO) ) _hm.h1[_hm.lpfo_reco_K_mom]->Fill(iLPFO.p_mag);
    }

    for (int i=0; i<2; i++){
      if( abs(_stats_lpfo.lpfo_pdgcheat[i]) == 321 ) _hm.h1[_hm.lpfo_gen_K_mom]->Fill(pfot.KLPFO[i].p_mag);
    }
  }

  _hTree->Fill();

  ClearStructs();
  for (int i=0; i < 3; i++){
    delete h_n_reco_particles[i];
    delete h_n_gen_particles[i];
  }

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

Int_t *EventAnalyzer::Gen_Reco_Stats_Stable( PFOTools mct, PFOTools pfot, Float_t cos_min, Float_t cos_max )
{
  std::vector<PFO_Info> PFO_Collection = pfot.Get_Valid_PFOs();

  Float_t p_min = _anCfg.PFO_p_min;

  std::vector<PFO_Info> PFO_K_Collection;
  std::vector<PFO_Info> Gen_K_Collection;
  for ( auto iPFO : PFO_Collection ){
    Bool_t cos_range = (cos_min < iPFO.cos && iPFO.cos < cos_max );
    Bool_t p_range   = p_min < iPFO.p_mag;
    if( PFOTools::isKaon(iPFO) && cos_range && p_range )   PFO_K_Collection.push_back(iPFO);
    if( abs(iPFO.pfo_pdgcheat) == 321 && cos_range && p_range ) Gen_K_Collection.push_back(iPFO);
  }

  Int_t N_K_corr  = 0;

  // Float_t cos_r = 0.02;
  Float_t cos_r = 0.37;
  std::vector<PFO_Info> PFO_K_Remain = PFO_K_Collection;

  for ( auto igen : Gen_K_Collection ){
    Float_t min_cos_diff   = 1000.0;
    Int_t   i_min_cos_diff = -1;

    Int_t counter = 0;
    for ( auto iremain : PFO_K_Remain ){
      // Float_t cos_diff = igen.cos - iremain.cos;
      Float_t cos_diff = 1.0 - std::cos( igen.vt.v3().Theta() - iremain.vt.v3().Theta() );
      if( cos_diff < min_cos_diff ) {
        min_cos_diff = cos_diff;
        i_min_cos_diff = counter;
      }
      counter++;
    }

    // Fill min_cos_diff
    _hm.h1[_hm.gen_reco_K_sep_cos]->Fill(min_cos_diff);

    if( min_cos_diff < cos_r ) {
      N_K_corr++;
      PFO_K_Remain.erase( PFO_K_Remain.begin() + i_min_cos_diff );
    }

  }

  static Int_t N_array[3] = {0};
  N_array[0] = Gen_K_Collection.size();
  N_array[1] = PFO_K_Collection.size();
  N_array[2] = N_K_corr;

  return N_array;

}

Int_t *EventAnalyzer::Gen_Reco_Stats_Cheat( PFOTools mct, PFOTools pfot, Float_t cos_min, Float_t cos_max )
{
  std::vector<PFO_Info> PFO_Collection;
  std::vector<PFO_Info> jet[2] = { pfot.GetJet(0), pfot.GetJet(1) };

  PFO_Collection.reserve( jet[0].size() + jet[1].size() );
  PFO_Collection.insert( PFO_Collection.begin(), jet[0].begin(), jet[0].end() );
  PFO_Collection.insert( PFO_Collection.end(), jet[1].begin(), jet[1].end() );

  Float_t p_min = _anCfg.PFO_p_min;

  std::vector<PFO_Info> PFO_K_Collection;
  for ( auto iPFO : PFO_Collection ){
    Bool_t cos_range = (cos_min < iPFO.cos && iPFO.cos < cos_max );
    Bool_t p_range   = p_min < iPFO.p_mag;
    if( PFOTools::isKaon(iPFO) && cos_range && p_range ) PFO_K_Collection.push_back(iPFO);
  }

  std::vector<PFO_Info> PFO_Cheat_K_Collection;
  for ( auto iPFO : PFO_Collection ){
    Bool_t cos_range = (cos_min < iPFO.cos && iPFO.cos < cos_max );
    Bool_t p_range   = p_min < iPFO.p_mag;
    Bool_t k_cheat   = ( abs(iPFO.pfo_pdgcheat) == 321);
    if( k_cheat && cos_range && p_range ) PFO_Cheat_K_Collection.push_back(iPFO);
  }

  Int_t N_K_corr  = 0;

  // Float_t cos_r = 0.02;
  Float_t cos_r = 0.37;
  std::vector<PFO_Info> PFO_K_Remain = PFO_K_Collection;

  for ( auto igen : PFO_Cheat_K_Collection ){
    Float_t min_cos_diff   = 1000.0;
    Int_t   i_min_cos_diff = -1;

    Int_t counter = 0;
    for ( auto iremain : PFO_K_Remain ){
      // Float_t cos_diff = igen.cos - iremain.cos;
      Float_t cos_diff = 1.0 - std::cos( igen.vt.v3().Theta() - iremain.vt.v3().Theta() );
      if( cos_diff < min_cos_diff ) {
        min_cos_diff = cos_diff;
        i_min_cos_diff = counter;
      }
      counter++;
    }

    // Fill min_cos_diff
    _hm.h1[_hm.gen_reco_K_sep_cos]->Fill(min_cos_diff);

    if( min_cos_diff < cos_r ) {
      N_K_corr++;
      PFO_K_Remain.erase( PFO_K_Remain.begin() + i_min_cos_diff );
    }

  }

  static Int_t N_array[3] = {0};
  N_array[0] = PFO_Cheat_K_Collection.size();
  N_array[1] = PFO_K_Collection.size();
  N_array[2] = N_K_corr;

  return N_array;

}

Float_t *EventAnalyzer::Get_Stable_Purity( Int_t *N_Ks )
{
  // 0: N Gen Kaon, 1: N PFO Kaon, 2: N Correct Kaon
  static Float_t SP_array[2] = {-1,-1};
  if( N_Ks[0] )  SP_array[0] = (Float_t) N_Ks[2] / (Float_t) N_Ks[0];
  if( N_Ks[1] )  SP_array[1] = (Float_t) N_Ks[2] / (Float_t) N_Ks[1];

  return SP_array;

}

void EventAnalyzer::Count_Particle( PFO_Info ipfo, Int_t pdg, TH1F *h_n_reco, TH1F *h_n_gen )
{
  if ( abs(ipfo.pfo_pdgcheat) == pdg ) h_n_gen->Fill(ipfo.qcos);

  switch ( pdg ){
    case 321:
      if ( PFOTools::isKaon(ipfo) ) h_n_reco->Fill(ipfo.qcos);
      break;
    case 211:
      if ( PFOTools::isPion(ipfo) ) h_n_reco->Fill(ipfo.qcos);
      break;
    case 2212:
      if ( PFOTools::isProton(ipfo) ) h_n_reco->Fill(ipfo.qcos);
      break;
    default:
      break;
  }
}

Float_t *EventAnalyzer::Particle_Ratios( TH1F *h_n_particles[], Int_t mode )
{
  Int_t N_Particles[3]  = {0};
  Int_t N_sum = 0;
  for ( int i=0; i < 3; i++){
    N_Particles[i] = h_n_particles[i]->GetEntries();
    N_sum += N_Particles[i];
  }

  TH1F *h_sum = (TH1F*) h_n_particles[0]->Clone();
  h_sum->Add(h_n_particles[1]);
  h_sum->Add(h_n_particles[2]);


  // cos theta ratio
  for ( int i=0; i < 3; i++){
    TH1F *h_div = (TH1F*) h_n_particles[i]->Clone();
    h_div->Divide(h_sum);

    for ( int ibin=1; ibin <= h_div->GetNbinsX(); ibin++ ){
      Float_t eff = h_div->GetBinContent(ibin);
      Float_t cos = h_div->GetXaxis()->GetBinCenter(ibin);

      switch (i) {
        case 0:
          if(mode==0) _hm.h2_particle_ratio_cos[_hm.K_rate_cos_reco]->Fill(cos,eff);
          if(mode==1) _hm.h2_particle_ratio_cos[_hm.K_rate_cos_gen]->Fill(cos,eff);
          break;
        case 1:
          if(mode==0) _hm.h2_particle_ratio_cos[_hm.pi_rate_cos_reco]->Fill(cos,eff);
          if(mode==1) _hm.h2_particle_ratio_cos[_hm.pi_rate_cos_gen]->Fill(cos,eff);
          break;
        case 2:
          if(mode==0) _hm.h2_particle_ratio_cos[_hm.p_rate_cos_reco]->Fill(cos,eff);
          if(mode==1) _hm.h2_particle_ratio_cos[_hm.p_rate_cos_gen]->Fill(cos,eff);
          break;

        default:
          break;
      }

    }

    delete h_div;

  }

  delete h_sum;

  // overall ratio
  static Float_t ratio_arr[3]  = {-1,-1,-1};
  static Float_t ratio_arr0[3] = {-1,-1,-1};

  if( N_sum == 0 ) return ratio_arr0;
  for( int i=0; i<3; i++ ){
    ratio_arr[i] = (Float_t) N_Particles[i] / (Float_t) N_sum;
  }

  return ratio_arr;

}

void EventAnalyzer::PolarAngleGen(PFOTools mct)
{
  // Gen QQbar
  for ( auto iq : mct.mc_quark ){
    _hm.h1[_hm.gen_q_cos]->Fill(iq.cos);
    _hm.h1[_hm.gen_q_qcos]->Fill(iq.qcos);
  }

  // Gen K
  for ( int istable=0; istable < _mc.mc_stable_n; istable++ ){
    if(abs(_mc.mc_stable_pdg[istable]) == 321) {
      _hm.h1[_hm.gen_K_cos]->Fill(mct.mc_stable[istable].cos);
      _hm.h1[_hm.gen_K_qcos]->Fill(mct.mc_stable[istable].qcos);
    }
  }

}

void EventAnalyzer::Mom_Polar_Gen(PFOTools mct, PFOTools pfot)
{
  Int_t cnt_gen_K = 0;
  Int_t cnt_reco_K = 0;

  // Gen K
  for ( int istable=0; istable < _mc.mc_stable_n; istable++ ){

    if(abs(_mc.mc_stable_pdg[istable]) == 321 && 20 < mct.mc_stable[istable].p_mag) {
      cnt_gen_K++;
      _hm.h2[_hm.gen_K_p_cos]->Fill(mct.mc_stable[istable].cos,mct.mc_stable[istable].p_mag);
    }
  }

  // PFO cheat K
  for ( int ipfo=0; ipfo < _pfo.pfo_n; ipfo++ ){
    VectorTools vt(_pfo.pfo_px[ipfo], _pfo.pfo_py[ipfo], _pfo.pfo_pz[ipfo], _pfo.pfo_E[ipfo]);
    Float_t pfo_p_mag = (Float_t) vt.v3().R();
    Float_t pfo_cos   = std::cos(vt.v3().Theta());

    if(abs(_pfo.pfo_pdgcheat[ipfo]) == 321 && 20 < pfo_p_mag) {
      cnt_reco_K++;
      _hm.h2[_hm.reco_K_p_cos]->Fill(pfo_cos,pfo_p_mag);
    }
  }

}

void EventAnalyzer::ProcessDoubleTag(PFOTools pfot, PFOTools mct, vector<Bool_t> cuts, PDGConfig double_tag)
{
  Bool_t LPFO_checks = true;
  Int_t vec_size = cuts.size();
  for (int i=0; i<vec_size-1; i++){

    if (!cuts.at(i)){
      LPFO_checks = false;
      break;
    }

  }

  Bool_t sign_check = cuts.back();

  // Reco K_K
  if(LPFO_checks){

    Int_t ineg = -1;
      
    switch ( double_tag ){
      
      case K_K:

        if( pfot.KLPFO[0].pfo_charge < 0 ){
          ineg = 0;
        }else{
          ineg = 1;
        }

        if(sign_check){

          _hm.h1[_hm.reco_K_cos]->Fill( pfot.KLPFO[ineg].cos );
          _hm.h1[_hm.reco_K_qcos]->Fill( pfot.KLPFO[ineg].qcos );
          _hm.h1[_hm.reco_K_scos]->Fill( abs(pfot.KLPFO[ineg].cos) * sgn( -_mc.mc_quark_charge[0] ) * mct.mc_quark[0].cos / abs(mct.mc_quark[0].cos) );

          // cheat
          switch ( abs(pfot.KLPFO[ineg].pfo_pdgcheat) ) {
            case 321:
              _hm.h2_dEdx[_hm.gen_K_reco_K_dEdx_p]->Fill(pfot.KLPFO[ineg].p_mag,pfot.KLPFO[ineg].pfo_dedx);
              _hm.h2_dEdx[_hm.gen_K_reco_K_KdEdx_dist_cos]->Fill(pfot.KLPFO[ineg].cos,pfot.KLPFO[ineg].pfo_piddedx_k_dedxdist);
              break;
            case 211:
              _hm.h2_dEdx[_hm.gen_pi_reco_K_dEdx_p]->Fill(pfot.KLPFO[ineg].p_mag,pfot.KLPFO[ineg].pfo_dedx);
              _hm.h2_dEdx[_hm.gen_pi_reco_K_KdEdx_dist_cos]->Fill(pfot.KLPFO[ineg].cos,pfot.KLPFO[ineg].pfo_piddedx_k_dedxdist);
              break;
            case 2212:
              _hm.h2_dEdx[_hm.gen_p_reco_K_dEdx_p]->Fill(pfot.KLPFO[ineg].p_mag,pfot.KLPFO[ineg].pfo_dedx);
              _hm.h2_dEdx[_hm.gen_p_reco_K_KdEdx_dist_cos]->Fill(pfot.KLPFO[ineg].cos,pfot.KLPFO[ineg].pfo_piddedx_k_dedxdist);
              break;
            default:
              break;
          }

          _hm.h1_pq[_hm.acc_KK]->Fill( pfot.KLPFO[ineg].qcos );

        }else{

          _hm.h1_pq[_hm.rej_KK]->Fill( pfot.KLPFO[ineg].cos );
          _hm.h1_pq[_hm.rej_KK]->Fill( -pfot.KLPFO[1-ineg].cos );
        
        }

        break;

      case Pi_Pi:

        if( pfot.PiLPFO[0].pfo_charge < 0 ){
          ineg = 0;
        }else{
          ineg = 1;
        }

        if(sign_check){

          _hm.h1[_hm.reco_Pi_cos]->Fill( pfot.PiLPFO[ineg].cos );
          _hm.h1[_hm.reco_Pi_qcos]->Fill( pfot.PiLPFO[ineg].qcos );
          _hm.h1[_hm.reco_Pi_scos]->Fill( abs(pfot.PiLPFO[ineg].cos) * sgn( -_mc.mc_quark_charge[1] ) * mct.mc_quark[1].cos / abs(mct.mc_quark[1].cos) );

          // cheat
          switch ( abs(pfot.PiLPFO[ineg].pfo_pdgcheat) ) {
            case 321:
              _hm.h2_dEdx[_hm.gen_K_reco_Pi_dEdx_p]->Fill(pfot.PiLPFO[ineg].p_mag,pfot.PiLPFO[ineg].pfo_dedx);
              _hm.h2_dEdx[_hm.gen_K_reco_Pi_PidEdx_dist_cos]->Fill(pfot.PiLPFO[ineg].cos,pfot.PiLPFO[ineg].pfo_piddedx_pi_dedxdist);
              break;
            case 211:
              _hm.h2_dEdx[_hm.gen_pi_reco_Pi_dEdx_p]->Fill(pfot.PiLPFO[ineg].p_mag,pfot.PiLPFO[ineg].pfo_dedx);
              _hm.h2_dEdx[_hm.gen_pi_reco_Pi_PidEdx_dist_cos]->Fill(pfot.PiLPFO[ineg].cos,pfot.PiLPFO[ineg].pfo_piddedx_pi_dedxdist);
              break;
            case 2212:
              _hm.h2_dEdx[_hm.gen_p_reco_Pi_dEdx_p]->Fill(pfot.PiLPFO[ineg].p_mag,pfot.PiLPFO[ineg].pfo_dedx);
              _hm.h2_dEdx[_hm.gen_p_reco_Pi_PidEdx_dist_cos]->Fill(pfot.PiLPFO[ineg].cos,pfot.PiLPFO[ineg].pfo_piddedx_pi_dedxdist);
              break;
            default:
              break;
          }

          _hm.h1_pq[_hm.acc_PiPi]->Fill( pfot.PiLPFO[ineg].qcos );

        }else{

          _hm.h1_pq[_hm.rej_PiPi]->Fill( pfot.PiLPFO[ineg].cos );
          _hm.h1_pq[_hm.rej_PiPi]->Fill( -pfot.PiLPFO[1-ineg].cos );
        
        }

        break;

      default:
        break;
    }
    
  }

}

void EventAnalyzer::PolarAngle_acc_rej(PFOTools pfot, vector<Bool_t> cuts, Bool_t ss_config)
{
  if (cuts.size()!=4) return;

  // LPFO present, PFO Quality Good
  for (int icut=0; icut < 3; icut++){
    if (!cuts.at(icut)) return;
  }

  // K_K config
  if(!ss_config) return;

  // Last element of cuts vector is charge comparison
  if(cuts.back()){
    _hm.h1_pq[_hm.acc_KK]->Fill( pfot.KLPFO[0].qcos );
  }else{
    _hm.h1_pq[_hm.rej_KK]->Fill( pfot.KLPFO[0].cos );
    _hm.h1_pq[_hm.rej_KK]->Fill( -pfot.KLPFO[0].cos );
  }

}

void EventAnalyzer::Jet_sum_n_acol()
{
  VectorTools jetvt[2];
  for (int i=0; i<2; i++){
    jetvt[i].SetCoordinates(_jet.jet_px[i],_jet.jet_py[i],_jet.jet_pz[i],_jet.jet_E[i]);
  }
  Float_t cosacol = std::cos( VectorTools::GetThetaBetween( jetvt[0].v3(), jetvt[1].v3() ) );
  Float_t jet_cos[2]  = { std::cos( jetvt[0].v3().theta() ), std::cos( jetvt[1].v3().theta() ) };

  _data.sum_jet_E = _jet.jet_E[0] + _jet.jet_E[1];
  _data.jet_acol  = cosacol;

  _hm.h1[_hm.reco_sum_jetE]->Fill( _data.sum_jet_E );
  _hm.h1[_hm.reco_jet_sep]->Fill( _data.jet_acol );

  _hm.h2_jet[_hm.jet_mult_cos]->Fill( jet_cos[0], _pfo.pfo_n_j1 );
  _hm.h2_jet[_hm.jet_mult_cos]->Fill( jet_cos[1], _pfo.pfo_n_j2 );

}