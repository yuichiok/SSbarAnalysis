
/*------------------------------------------------------------------------------
EventAnalyzer.cpp
 Created : 2022-09-05  okugawa
------------------------------------------------------------------------------*/
#include "EventAnalyzer.hh"

#include <TBranch.h>
#include <TLeaf.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <utility>

using std::cout;   using std::endl;
typedef unsigned int Index;

namespace QQbarAnalysis
{
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
    reader.InitializeVTXReadTree(fChain, _vtx, _branch);
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

  void EventAnalyzer::InitHists()
  {
    _hm.InitializeHists();
  }

  void EventAnalyzer::WriteFile()
  {
    // Write Tree
    _hm.WriteLists(_hfile);
    _hfile->Close();
  }

  void EventAnalyzer::AnalyzeGen()
  {
    PFOTools mct( &_mc, _config );
    PolarAngleGen(mct);
  }

  void EventAnalyzer::AnalyzeReco(Long64_t entry)
  {
    // MC, PFO Analysis
    PFOTools mct( &_mc, _config );
    PFOTools pfot( &_mc, &_pfo, _config );

    ientry = entry;
    // cout << "evt: " << entry << endl;

    // jet info
    VectorTools jetvt[2];
    for (int i=0; i<2; i++){
      jetvt[i].SetCoordinates(_jet.jet_px[i],_jet.jet_py[i],_jet.jet_pz[i],_jet.jet_E[i]);
    }
    Float_t jet_cos[2]  = { std::cos( jetvt[0].v3().theta() ), std::cos( jetvt[1].v3().theta() ) };
    _hm.h2_jet[_hm.jet_mult_cos_noISR]->Fill( jet_cos[0], _jet.jet_npfo[0] );
    _hm.h2_jet[_hm.jet_mult_cos_noISR]->Fill( jet_cos[1], _jet.jet_npfo[1] );

    ////////////////
    // Selections //
    ////////////////

    vector<Bool_t> CutTrigger[3];    
    unordered_map< TString, unordered_map<TString, Bool_t> > CutTriggerMap; // [particle][cutname]

    for( const auto &[i_lmode, val_LPFO]: pfot.LPFO_ ){

      unordered_map<int, PFO_Info> LPFOs = val_LPFO;

      // check momentum
      CutTriggerMap[i_lmode]["momentum"] = pfot.is_momentum( LPFOs.at(0), _anCfg.LPFO_p_min, _anCfg.LPFO_p_max ) &&
                                           pfot.is_momentum( LPFOs.at(1), _anCfg.LPFO_p_min, _anCfg.LPFO_p_max );

      // check tpc hits
      CutTriggerMap[i_lmode]["tpc_hits"] = pfot.is_tpc_hits( LPFOs.at(0), _anCfg.PFO_TPCHits_min ) &&
                                           pfot.is_tpc_hits( LPFOs.at(1), _anCfg.PFO_TPCHits_min );

      // check offsets
      CutTriggerMap[i_lmode]["offset"] = pfot.is_offset_small( LPFOs.at(0), _anCfg.PFO_offset_max ) &&
                                         pfot.is_offset_small( LPFOs.at(1), _anCfg.PFO_offset_max );
      
      // SPFO opposite check
      vector<Bool_t> is_SPFO_charge_opposite(2,false);
      for ( int ijet=0; ijet<2; ijet++ ){
        if( pfot.PFO_jet[ijet].size() > 1 ){
          for ( auto iSPFO : pfot.SPFOs_.at(i_lmode).at(ijet) ){
            Bool_t charge_opposite = iSPFO.pfo_charge * LPFOs.at(ijet).pfo_charge < 0;
            Bool_t momentum_above  = iSPFO.p_mag > 10;
            if( charge_opposite && momentum_above ) is_SPFO_charge_opposite.at(ijet) = true;
          }
        }
      }
      CutTriggerMap[i_lmode]["SPFO"] = std::none_of(is_SPFO_charge_opposite.begin(), is_SPFO_charge_opposite.end(), [](bool v) { return v; });

      // Charge opposite check
      CutTriggerMap[i_lmode]["charge"] = pfot.is_charge_config(pfot.kOpposite,LPFOs.at(0).pfo_charge,LPFOs.at(1).pfo_charge);

      // Particle ID both sides
      CutTriggerMap[i_lmode]["PID"]    = pfot.is_PID_config( i_lmode );

      // Higher LPFO momentum
      CutTriggerMap[i_lmode]["LPFO_higher_p"] = pfot.is_high_LPFO( i_lmode );

    }

///////////////////////////////////////////
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
    PDGConfig dEdx_pdg_match[3] = {noKPi, noKPi, noKPi};
    
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
///////////////////////////////////////////

    // Try Stability and Purity Calculation here.
    // Kaon Efficiency
    Int_t nbins_cos_K = _hm.h2[_hm.stable_K_cos]->GetNbinsX();
    TAxis *xaxis_K    = _hm.h2[_hm.stable_K_cos]->GetXaxis();
    for ( int ibin=1; ibin<=nbins_cos_K; ibin++ ){
      Float_t bin_center = xaxis_K->GetBinCenter(ibin);
      Float_t bin_width  = xaxis_K->GetBinWidth(ibin);
      Float_t cos_min    = xaxis_K->GetBinLowEdge(ibin);
      Float_t cos_max    = cos_min + bin_width;
      Int_t   *dN_Ks     = Gen_Reco_Stats_Stable( mct, pfot, kKaon, cos_min, cos_max );

      _hm.h1[_hm.gen_N_K_cos]->Fill( bin_center ,dN_Ks[0]);
      _hm.h1[_hm.reco_N_K_cos]->Fill( bin_center ,dN_Ks[1]);
      _hm.h1[_hm.N_K_corr_cos]->Fill( bin_center ,dN_Ks[2]);
    }

    // Pion Efficiency
    Int_t nbins_cos_Pi = _hm.h2[_hm.stable_Pi_cos]->GetNbinsX();
    TAxis *xaxis_Pi    = _hm.h2[_hm.stable_Pi_cos]->GetXaxis();
    for ( int ibin=1; ibin<=nbins_cos_K; ibin++ ){
      Float_t bin_center = xaxis_Pi->GetBinCenter(ibin);
      Float_t bin_width  = xaxis_Pi->GetBinWidth(ibin);
      Float_t cos_min    = xaxis_Pi->GetBinLowEdge(ibin);
      Float_t cos_max    = cos_min + bin_width;
      Int_t   *dN_Pis    = Gen_Reco_Stats_Stable( mct, pfot, kPion, cos_min, cos_max );

      _hm.h1[_hm.gen_N_Pi_cos]->Fill( bin_center ,dN_Pis[0]);
      _hm.h1[_hm.reco_N_Pi_cos]->Fill( bin_center ,dN_Pis[1]);
      _hm.h1[_hm.N_Pi_corr_cos]->Fill( bin_center ,dN_Pis[2]);
    }

    // Fill Hists can make another class called histogram extractor?
    // CutTrigger = [ Valid_LPFO, Quality, Not_Gluon_K, charge_check ]
    //                          * Quality = {momentum, tpc hits, offset}

    ProcessDoubleTag(pfot,mct,CutTriggerMap);
    ProcessDoubleTag(pfot,mct,CutTrigger,dEdx_pdg_match);

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

    std::vector<PFO_Info> PFO_Collection = pfot.Valid_PFOs;
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
          _hm.h2_dEdx_nomap[_hm.gen_ipart_dEdx_p][kKaon]->Fill(ipfo.p_mag,ipfo.pfo_dedx);
          _hm.h2_dEdx_nomap[_hm.gen_ipart_KdEdx_dist_cos][kKaon]->Fill(ipfo.cos,ipfo.pfo_piddedx_k_dedxdist);
          _hm.h2_dEdx_nomap[_hm.gen_ipart_PidEdx_dist_cos][kKaon]->Fill(ipfo.cos,ipfo.pfo_piddedx_pi_dedxdist);
          break;
        case 211:
          _hm.h2_dEdx_nomap[_hm.gen_ipart_dEdx_p][kPion]->Fill(ipfo.p_mag,ipfo.pfo_dedx);
          _hm.h2_dEdx_nomap[_hm.gen_ipart_KdEdx_dist_cos][kPion]->Fill(ipfo.cos,ipfo.pfo_piddedx_k_dedxdist);
          _hm.h2_dEdx_nomap[_hm.gen_ipart_PidEdx_dist_cos][kPion]->Fill(ipfo.cos,ipfo.pfo_piddedx_pi_dedxdist);
          break;
        case 2212:
          _hm.h2_dEdx_nomap[_hm.gen_ipart_dEdx_p][kProton]->Fill(ipfo.p_mag,ipfo.pfo_dedx);
          _hm.h2_dEdx_nomap[_hm.gen_ipart_KdEdx_dist_cos][kProton]->Fill(ipfo.cos,ipfo.pfo_piddedx_k_dedxdist);
          _hm.h2_dEdx_nomap[_hm.gen_ipart_PidEdx_dist_cos][kProton]->Fill(ipfo.cos,ipfo.pfo_piddedx_pi_dedxdist);
          break;
        case 11:
          _hm.h2_dEdx_nomap[_hm.gen_ipart_dEdx_p][kElectron]->Fill(ipfo.p_mag,ipfo.pfo_dedx);
          _hm.h2_dEdx_nomap[_hm.gen_ipart_KdEdx_dist_cos][kElectron]->Fill(ipfo.cos,ipfo.pfo_piddedx_k_dedxdist);
          _hm.h2_dEdx_nomap[_hm.gen_ipart_PidEdx_dist_cos][kElectron]->Fill(ipfo.cos,ipfo.pfo_piddedx_pi_dedxdist);
          break;
        case 13:
          _hm.h2_dEdx_nomap[_hm.gen_ipart_dEdx_p][kMuon]->Fill(ipfo.p_mag,ipfo.pfo_dedx);
          _hm.h2_dEdx_nomap[_hm.gen_ipart_KdEdx_dist_cos][kMuon]->Fill(ipfo.cos,ipfo.pfo_piddedx_k_dedxdist);
          _hm.h2_dEdx_nomap[_hm.gen_ipart_PidEdx_dist_cos][kMuon]->Fill(ipfo.cos,ipfo.pfo_piddedx_pi_dedxdist);
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

        Int_t ineg = -1;

        if( pfot.cheat_KLPFO[0].pfo_pdgcheat < 0 ){
          ineg = 0;
        }else{
          ineg = 1;
        }

        Float_t cheat_gen_angle = abs(pfot.cheat_KLPFO[ineg].cos) * sgn( -_mc.mc_quark_charge[0] ) * mct.mc_quark[0].cos / abs(mct.mc_quark[0].cos);
        _hm.h1[_hm.cheat_K_cos]->Fill(pfot.cheat_KLPFO[ineg].cos);
        _hm.h1[_hm.cheat_K_qcos]->Fill(cheat_gen_angle);

      }

    }

    // Access cheated Pion information
    if( pfot.PFO_cheat_Pis[0].size() && pfot.PFO_cheat_Pis[1].size() ){

      vector<Bool_t> Cheat_Pi_CutTrigger;
      // quality check
      Bool_t cheat_Pi_double_quality    = true;
      for ( auto iLPFO : pfot.cheat_PiLPFO ){
        if( !pfot.LPFO_Quality_checks(iLPFO) ){
          cheat_Pi_double_quality = false;
          break;
        }
      }
      Cheat_Pi_CutTrigger.push_back(cheat_Pi_double_quality);

      // SPFO opposite check
      Bool_t is_cheat_gluon_Pi[2] = {0};
      Bool_t is_there_a_cheat_gluon_Pi = false;
      for ( int ijet=0; ijet<2; ijet++){
        if( pfot.SPFOs_cheat_Pi[ijet].size() ){
          for ( auto iSPFO_Pi : pfot.SPFOs_cheat_Pi[ijet] ){
            Bool_t charge_opposite = iSPFO_Pi.pfo_pdgcheat * pfot.cheat_PiLPFO[ijet].pfo_pdgcheat < 0;
            Bool_t momentum_above  = iSPFO_Pi.p_mag > 10;
            if( charge_opposite && momentum_above ) is_cheat_gluon_Pi[ijet] = true;
          }
        }
      }
      for ( auto ibool : is_cheat_gluon_Pi ){
        if( ibool ) is_there_a_cheat_gluon_Pi = true;
      }
      Cheat_Pi_CutTrigger.push_back(!is_there_a_cheat_gluon_Pi);

      // charge check
      Bool_t cheat_Pi_charge_check = pfot.is_charge_config(pfot.kOpposite,pfot.cheat_PiLPFO[0].pfo_pdgcheat,pfot.cheat_PiLPFO[1].pfo_pdgcheat);
      Cheat_Pi_CutTrigger.push_back(cheat_Pi_charge_check);

      Bool_t cheat_Pi_all_checks = true;
      for (auto ibool : Cheat_Pi_CutTrigger){
        if (!ibool) {
          cheat_Pi_all_checks = false;
          break;
        }
      }

      if ( cheat_Pi_all_checks ){

        Int_t ineg = -1;

        if( pfot.cheat_PiLPFO[0].pfo_pdgcheat < 0 ){
          ineg = 0;
        }else{
          ineg = 1;
        }

        Float_t cheat_gen_angle = abs(pfot.cheat_PiLPFO[ineg].cos) * sgn( -_mc.mc_quark_charge[1] ) * mct.mc_quark[1].cos / abs(mct.mc_quark[1].cos);
        _hm.h1[_hm.cheat_Pi_cos]->Fill(pfot.cheat_PiLPFO[ineg].cos);
        _hm.h1[_hm.cheat_Pi_qcos]->Fill(cheat_gen_angle);

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
    Float_t cosacol = VectorTools::GetCosBetween( v[0].v3(), v[1].v3() );

    return (cosacol < -0.95);
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

  Int_t *EventAnalyzer::Gen_Reco_Stats_Stable( PFOTools mct, PFOTools pfot, SelectID pid, Float_t cos_min, Float_t cos_max )
  {
    std::vector<PFO_Info> PFO_Collection = pfot.Valid_PFOs;

    Float_t p_min = _anCfg.PFO_p_min;

    std::vector<PFO_Info> PFO_Hadron_Collection;
    std::vector<PFO_Info> Gen_Hadron_Collection;
    for ( auto iPFO : PFO_Collection ){
      Bool_t cos_range = (cos_min < iPFO.cos && iPFO.cos < cos_max );
      Bool_t p_range   = p_min < iPFO.p_mag;
      Bool_t is_cos_p  = cos_range && p_range;
      
      if ( !is_cos_p ) continue;

      switch ( pid ) {
        case kKaon:
          if( PFOTools::isKaon(iPFO) )        PFO_Hadron_Collection.push_back(iPFO);
          if( abs(iPFO.pfo_pdgcheat) == 321 ) Gen_Hadron_Collection.push_back(iPFO);
          break;

        case kPion:
          if( PFOTools::isPion(iPFO) )        PFO_Hadron_Collection.push_back(iPFO);
          if( abs(iPFO.pfo_pdgcheat) == 211 ) Gen_Hadron_Collection.push_back(iPFO);
          break;

        default:
          break;
      }

    }

    Int_t N_K_corr  = 0;

    // Float_t cos_r = 0.02;
    Float_t cos_r = 0.37;
    std::vector<PFO_Info> PFO_Hadron_Remain = PFO_Hadron_Collection;

    for ( auto igen : Gen_Hadron_Collection ){
      Float_t min_cos_diff   = 1000.0;
      Int_t   i_min_cos_diff = -1;

      Int_t counter = 0;
      for ( auto iremain : PFO_Hadron_Remain ){
        Float_t cos_diff = 1.0 - VectorTools::GetCosBetween( igen.vt.v3(), iremain.vt.v3() );

        if( cos_diff < min_cos_diff ) {
          min_cos_diff = cos_diff;
          i_min_cos_diff = counter;
        }
        counter++;
      }

      if( min_cos_diff < cos_r ) {
        N_K_corr++;
        PFO_Hadron_Remain.erase( PFO_Hadron_Remain.begin() + i_min_cos_diff );
      }

    }

    static Int_t N_array[3] = {0};
    N_array[0] = Gen_Hadron_Collection.size();
    N_array[1] = PFO_Hadron_Collection.size();
    N_array[2] = N_K_corr;

    return N_array;

  }

  Int_t *EventAnalyzer::Gen_Reco_Stats_Cheat( PFOTools mct, PFOTools pfot, SelectID pid, Float_t cos_min, Float_t cos_max )
  {
    std::vector<PFO_Info> PFO_Collection;
    std::vector<PFO_Info> jet[2] = { pfot.GetJet(0), pfot.GetJet(1) };

    PFO_Collection.reserve( jet[0].size() + jet[1].size() );
    PFO_Collection.insert( PFO_Collection.begin(), jet[0].begin(), jet[0].end() );
    PFO_Collection.insert( PFO_Collection.end(), jet[1].begin(), jet[1].end() );

    Float_t p_min = _anCfg.PFO_p_min;

    std::vector<PFO_Info> PFO_Hadron_Collection;
    for ( auto iPFO : PFO_Collection ){
      Bool_t cos_range = (cos_min < iPFO.cos && iPFO.cos < cos_max );
      Bool_t p_range   = p_min < iPFO.p_mag;
      if( PFOTools::isKaon(iPFO) && cos_range && p_range ) PFO_Hadron_Collection.push_back(iPFO);
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
    std::vector<PFO_Info> PFO_Hadron_Remain = PFO_Hadron_Collection;

    for ( auto igen : PFO_Cheat_K_Collection ){
      Float_t min_cos_diff   = 1000.0;
      Int_t   i_min_cos_diff = -1;

      Int_t counter = 0;
      for ( auto iremain : PFO_Hadron_Remain ){
        Float_t cos_diff = 1.0 - VectorTools::GetCosBetween( igen.vt.v3(), iremain.vt.v3() );
        if( cos_diff < min_cos_diff ) {
          min_cos_diff = cos_diff;
          i_min_cos_diff = counter;
        }
        counter++;
      }

      if( min_cos_diff < cos_r ) {
        N_K_corr++;
        PFO_Hadron_Remain.erase( PFO_Hadron_Remain.begin() + i_min_cos_diff );
      }

    }

    static Int_t N_array[3] = {0};
    N_array[0] = PFO_Cheat_K_Collection.size();
    N_array[1] = PFO_Hadron_Collection.size();
    N_array[2] = N_K_corr;

    return N_array;

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
    _hm.h1[_hm.gen_q_cos]->Fill(mct.mc_quark[0].cos);
    _hm.h1[_hm.gen_q_qcos]->Fill(mct.mc_quark[0].qcos);

    // Gen K
    for ( int istable=0; istable < _mc.mc_stable_n; istable++ ){
      if(abs(_mc.mc_stable_pdg[istable]) == 321) {
        _hm.h1[_hm.gen_K_cos]->Fill(mct.mc_stable[istable].cos);
        _hm.h1[_hm.gen_K_qcos]->Fill(mct.mc_stable[istable].qcos);
      }
    }

  }

  void EventAnalyzer::ProcessDoubleTag(PFOTools pfot, PFOTools mct, unordered_map< TString, unordered_map<TString, Bool_t> > cuts)
  {
    unordered_map< TString, vector<Bool_t> > is_pass;
    for( const auto &[i_lmode, cuts_for_lmode] : cuts ){
      for ( const auto &[cut_name, cut_val] : cuts_for_lmode ) {
        if( cut_name != "charge" ) is_pass[i_lmode].push_back(cut_val);
      }
    }

    for( const auto &[i_lmode, val_LPFO] : pfot.LPFO_ ){

      Bool_t isPass = std::all_of(is_pass.at(i_lmode).begin(), is_pass.at(i_lmode).end(), [](bool v) { return v; });

      if( isPass ){

        Int_t ineg             = val_LPFO.at(0).pfo_charge > 0 ;
        PFO_Info LPFO          = val_LPFO.at(ineg);
        PFO_Info LPFO_opposite = val_LPFO.at(1-ineg);

        if( cuts.at(i_lmode).at("charge") ){

          _hm.h1_cos.at(i_lmode).at("cos")->Fill( LPFO.cos );
          _hm.h1_cos.at(i_lmode).at("qcos")->Fill( LPFO.qcos );
          _hm.h1_cos.at(i_lmode).at("acc_cos")->Fill( LPFO.qcos );

          Float_t scos = abs(LPFO.cos) * sgn( -_mc.mc_quark_charge[0] ) * mct.mc_quark[0].cos / abs(mct.mc_quark[0].cos);
          _hm.h1_cos.at(i_lmode).at("scos")->Fill( scos );

          auto it = _pt.PFO_type_map.find(LPFO.pfo_pdgcheat);
          if( it != _pt.PFO_type_map.end() ){
            TString type = it->second;
            _hm.h2_dEdx.at(i_lmode).at(type).at("dEdx_p")->Fill( LPFO.p_mag, LPFO.pfo_dedx );
            _hm.h2_dEdx.at(i_lmode).at(type).at("dEdx_dist_cos")->Fill( LPFO.p_mag, pfot.Get_dEdx_dist(LPFO, i_lmode) );
          }

        }else{

          _hm.h1_cos.at(i_lmode).at("rej_cos")->Fill( LPFO.cos );
          _hm.h1_cos.at(i_lmode).at("rej_cos")->Fill( -LPFO_opposite.cos );

        } // charge consistency check

      } // LPFO event selection

    } // LPFO mode loop

  }

  void EventAnalyzer::ProcessDoubleTag(PFOTools pfot, PFOTools mct, vector<Bool_t> cuts[3], PDGConfig double_tag[3])
  {
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

    // Reco K_K

    if ( LPFO_checks[kKaon] && pfot.is_ss() && double_tag[kKaon] == K_K ){

      Int_t ineg = -1;

      if( pfot.KLPFO[0].pfo_charge < 0 ){
        ineg = 0;
      }else{
        ineg = 1;
      }

      if(sign_check[kKaon]){

        Float_t gen_reco_K_sep_cos  = VectorTools::GetCosBetween(pfot.KLPFO[ineg].vt.v3(), mct.mc_quark[0].vt.v3());

        _hm.h1[_hm.reco_K_cos]->Fill( pfot.KLPFO[ineg].cos );
        _hm.h1[_hm.reco_K_qcos]->Fill( pfot.KLPFO[ineg].qcos );
        _hm.h1[_hm.reco_K_scos]->Fill( abs(pfot.KLPFO[ineg].cos) * sgn( -_mc.mc_quark_charge[0] ) * mct.mc_quark[0].cos / abs(mct.mc_quark[0].cos) );
        _hm.h1[_hm.reco_K_mom]->Fill( pfot.KLPFO[ineg].p_mag );
        _hm.h1[_hm.gen_reco_K_sep_cos]->Fill( gen_reco_K_sep_cos );

        // cheat
        switch ( abs(pfot.KLPFO[ineg].pfo_pdgcheat) ) {
          case 321:
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_K_dEdx_p][kKaon]->Fill(pfot.KLPFO[ineg].p_mag,pfot.KLPFO[ineg].pfo_dedx);
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_K_KdEdx_dist_cos][kKaon]->Fill(pfot.KLPFO[ineg].cos,pfot.KLPFO[ineg].pfo_piddedx_k_dedxdist);
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 1 );
            break;
          case 211:
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_K_dEdx_p][kPion]->Fill(pfot.KLPFO[ineg].p_mag,pfot.KLPFO[ineg].pfo_dedx);
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_K_KdEdx_dist_cos][kPion]->Fill(pfot.KLPFO[ineg].cos,pfot.KLPFO[ineg].pfo_piddedx_k_dedxdist);
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 0 );
            break;
          case 2212:
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_K_dEdx_p][kProton]->Fill(pfot.KLPFO[ineg].p_mag,pfot.KLPFO[ineg].pfo_dedx);
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_K_KdEdx_dist_cos][kProton]->Fill(pfot.KLPFO[ineg].cos,pfot.KLPFO[ineg].pfo_piddedx_k_dedxdist);
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 2 );
            break;
          case 11:
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_K_dEdx_p][kElectron]->Fill(pfot.KLPFO[ineg].p_mag,pfot.KLPFO[ineg].pfo_dedx);
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_K_KdEdx_dist_cos][kElectron]->Fill(pfot.KLPFO[ineg].cos,pfot.KLPFO[ineg].pfo_piddedx_k_dedxdist);
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 3 );
            break;
          case 13:
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_K_dEdx_p][kMuon]->Fill(pfot.KLPFO[ineg].p_mag,pfot.KLPFO[ineg].pfo_dedx);
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_K_KdEdx_dist_cos][kMuon]->Fill(pfot.KLPFO[ineg].cos,pfot.KLPFO[ineg].pfo_piddedx_k_dedxdist);
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 4 );
            break;
          default:
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 5 );
            break;
        }

        switch ( abs(pfot.KLPFO[1-ineg].pfo_pdgcheat) ) {
          case 321:
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 1 );
            break;
          case 211:
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 0 );
            break;
          case 2212:
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 2 );
            break;
          case 11:
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 3 );
            break;
          case 13:
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 4 );
            break;
          default:
            _hm.h1[_hm.reco_K_pdgcheat]->Fill( 5 );
            break;
        }

        _hm.h1_pq[_hm.acc_KK]->Fill( pfot.KLPFO[ineg].qcos );

      }else{

        _hm.h1_pq[_hm.rej_KK]->Fill( pfot.KLPFO[ineg].cos );
        _hm.h1_pq[_hm.rej_KK]->Fill( -pfot.KLPFO[1-ineg].cos );
      
      }

    }

    // case Pi_Pi:
    if ( LPFO_checks[kPion] && pfot.is_uu_dd() && double_tag[kPion] == Pi_Pi ){

      Int_t ineg = -1;

      if( pfot.PiLPFO[0].pfo_charge < 0 ){
        ineg = 0;
      }else{
        ineg = 1;
      }

      if(sign_check[kPion]){

        Float_t gen_reco_Pi_sep_cos  = VectorTools::GetCosBetween(pfot.PiLPFO[ineg].vt.v3(), mct.mc_quark[0].vt.v3());
        Float_t vtx_endpt = sqrt(pfot.PiLPFO[ineg].pfo_endpt[0] * pfot.PiLPFO[ineg].pfo_endpt[0] + pfot.PiLPFO[ineg].pfo_endpt[1] * pfot.PiLPFO[ineg].pfo_endpt[1] + pfot.PiLPFO[ineg].pfo_endpt[2] * pfot.PiLPFO[ineg].pfo_endpt[2]);

        if( gen_reco_Pi_sep_cos > 0 ){
          _hm.h1[_hm.good_reco_Pi_endpt]->Fill( vtx_endpt );
          _hm.h1[_hm.good_reco_Pi_tpchits]->Fill( pfot.PiLPFO[ineg].pfo_tpc_hits );
          _hm.h1[_hm.good_reco_Pi_pidedx_dist]->Fill( pfot.PiLPFO[ineg].pfo_piddedx_pi_dedxdist );
          _hm.h1[_hm.good_reco_Pi_kdedx_dist]->Fill( pfot.PiLPFO[ineg].pfo_piddedx_k_dedxdist );
        }else{
          _hm.h1[_hm.bad_reco_Pi_endpt]->Fill( vtx_endpt );
          _hm.h1[_hm.bad_reco_Pi_tpchits]->Fill( pfot.PiLPFO[ineg].pfo_tpc_hits );
          _hm.h1[_hm.bad_reco_Pi_pidedx_dist]->Fill( pfot.PiLPFO[ineg].pfo_piddedx_pi_dedxdist );
          _hm.h1[_hm.bad_reco_Pi_kdedx_dist]->Fill( pfot.PiLPFO[ineg].pfo_piddedx_k_dedxdist );
        }

        _hm.h1[_hm.reco_Pi_cos]->Fill( pfot.PiLPFO[ineg].cos );
        _hm.h1[_hm.reco_Pi_qcos]->Fill( pfot.PiLPFO[ineg].qcos );
        _hm.h1[_hm.reco_Pi_scos]->Fill( abs(pfot.PiLPFO[ineg].cos) * sgn( -_mc.mc_quark_charge[1] ) * mct.mc_quark[1].cos / abs(mct.mc_quark[1].cos) );
        _hm.h1[_hm.reco_Pi_mom]->Fill( pfot.PiLPFO[ineg].p_mag );
        _hm.h1[_hm.gen_reco_Pi_sep_cos]->Fill( gen_reco_Pi_sep_cos );

        // cheat
        switch ( abs(pfot.PiLPFO[ineg].pfo_pdgcheat) ) {
          case 321:
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_Pi_dEdx_p][kKaon]->Fill(pfot.PiLPFO[ineg].p_mag,pfot.PiLPFO[ineg].pfo_dedx);
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_Pi_PidEdx_dist_cos][kKaon]->Fill(pfot.PiLPFO[ineg].cos,pfot.PiLPFO[ineg].pfo_piddedx_pi_dedxdist);
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 1 );
            break;
          case 211:
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_Pi_dEdx_p][kPion]->Fill(pfot.PiLPFO[ineg].p_mag,pfot.PiLPFO[ineg].pfo_dedx);
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_Pi_PidEdx_dist_cos][kPion]->Fill(pfot.PiLPFO[ineg].cos,pfot.PiLPFO[ineg].pfo_piddedx_pi_dedxdist);
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 0 );
            break;
          case 2212:
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_Pi_dEdx_p][kProton]->Fill(pfot.PiLPFO[ineg].p_mag,pfot.PiLPFO[ineg].pfo_dedx);
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_Pi_PidEdx_dist_cos][kProton]->Fill(pfot.PiLPFO[ineg].cos,pfot.KLPFO[ineg].pfo_piddedx_pi_dedxdist);
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 2 );
            break;
          case 11:
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_Pi_dEdx_p][kElectron]->Fill(pfot.PiLPFO[ineg].p_mag,pfot.PiLPFO[ineg].pfo_dedx);
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_Pi_PidEdx_dist_cos][kElectron]->Fill(pfot.PiLPFO[ineg].cos,pfot.PiLPFO[ineg].pfo_piddedx_pi_dedxdist);
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 3 );
            break;
          case 13:
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_Pi_dEdx_p][kMuon]->Fill(pfot.PiLPFO[ineg].p_mag,pfot.PiLPFO[ineg].pfo_dedx);
            _hm.h2_dEdx_nomap[_hm.gen_ipart_reco_Pi_PidEdx_dist_cos][kMuon]->Fill(pfot.PiLPFO[ineg].cos,pfot.PiLPFO[ineg].pfo_piddedx_pi_dedxdist);
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 4 );
            break;
          default:
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 5 );
            break;
        }

        switch ( abs(pfot.PiLPFO[1-ineg].pfo_pdgcheat) ) {
          case 321:
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 1 );
            break;
          case 211:
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 0 );
            break;
          case 2212:
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 2 );
            break;
          case 11:
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 3 );
            break;
          case 13:
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 4 );
            break;
          default:
            _hm.h1[_hm.reco_Pi_pdgcheat]->Fill( 5 );
            break;
        }

        _hm.h1_pq[_hm.acc_PiPi]->Fill( pfot.PiLPFO[ineg].qcos );

      }else{

        _hm.h1_pq[_hm.rej_PiPi]->Fill( pfot.PiLPFO[ineg].cos );
        _hm.h1_pq[_hm.rej_PiPi]->Fill( -pfot.PiLPFO[1-ineg].cos );
      
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
    Float_t cosacol = VectorTools::GetCosBetween( jetvt[0].v3(), jetvt[1].v3() );
    Float_t jet_cos[2]  = { std::cos( jetvt[0].v3().theta() ), std::cos( jetvt[1].v3().theta() ) };

    _data.sum_jet_E = _jet.jet_E[0] + _jet.jet_E[1];
    _data.jet_acol  = cosacol;
    _data.jet_theta_diff = std::abs( jetvt[0].v3().theta() - jetvt[1].v3().theta() );

    _hm.h1[_hm.reco_sum_jetE]->Fill( _data.sum_jet_E );
    _hm.h1[_hm.reco_jet_sep]->Fill( _data.jet_acol );

    _hm.h2_jet[_hm.jet_mult_cos]->Fill( jet_cos[0], _jet.jet_npfo[0] );
    _hm.h2_jet[_hm.jet_mult_cos]->Fill( jet_cos[1], _jet.jet_npfo[1] );

  }
}