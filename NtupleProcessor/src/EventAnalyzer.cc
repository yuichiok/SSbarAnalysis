
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

    ////////////////
    // Selections //
    ////////////////

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

    //Double Tagging
    ProcessDoubleTag(pfot,mct,CutTriggerMap);

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

  }
}