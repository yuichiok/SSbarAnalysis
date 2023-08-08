
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

#define DEBUG_PRINT std::cout << "debug: " << _check_pt++ << std::endl;

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

  void EventAnalyzer::AnalyzeGen(Long64_t entry)
  {
    ientry = entry;
    // cout << "evt: " << entry << endl;

    PFOTools mct( &_mc, _config );
    PFOTools pfot( &_mc, &_pfo, _config );

    // Gen QQbar
    _hm.h1_gen_cos.at("cos")->Fill(mct.mc_quark[0].cos);
    _hm.h1_gen_cos.at("cos")->Fill(mct.mc_quark[1].cos);

    _hm.h1_gen_cos.at("qcos")->Fill(mct.mc_quark[0].qcos);

    Bool_t is_noJet = pfot.PFO_sorted_jet[0].size() == 0 || pfot.PFO_sorted_jet[1].size() == 0;
    if( is_noJet ) return;
    unordered_map< int, vector<PFO_Info> > hadronJet;
    for ( int ijet = 0; ijet < 2; ijet++ ){
      for ( const auto& ipfo : pfot.PFO_sorted_jet[ijet] ){
        if( abs(ipfo.pfo_pdgcheat) != 11 && abs(ipfo.pfo_pdgcheat) != 13 ){
          hadronJet[ijet].push_back(ipfo);
        }
      }
    }
    Bool_t is_LPFO_hadron = hadronJet[0].size() == 0 || hadronJet[1].size() == 0;
    if( is_LPFO_hadron ) return;

    unordered_map< TString, unordered_map<TString, Bool_t> > CutTriggerMap; // [particle][cutname]
    for ( const auto i_lmode : _pt.PFO_mode ){

      CutTriggerMap[i_lmode] = TriggerMap( pfot, i_lmode, hadronJet, "gen" );

    }

    //Double Tagging Efficiency
    ProcessDoubleTagEfficiency(pfot,mct,CutTriggerMap,hadronJet,"gen");

    ClearStructs();

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

    Bool_t is_noJet = pfot.PFO_sorted_jet[0].size() == 0 || pfot.PFO_sorted_jet[1].size() == 0;
    if( is_noJet ) return;

    unordered_map< int, vector<PFO_Info> > hadronJet;
    for ( int ijet = 0; ijet < 2; ijet++ ){
      for ( const auto& ipfo : pfot.PFO_sorted_jet[ijet] ){
        if( abs(ipfo.pfo_pid) != 11 && abs(ipfo.pfo_pid) != 13 ){
          hadronJet[ijet].push_back(ipfo);
        }
      }
    }
    Bool_t is_LPFO_hadron = hadronJet[0].size() == 0 || hadronJet[1].size() == 0;
    if( is_LPFO_hadron ) return;

    unordered_map< TString, unordered_map<TString, Bool_t> > CutTriggerMap; // [particle][cutname]
    for ( const auto i_lmode : _pt.PFO_mode ){

      // CutTriggerMap[i_lmode] = TriggerMap( pfot, i_lmode, pfot.PFO_sorted_jet, "reco" );
      CutTriggerMap[i_lmode] = TriggerMap( pfot, i_lmode, hadronJet, "reco" );

    }

    //Double Tagging Efficiency
    ProcessDoubleTagEfficiency(pfot,mct,CutTriggerMap,hadronJet,"reco");

    //Double Tagging
    ProcessDoubleTag(pfot,mct,CutTriggerMap,hadronJet);

    // Resolution Analysis
    ResolutionAnalysis(pfot,mct);

    ClearStructs();

  }

  void EventAnalyzer::ClearStructs()
  {
    _check_pt = 0;
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

  unordered_map<TString, Bool_t> EventAnalyzer::TriggerMap( PFOTools pfot, TString i_lmode, unordered_map< int, vector<PFO_Info> > subjet_pair, TString gen_reco )
  {
    unordered_map<TString, Bool_t> outMap;
    vector<PFO_Info> LPFOs = {subjet_pair.at(0).at(0), subjet_pair.at(1).at(0)};

    // check momentum
    outMap["momentum"] = pfot.is_momentum( LPFOs.at(0), _anCfg.LPFO_p_min, _anCfg.LPFO_p_max ) &&
                         pfot.is_momentum( LPFOs.at(1), _anCfg.LPFO_p_min, _anCfg.LPFO_p_max );

    // check tpc hits
    outMap["tpc_hits"] = pfot.is_tpc_hits( LPFOs.at(0), _anCfg.PFO_TPCHits_min ) &&
                         pfot.is_tpc_hits( LPFOs.at(1), _anCfg.PFO_TPCHits_min );

    // check offsets
    outMap["offset"] = pfot.is_offset_small( LPFOs.at(0), _anCfg.PFO_offset_max ) &&
                       pfot.is_offset_small( LPFOs.at(1), _anCfg.PFO_offset_max );

    // check dEdx dist
    if( gen_reco == "reco" ){
      outMap["PID"] = pfot.is_PID( i_lmode, LPFOs.at(0) ) &&
                      pfot.is_PID( i_lmode, LPFOs.at(1) );
    }else{
      outMap["PID"] = ( abs(LPFOs.at(0).pfo_pdgcheat) == pfot.PFO_type_map_rev.at(i_lmode) ) &&
                      ( abs(LPFOs.at(1).pfo_pdgcheat) == pfot.PFO_type_map_rev.at(i_lmode) );
    }
    
    // SPFO opposite check
    vector<Bool_t> is_SPFO_charge_opposite(2,false);
    for ( int ijet=0; ijet<2; ijet++ ){
      vector<PFO_Info> subjet = subjet_pair.at(ijet);
      if( subjet.size() == 0 ) continue;
      for ( auto iSPFO : subjet ){
        if(iSPFO.ipfo == LPFOs.at(ijet).ipfo) continue;
        Bool_t charge_opposite = iSPFO.pfo_charge * LPFOs.at(ijet).pfo_charge < 0;
        Bool_t momentum_above  = iSPFO.p_mag > 10;
        if( charge_opposite && momentum_above ) is_SPFO_charge_opposite.at(ijet) = true;
      }
    }
    outMap["SPFO"] = std::none_of(is_SPFO_charge_opposite.begin(), is_SPFO_charge_opposite.end(), [](bool v) { return v; });

    // Charge opposite check
    outMap["charge"] = pfot.is_charge_config(pfot.kOpposite,LPFOs.at(0).pfo_charge,LPFOs.at(1).pfo_charge);

    return outMap;

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

  void EventAnalyzer::ResolutionAnalysis( PFOTools pfot, PFOTools mct )
  {
    for( auto i_lmode : pfot.PFO_mode ){
      Int_t nbins_cos = _hm.nbins_cos;
      TAxis *xaxis    = _hm.h1_resolution.at(i_lmode).at("gen_N_cos")->GetXaxis();
      for ( int ibin=1; ibin<=nbins_cos; ibin++ ){
        Float_t bin_center = xaxis->GetBinCenter(ibin);
        Float_t bin_width  = xaxis->GetBinWidth(ibin);
        Float_t cos_min    = xaxis->GetBinLowEdge(ibin);
        Float_t cos_max    = cos_min + bin_width;
        unordered_map<TString, Int_t> N_particles = Gen_Reco_Stats_Stable( pfot, mct, i_lmode, cos_min, cos_max );

        for( auto iname : _hm.hres_name ){
          _hm.h1_resolution.at(i_lmode).at(iname)->Fill( bin_center ,N_particles.at(iname));
        }
      }
    }
  }

  unordered_map<TString, Int_t> EventAnalyzer::Gen_Reco_Stats_Stable( PFOTools pfot, PFOTools mct, TString lmode, Float_t cos_min, Float_t cos_max )
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

      if( pfot.is_PID(lmode,iPFO) )                                   PFO_Hadron_Collection.push_back(iPFO);
      if( abs(iPFO.pfo_pdgcheat) == pfot.PFO_type_map_rev.at(lmode) ) Gen_Hadron_Collection.push_back(iPFO);


    }

    Int_t N_K_corr  = 0;

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


    unordered_map<TString, Int_t> N_map;
    N_map[_hm.hres_name.at(0)] = Gen_Hadron_Collection.size();
    N_map[_hm.hres_name.at(1)] = PFO_Hadron_Collection.size();
    N_map[_hm.hres_name.at(2)] = N_K_corr;

    return N_map;

  }

  Bool_t EventAnalyzer::is_high_LPFO( PFOTools pfot, TString lmode )
  {
    TString other_mode = (lmode == "K") ? "Pi" : "K";

    Int_t check = 0;
    for ( int ijet=0; ijet<2; ijet++ ){
      vector<PFO_Info> subjet_lmode = pfot.PFO_subjet.at(lmode).at(ijet);
      vector<PFO_Info> subjet_other = pfot.PFO_subjet.at(other_mode).at(ijet);
      if( subjet_lmode.size() == 0 || subjet_other.size() == 0 ){
        check++;
        continue;
      }

      PFO_Info LPFO_lmode = subjet_lmode.at(0);
      PFO_Info LPFO_other = subjet_other.at(0);
      if( LPFO_lmode.p_mag > LPFO_other.p_mag ) check++;
    }

    return check == 2;

  }

  void EventAnalyzer::ProcessDoubleTag(PFOTools pfot, PFOTools mct, unordered_map< TString, unordered_map<TString, Bool_t> > cuts, unordered_map< int, vector<PFO_Info> > subjet_pair)
  {
    unordered_map< TString, vector<Bool_t> > is_pass;
    for( const auto &[i_lmode, cuts_for_lmode] : cuts ){
      for ( const auto &[cut_name, cut_val] : cuts_for_lmode ) {
        if( cut_name != "charge" ) is_pass[i_lmode].push_back(cut_val);
      }
    }

    // vector<PFO_Info> LPFOs = pfot.LPFO;
    for ( const auto &[i_lmode, cut_vec] : is_pass ){

      Bool_t isPass = std::all_of(cut_vec.begin(), cut_vec.end(), [](bool v) { return v; });

      if( isPass ){

        vector<PFO_Info> LPFOs = {subjet_pair.at(0).at(0), subjet_pair.at(1).at(0)};

        Int_t ineg             = LPFOs.at(0).pfo_charge > 0 ;
        PFO_Info LPFO          = LPFOs.at(ineg);
        PFO_Info LPFO_opposite = LPFOs.at(1-ineg);

        if( cuts.at(i_lmode).at("charge") ){

          _hm.h1_cos.at(i_lmode).at("cos")->Fill( LPFO.cos );
          _hm.h1_cos.at(i_lmode).at("cos")->Fill( LPFO_opposite.cos );

          _hm.h1_cos.at(i_lmode).at("qcos")->Fill( LPFO.qcos );
          _hm.h1_cos.at(i_lmode).at("acc_cos")->Fill( LPFO.qcos );

          Float_t scos = abs(LPFO.cos) * sgn( -_mc.mc_quark_charge[0] ) * mct.mc_quark[0].cos / abs(mct.mc_quark[0].cos);
          _hm.h1_cos.at(i_lmode).at("scos")->Fill( scos );

          if( _pt.PFO_type_map.find(abs(LPFO.pfo_pdgcheat)) != _pt.PFO_type_map.end() ){

            TString type = _pt.PFO_type_map.at(abs(LPFO.pfo_pdgcheat));

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

  void EventAnalyzer::ProcessDoubleTagEfficiency(PFOTools pfot, PFOTools mct, unordered_map< TString, unordered_map<TString, Bool_t> > cuts, unordered_map< int, vector<PFO_Info> > subjet_pair, TString gen_reco)
  {
    Bool_t isLPFO = subjet_pair[0].size() == 0 || subjet_pair[1].size() == 0;
    if( isLPFO ) return;
    vector<PFO_Info> LPFOs = {subjet_pair.at(0).at(0), subjet_pair.at(1).at(0)};

    for( const auto &[i_lmode, cuts_for_lmode] : cuts ){

      Bool_t selection = true;
      for( auto icut_name : _hm.heff_name ){
        selection = selection && cuts_for_lmode.at(icut_name);
        if(selection) {
          _hm.h1_cos_eff.at(gen_reco).at(i_lmode).at(icut_name)->Fill( LPFOs.at(0).cos );
          _hm.h1_cos_eff.at(gen_reco).at(i_lmode).at(icut_name)->Fill( LPFOs.at(1).cos );
        }
      }
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

  }
}