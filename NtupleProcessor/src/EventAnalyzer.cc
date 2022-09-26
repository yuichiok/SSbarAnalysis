
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
#include "TreeWriter.hh"

using std::cout;   using std::endl;
typedef unsigned int Index;

EventAnalyzer::EventAnalyzer(TString o)
: options(o)
{
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

void EventAnalyzer::InitWriteTree()
{

  // Initialize Write Tree
    TreeWriter writer;
    _hfile = new TFile( _hfilename, "RECREATE", _hfilename ) ;
    
    _hTree_LPFO     = new TTree( "LPFO", "tree" );
    _hTree_LPFO_KK  = new TTree( "LPFO_KK", "tree" );
    _hTree_LPFO_KPi = new TTree( "LPFO_KPi", "tree" );

    writer.InitializeLPFOTree(_hTree_LPFO, _tree_lpfo);
    writer.InitializeLPFOTree(_hTree_LPFO_KK, _tree_lpfo_kk);
    writer.InitializeLPFOTree(_hTree_LPFO_KPi, _tree_lpfo_kpi);

}

void EventAnalyzer::WriteFile()
{

  // Write Tree
    _hfile->Write();
    _hfile->Close();

}

void EventAnalyzer::WriteLPFOVariables(PFOTools pt, Tree_SSbar *data)
{
  Int_t iLPFOs[2] = {pt.LPFO[0].ipfo, pt.LPFO[1].ipfo};
  int i = 0;
  for (auto iLPFO : iLPFOs){
    data->lpfo_match                       [i] = _pfo.pfo_match                       [iLPFO];
    data->lpfo_truejet_pdg                 [i] = _pfo.pfo_truejet_pdg                 [iLPFO];
    data->lpfo_truejet_type                [i] = _pfo.pfo_truejet_type                [iLPFO];
    data->lpfo_pdgcheat                    [i] = _pfo.pfo_pdgcheat                    [iLPFO];
    data->lpfo_nparents                    [i] = _pfo.pfo_nparents                    [iLPFO];
    // data->lpfo_pdgcheat_parent             [i] = _pfo.pfo_pdgcheat_parent             [iLPFO];
    data->lpfo_E                           [i] = _pfo.pfo_E                           [iLPFO];
    data->lpfo_px                          [i] = _pfo.pfo_px                          [iLPFO];
    data->lpfo_py                          [i] = _pfo.pfo_py                          [iLPFO];
    data->lpfo_pz                          [i] = _pfo.pfo_pz                          [iLPFO];
    data->lpfo_m                           [i] = _pfo.pfo_m                           [iLPFO];
    data->lpfo_type                        [i] = _pfo.pfo_type                        [iLPFO];
    data->lpfo_isoverlay                   [i] = _pfo.pfo_isoverlay                   [iLPFO];
    data->lpfo_isisr                       [i] = _pfo.pfo_isisr                       [iLPFO];
    data->lpfo_vtx                         [i] = _pfo.pfo_vtx                         [iLPFO];
    data->lpfo_charge                      [i] = _pfo.pfo_charge                      [iLPFO];
    data->lpfo_ntracks                     [i] = _pfo.pfo_ntracks                     [iLPFO];
    data->lpfo_tpc_hits                    [i] = _pfo.pfo_tpc_hits                    [iLPFO];
    data->lpfo_dedx                        [i] = _pfo.pfo_dedx                        [iLPFO];
    data->lpfo_dedxerror                   [i] = _pfo.pfo_dedxerror                   [iLPFO];
    data->lpfo_d0                          [i] = _pfo.pfo_d0                          [iLPFO];
    data->lpfo_d0error                     [i] = _pfo.pfo_d0error                     [iLPFO];
    data->lpfo_z0                          [i] = _pfo.pfo_z0                          [iLPFO];
    data->lpfo_z0error                     [i] = _pfo.pfo_z0error                     [iLPFO];
    data->lpfo_phi                         [i] = _pfo.pfo_phi                         [iLPFO];
    data->lpfo_phierror                    [i] = _pfo.pfo_phierror                    [iLPFO];
    data->lpfo_omega                       [i] = _pfo.pfo_omega                       [iLPFO];
    data->lpfo_omegaerror                  [i] = _pfo.pfo_omegaerror                  [iLPFO];
    data->lpfo_tanlambda                   [i] = _pfo.pfo_tanlambda                   [iLPFO];
    data->lpfo_tanlambdaerror              [i] = _pfo.pfo_tanlambdaerror              [iLPFO];
    data->lpfo_chi2                        [i] = _pfo.pfo_chi2                        [iLPFO];
    data->lpfo_ndf                         [i] = _pfo.pfo_ndf                         [iLPFO];
    // data->lpfo_vtxpt                       [i] = _pfo.pfo_vtxpt                       [iLPFO];
    // data->lpfo_endpt                       [i] = _pfo.pfo_endpt                       [iLPFO];
    data->lpfo_pid                         [i] = _pfo.pfo_pid                         [iLPFO];
    data->lpfo_pid_likelihood              [i] = _pfo.pfo_pid_likelihood              [iLPFO];
    data->lpfo_pid_eprob                   [i] = _pfo.pfo_pid_eprob                   [iLPFO];
    data->lpfo_pid_muprob                  [i] = _pfo.pfo_pid_muprob                  [iLPFO];
    data->lpfo_pid_piprob                  [i] = _pfo.pfo_pid_piprob                  [iLPFO];
    data->lpfo_pid_kprob                   [i] = _pfo.pfo_pid_kprob                   [iLPFO];
    data->lpfo_pid_pprob                   [i] = _pfo.pfo_pid_pprob                   [iLPFO];
    data->lpfo_pid_hprob                   [i] = _pfo.pfo_pid_hprob                   [iLPFO];
    data->lpfo_piddedx                     [i] = _pfo.pfo_piddedx                     [iLPFO];
    data->lpfo_piddedx_likelihood          [i] = _pfo.pfo_piddedx_likelihood          [iLPFO];
    data->lpfo_piddedx_eprob               [i] = _pfo.pfo_piddedx_eprob               [iLPFO];
    data->lpfo_piddedx_muprob              [i] = _pfo.pfo_piddedx_muprob              [iLPFO];
    data->lpfo_piddedx_piprob              [i] = _pfo.pfo_piddedx_piprob              [iLPFO];
    data->lpfo_piddedx_kprob               [i] = _pfo.pfo_piddedx_kprob               [iLPFO];
    data->lpfo_piddedx_pprob               [i] = _pfo.pfo_piddedx_pprob               [iLPFO];
    data->lpfo_piddedx_hprob               [i] = _pfo.pfo_piddedx_hprob               [iLPFO];
    data->lpfo_piddedx_e_dedxdist          [i] = _pfo.pfo_piddedx_e_dedxdist          [iLPFO];
    data->lpfo_piddedx_mu_dedxdist         [i] = _pfo.pfo_piddedx_mu_dedxdist         [iLPFO];
    data->lpfo_piddedx_pi_dedxdist         [i] = _pfo.pfo_piddedx_pi_dedxdist         [iLPFO];
    data->lpfo_piddedx_k_dedxdist          [i] = _pfo.pfo_piddedx_k_dedxdist          [iLPFO];
    data->lpfo_piddedx_p_dedxdist          [i] = _pfo.pfo_piddedx_p_dedxdist          [iLPFO];
    data->lpfo_piddedx_e_lkhood            [i] = _pfo.pfo_piddedx_e_lkhood            [iLPFO];
    data->lpfo_piddedx_mu_lkhood           [i] = _pfo.pfo_piddedx_mu_lkhood           [iLPFO];
    data->lpfo_piddedx_pi_lkhood           [i] = _pfo.pfo_piddedx_pi_lkhood           [iLPFO];
    data->lpfo_piddedx_k_lkhood            [i] = _pfo.pfo_piddedx_k_lkhood            [iLPFO];
    data->lpfo_piddedx_p_lkhood            [i] = _pfo.pfo_piddedx_p_lkhood            [iLPFO];
    data->lpfo_pidtof_p_at_calo            [i] = _pfo.pfo_pidtof_p_at_calo            [iLPFO];
    data->lpfo_pidtof_closest_beta_0ps     [i] = _pfo.pfo_pidtof_closest_beta_0ps     [iLPFO];
    data->lpfo_pidtof_closest_beta_10ps    [i] = _pfo.pfo_pidtof_closest_beta_10ps    [iLPFO];
    data->lpfo_pidtof_closest_beta_50ps    [i] = _pfo.pfo_pidtof_closest_beta_50ps    [iLPFO];
    data->lpfo_pidtof_closest_beta_100ps   [i] = _pfo.pfo_pidtof_closest_beta_100ps   [iLPFO];
    data->lpfo_pidtof_fastest_beta_0ps     [i] = _pfo.pfo_pidtof_fastest_beta_0ps     [iLPFO];
    data->lpfo_pidtof_fastest_beta_10ps    [i] = _pfo.pfo_pidtof_fastest_beta_10ps    [iLPFO];
    data->lpfo_pidtof_fastest_beta_50ps    [i] = _pfo.pfo_pidtof_fastest_beta_50ps    [iLPFO];
    data->lpfo_pidtof_fastest_beta_100ps   [i] = _pfo.pfo_pidtof_fastest_beta_100ps   [iLPFO];
    data->lpfo_pidtof_cylfit_beta_0ps      [i] = _pfo.pfo_pidtof_cylfit_beta_0ps      [iLPFO];
    data->lpfo_pidtof_cylfit_beta_10ps     [i] = _pfo.pfo_pidtof_cylfit_beta_10ps     [iLPFO];
    data->lpfo_pidtof_cylfit_beta_50ps     [i] = _pfo.pfo_pidtof_cylfit_beta_50ps     [iLPFO];
    data->lpfo_pidtof_cylfit_beta_100ps    [i] = _pfo.pfo_pidtof_cylfit_beta_100ps    [iLPFO];
    data->lpfo_pidtof_closestfit_beta_0ps  [i] = _pfo.pfo_pidtof_closestfit_beta_0ps  [iLPFO];
    data->lpfo_pidtof_closestfit_beta_10ps [i] = _pfo.pfo_pidtof_closestfit_beta_10ps [iLPFO];
    data->lpfo_pidtof_closestfit_beta_50ps [i] = _pfo.pfo_pidtof_closestfit_beta_50ps [iLPFO];
    data->lpfo_pidtof_closestfit_beta_100ps[i] = _pfo.pfo_pidtof_closestfit_beta_100ps[iLPFO];
    i++;
  }

}

void EventAnalyzer::Analyze(Long64_t entry)
{

  Bool_t debug = (entry == 7515);

  // PFO Analysis
    PFOTools pfot( &_pfo );
    if ( !pfot.ValidPFO() ) return;

  // Fill raw LPFO info
    // if(debug) {
      WriteLPFOVariables(pfot,&_tree_lpfo);      
    // }
    


  // Selections
    vector<Bool_t> CutTrigger;

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
        Bool_t charge_opposite = iSPFO_K.charge * pfot.LPFO[ijet].charge < 0;
        Bool_t momentum_above  = iSPFO_K.p_mag > 10;
        if( charge_opposite && momentum_above ) is_gluon_K[ijet] = true;
      }
    }
    for ( auto ibool : is_gluon_K ){
      if( ibool ) is_there_a_gluon_K = true;
    }
    CutTrigger.push_back(!is_there_a_gluon_K);

  // dEdx dist PDG check
    enum PDGConfig { K_K, K_Pi, Pi_Pi };
    Int_t dEdx_pdg_check = -1;
    
    if(   pfot.isKaon(pfot.LPFO[0]) && pfot.isKaon(pfot.LPFO[1]) )   dEdx_pdg_check = K_K;
    if(   pfot.isPion(pfot.LPFO[0]) && pfot.isPion(pfot.LPFO[1]) )   dEdx_pdg_check = Pi_Pi;
    if( ( pfot.isKaon(pfot.LPFO[0]) && pfot.isPion(pfot.LPFO[1]) ) ||
        ( pfot.isKaon(pfot.LPFO[1]) && pfot.isPion(pfot.LPFO[0]) ) ) dEdx_pdg_check = K_Pi;

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





  // Kaon
    if( dEdx_pdg_check == K_K && all_checks ){
    // if( debug && dEdx_pdg_check == K_K ){
      cout << "K_K event!" << endl;
    }


    _hTree_LPFO->Fill();


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

/*
  // Some debug comments we might use later...
    if( debug ) {

      cout << entry << endl;
      cout << "check: " << LPFO_double_quality << "\n";

      for ( auto iLPFO : pfot.LPFO ){
        cout << "mom: " << pfot.is_momentum(iLPFO,20,60) << ", tpc: " << pfot.is_tpc_hits(iLPFO,210) << ", offset: " << pfot.is_offset_small(iLPFO,1.0) << endl;
      }
      cout << "is there a gluon K? -> " << is_there_a_gluon_K << endl;
      cout << "is charge check? -> " << charge_check << endl;
      cout << "all check? -> " << all_checks << endl;

    }
*/
