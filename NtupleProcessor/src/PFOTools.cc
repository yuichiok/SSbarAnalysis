
/*------------------------------------------------------------------------------
PFOTools.cpp
 Created : 2022-09-11  okugawa
------------------------------------------------------------------------------*/

#include <iostream>
#include <algorithm>
#include <TString.h>
#include <TFile.h> 

#include "PFOTools.hh"
#include "VectorTools.hh"

using std::cout;   using std::endl;

namespace QQbarAnalysis
{
  template<typename T>
  void pop_front(std::vector<T>& vec)
  {
      assert(!vec.empty());
      vec.front() = std::move(vec.back());
      vec.pop_back();
  }

  PFOTools::PFOTools() {}

  PFOTools::PFOTools( MC_QQbar *ptr, TString fnac )
  : mc_data(ptr)
  {
      _anCfg.SetConfig(fnac);
      InitializeMCTools( mc_data );
  }

  PFOTools::PFOTools( MC_QQbar *ptr_mc, PFO_QQbar *ptr, TString fnac )
  : mc_data(ptr_mc), data(ptr)
  {
      _anCfg.SetConfig(fnac);
      InitializePFOTools( mc_data, data );
  }

  void PFOTools::InitializeMCTools( MC_QQbar *mc_data )
  {
    for (int iqq=0; iqq < 2; iqq++){
      VectorTools mcqv(mc_data->mc_quark_px[iqq],mc_data->mc_quark_py[iqq],mc_data->mc_quark_pz[iqq],mc_data->mc_quark_E[iqq]);
      mc_quark[iqq].vt    = mcqv;
      mc_quark[iqq].p_mag = (Float_t) mcqv.v3().R();
      mc_quark[iqq].cos   = std::cos(mcqv.v3().Theta());
      mc_quark[iqq].qcos  = (mc_data->mc_quark_charge[iqq] < 0) ? mc_quark[iqq].cos : -mc_quark[iqq].cos;
    }
    for ( int istable=0; istable < mc_data->mc_stable_n; istable++){
      VectorTools mcsv(mc_data->mc_stable_px[istable],mc_data->mc_stable_py[istable],mc_data->mc_stable_pz[istable],mc_data->mc_stable_E[istable]);
      mc_stable[istable].vt    = mcsv;
      mc_stable[istable].p_mag = (Float_t) mcsv.v3().R();
      mc_stable[istable].cos   = std::cos(mcsv.v3().Theta());
      mc_stable[istable].qcos  = (mc_data->mc_stable_charge[istable] < 0) ? mc_stable[istable].cos : -mc_stable[istable].cos;
    }

  }

  void PFOTools::InitializePFOTools( MC_QQbar *mc_data, PFO_QQbar *data )
  {

    for (int ipfo=0; ipfo < data->pfo_n; ipfo++)
    {
      // Make sure PFO belongs to either jet0 or jet 1
      Bool_t no_match01 = (data->pfo_match[ipfo] != 0 && data->pfo_match[ipfo] != 1);

      // Make suer PFO has only one reconstructed track to avoid (lambda/sigma)
      Bool_t mult_track = (data->pfo_ntracks[ipfo] != 1);

      // if dEdx dist is 0
      Bool_t dEdx_dist_bad = is_dEdxdist_bad(data->pfo_piddedx_e_dedxdist[ipfo],
                          data->pfo_piddedx_mu_dedxdist[ipfo],
                          data->pfo_piddedx_pi_dedxdist[ipfo],
                          data->pfo_piddedx_k_dedxdist[ipfo],
                          data->pfo_piddedx_p_dedxdist[ipfo]);

      if( no_match01 || mult_track || dEdx_dist_bad ) continue;

      VectorTools vt(data->pfo_px[ipfo], data->pfo_py[ipfo], data->pfo_pz[ipfo], data->pfo_E[ipfo]);

      // Initialize & categorize PFO variables with different jets
      const int ijet = data->pfo_match[ipfo];

      PFO.ipfo  = ipfo;
      PFO.vt    = vt;
      PFO.p_mag = (Float_t) vt.v3().R();
      PFO.pv    = sqrt(data->pfo_d0[ipfo] * data->pfo_d0[ipfo] + data->pfo_z0[ipfo] * data->pfo_z0[ipfo]);

      PFO.pfo_match = data->pfo_match[ipfo];
      PFO.pfo_truejet_pdg = data->pfo_truejet_pdg[ipfo];
      PFO.pfo_truejet_type = data->pfo_truejet_type[ipfo];
      PFO.pfo_pdgcheat = data->pfo_pdgcheat[ipfo];
      PFO.pfo_pdgcheat_id = data->pfo_pdgcheat_id[ipfo];
      PFO.pfo_nparents = data->pfo_nparents[ipfo];
      for (int iparent=0; iparent<PFO.pfo_nparents; iparent++){
      PFO.pfo_pdgcheat_parent[iparent] = data->pfo_pdgcheat_parent[ipfo][iparent];
      }
      PFO.pfo_E = data->pfo_E[ipfo];
      PFO.pfo_px = data->pfo_px[ipfo];
      PFO.pfo_py = data->pfo_py[ipfo];
      PFO.pfo_pz = data->pfo_pz[ipfo];
      PFO.pfo_m = data->pfo_m[ipfo];
      PFO.pfo_type = data->pfo_type[ipfo];
      PFO.pfo_isoverlay = data->pfo_isoverlay[ipfo];
      PFO.pfo_isisr = data->pfo_isisr[ipfo];
      PFO.pfo_vtx = data->pfo_vtx[ipfo];
      PFO.pfo_charge = data->pfo_charge[ipfo];
      PFO.pfo_ntracks = data->pfo_ntracks[ipfo];
      PFO.pfo_tpc_hits = data->pfo_tpc_hits[ipfo];
      PFO.pfo_dedx = data->pfo_dedx[ipfo];
      PFO.pfo_dedxerror = data->pfo_dedxerror[ipfo];
      PFO.pfo_d0 = data->pfo_d0[ipfo];
      PFO.pfo_d0error = data->pfo_d0error[ipfo];
      PFO.pfo_z0 = data->pfo_z0[ipfo];
      PFO.pfo_z0error = data->pfo_z0error[ipfo];
      PFO.pfo_phi = data->pfo_phi[ipfo];
      PFO.pfo_phierror = data->pfo_phierror[ipfo];
      PFO.pfo_omega = data->pfo_omega[ipfo];
      PFO.pfo_omegaerror = data->pfo_omegaerror[ipfo];
      PFO.pfo_tanlambda = data->pfo_tanlambda[ipfo];
      PFO.pfo_tanlambdaerror = data->pfo_tanlambdaerror[ipfo];
      PFO.pfo_chi2 = data->pfo_chi2[ipfo];
      PFO.pfo_ndf = data->pfo_ndf[ipfo];
      for (int i=0; i<3; i++){
      PFO.pfo_vtxpt[i] = data->pfo_vtxpt[ipfo][i];
      PFO.pfo_endpt[i] = data->pfo_endpt[ipfo][i];
      }
      PFO.pfo_pid = data->pfo_pid[ipfo];
      PFO.pfo_pid_likelihood = data->pfo_pid_likelihood[ipfo];
      PFO.pfo_pid_eprob = data->pfo_pid_eprob[ipfo];
      PFO.pfo_pid_muprob = data->pfo_pid_muprob[ipfo];
      PFO.pfo_pid_piprob = data->pfo_pid_piprob[ipfo];
      PFO.pfo_pid_kprob = data->pfo_pid_kprob[ipfo];
      PFO.pfo_pid_pprob = data->pfo_pid_pprob[ipfo];
      PFO.pfo_pid_hprob = data->pfo_pid_hprob[ipfo];
      PFO.pfo_piddedx = data->pfo_piddedx[ipfo];
      PFO.pfo_piddedx_likelihood = data->pfo_piddedx_likelihood[ipfo];
      PFO.pfo_piddedx_eprob = data->pfo_piddedx_eprob[ipfo];
      PFO.pfo_piddedx_muprob = data->pfo_piddedx_muprob[ipfo];
      PFO.pfo_piddedx_piprob = data->pfo_piddedx_piprob[ipfo];
      PFO.pfo_piddedx_kprob = data->pfo_piddedx_kprob[ipfo];
      PFO.pfo_piddedx_pprob = data->pfo_piddedx_pprob[ipfo];
      PFO.pfo_piddedx_hprob = data->pfo_piddedx_hprob[ipfo];
      PFO.pfo_piddedx_e_dedxdist = data->pfo_piddedx_e_dedxdist[ipfo];
      PFO.pfo_piddedx_mu_dedxdist = data->pfo_piddedx_mu_dedxdist[ipfo];
      PFO.pfo_piddedx_pi_dedxdist = data->pfo_piddedx_pi_dedxdist[ipfo];
      PFO.pfo_piddedx_k_dedxdist = data->pfo_piddedx_k_dedxdist[ipfo];
      PFO.pfo_piddedx_p_dedxdist = data->pfo_piddedx_p_dedxdist[ipfo];
      PFO.pfo_piddedx_e_lkhood = data->pfo_piddedx_e_lkhood[ipfo];
      PFO.pfo_piddedx_mu_lkhood = data->pfo_piddedx_mu_lkhood[ipfo];
      PFO.pfo_piddedx_pi_lkhood = data->pfo_piddedx_pi_lkhood[ipfo];
      PFO.pfo_piddedx_k_lkhood = data->pfo_piddedx_k_lkhood[ipfo];
      PFO.pfo_piddedx_p_lkhood = data->pfo_piddedx_p_lkhood[ipfo];
      PFO.pfo_pidtof_p_at_calo = data->pfo_pidtof_p_at_calo[ipfo];
      PFO.pfo_pidtof_closest_beta_0ps = data->pfo_pidtof_closest_beta_0ps[ipfo];
      PFO.pfo_pidtof_closest_beta_10ps = data->pfo_pidtof_closest_beta_10ps[ipfo];
      PFO.pfo_pidtof_closest_beta_50ps = data->pfo_pidtof_closest_beta_50ps[ipfo];
      PFO.pfo_pidtof_closest_beta_100ps = data->pfo_pidtof_closest_beta_100ps[ipfo];
      PFO.pfo_pidtof_fastest_beta_0ps = data->pfo_pidtof_fastest_beta_0ps[ipfo];
      PFO.pfo_pidtof_fastest_beta_10ps = data->pfo_pidtof_fastest_beta_10ps[ipfo];
      PFO.pfo_pidtof_fastest_beta_50ps = data->pfo_pidtof_fastest_beta_50ps[ipfo];
      PFO.pfo_pidtof_fastest_beta_100ps = data->pfo_pidtof_fastest_beta_100ps[ipfo];
      PFO.pfo_pidtof_cylfit_beta_0ps = data->pfo_pidtof_cylfit_beta_0ps[ipfo];
      PFO.pfo_pidtof_cylfit_beta_10ps = data->pfo_pidtof_cylfit_beta_10ps[ipfo];
      PFO.pfo_pidtof_cylfit_beta_50ps = data->pfo_pidtof_cylfit_beta_50ps [ipfo]; 
      PFO.pfo_pidtof_cylfit_beta_100ps = data->pfo_pidtof_cylfit_beta_100ps[ipfo];
      PFO.pfo_pidtof_closestfit_beta_0ps = data->pfo_pidtof_closestfit_beta_0ps[ipfo];
      PFO.pfo_pidtof_closestfit_beta_10ps = data->pfo_pidtof_closestfit_beta_10ps[ipfo];
      PFO.pfo_pidtof_closestfit_beta_50ps = data->pfo_pidtof_closestfit_beta_50ps[ipfo];
      PFO.pfo_pidtof_closestfit_beta_100ps = data->pfo_pidtof_closestfit_beta_100ps[ipfo];

      PFO.dEdx_dist_pdg = Get_dEdx_dist_PID( PFO.pfo_piddedx_k_dedxdist, PFO.pfo_piddedx_pi_dedxdist, PFO.pfo_piddedx_p_dedxdist );
      PFO.cos           = std::cos(PFO.vt.v3().Theta());
      PFO.qcos          = (PFO.pfo_charge < 0) ? PFO.cos : -PFO.cos;

      Valid_PFOs.push_back(PFO);
      PFO_jet[ijet].push_back(PFO);

      // Cheated PFO
      for( auto i_lmode : PFO_mode ){
        if( abs(PFO.pfo_pdgcheat) == PFO_type_map_rev.at(i_lmode) ) PFO_unsorted_subjet_cheat[i_lmode][ijet].push_back(PFO);
      }

    }

    for (int ijet=0; ijet < 2; ijet++){
      
      // sort jet
      PFO_sorted_jet[ijet] = SortJet(PFO_jet[ijet]);
      if( PFO_sorted_jet[ijet].size() ) LPFO.push_back(PFO_sorted_jet[ijet].at(0));
      
      // get subjet
      for ( auto i_lmode : PFO_mode ){
        
        PFO_subjet[i_lmode][ijet] = GetSubjet(ijet, i_lmode);

        if( PFO_unsorted_subjet_cheat[i_lmode][ijet].size() ){
          PFO_subjet_cheat[i_lmode][ijet] = SortJet(PFO_unsorted_subjet_cheat.at(i_lmode).at(ijet));
        }
      }

      // cheat

    }

    fCorrection->SetParameters(pars);

  }

  vector<PFO_Info> PFOTools::SortJet( vector<PFO_Info> jet )
  {
      std::sort(jet.begin(), jet.end(),std::greater<PFO_Info>());
      return jet;
  }

  vector<PFO_Info> PFOTools::GetJet( int ijet )
  {
      return PFO_jet[ijet];
  }

  vector<PFO_Info> PFOTools::GetSortedJet( int ijet )
  {
      vector<PFO_Info> sorted_jet = PFO_jet[ijet];
      sorted_jet = SortJet( sorted_jet );
      return sorted_jet;
  }

  vector<PFO_Info> PFOTools::GetSubjet( int ijet, TString lmode )
  {
    vector<PFO_Info> subjet;
    for ( const auto iPFO : PFO_sorted_jet[ijet] ){
      if( is_PID( lmode, iPFO ) ) subjet.push_back( iPFO );
    }
    return subjet;
  }

  Int_t PFOTools::Get_dEdx_dist_PID( Float_t kdEdx_dist, Float_t pidEdx_dist, Float_t pdEdx_dist )
  {
      const static int nparticles = 3;
      Int_t        k_pi_p[nparticles] = { 321, 211, 2212 };
      Float_t k_pi_p_dEdx[nparticles] = { abs(kdEdx_dist), abs(pidEdx_dist), abs(pdEdx_dist) };

      Float_t  min_dEdx_dist = *std::min_element(k_pi_p_dEdx, k_pi_p_dEdx + nparticles );
      Int_t   imin_dEdx_dist =  std::find( k_pi_p_dEdx, k_pi_p_dEdx + nparticles, min_dEdx_dist ) - k_pi_p_dEdx;

      return k_pi_p[imin_dEdx_dist];

  }

  Float_t PFOTools::Get_dEdx_dist( PFO_Info iPFO, TString particle )
  {
         if( particle == "K" )  { return iPFO.pfo_piddedx_k_dedxdist;  }
    else if( particle == "Pi" ) { return iPFO.pfo_piddedx_pi_dedxdist; }
    else if( particle == "p" )  { return iPFO.pfo_piddedx_p_dedxdist;  }
    else if( particle == "e" )  { return iPFO.pfo_piddedx_e_dedxdist;  }
    else if( particle == "mu" )  { return iPFO.pfo_piddedx_mu_dedxdist;  }
    else {return -1000;}
  }

  Bool_t PFOTools::isKaon( PFO_Info iPFO )
  {
      return iPFO.dEdx_dist_pdg == 321;
  }

  Bool_t PFOTools::isPion( PFO_Info iPFO )
  {
    // return iPFO.dEdx_dist_pdg == 211;
    // return (iPFO.dEdx_dist_pdg == 211) && (0 < iPFO.pfo_piddedx_pi_dedxdist);
    // return (iPFO.dEdx_dist_pdg == 211) && (-2 < iPFO.pfo_piddedx_pi_dedxdist) && (iPFO.pfo_piddedx_pi_dedxdist < 2);
    Float_t dedx_pi_mean  = 1.82249E-1;
    Float_t dedx_pi_sigma = 9.60323E-3;
    Float_t dedx_pi_min   = dedx_pi_mean - dedx_pi_sigma;
    Float_t dedx_pi_max   = dedx_pi_mean + dedx_pi_sigma;
    return (dedx_pi_min < iPFO.pfo_piddedx_pi_dedxdist) && (iPFO.pfo_piddedx_pi_dedxdist < dedx_pi_max);

    // if(abs(iPFO.cos)>0.9) return (iPFO.dEdx_dist_pdg == 211);
    // else{
    //   Double_t correction = fCorrection->Eval(iPFO.cos);
    //   return (iPFO.pfo_piddedx_pi_dedxdist > correction);
    // }

  }

  Bool_t PFOTools::isProton( PFO_Info iPFO )
  {
    return iPFO.dEdx_dist_pdg == 2212;
  }

  Bool_t PFOTools::is_cheatNoOthers( PFO_Info iPFO )
  {
    return (abs(iPFO.pfo_pdgcheat) != 11) && (abs(iPFO.pfo_pdgcheat) != 13);
  }

  Bool_t PFOTools::is_PID( TString lmode, PFO_Info iPFO )
  {
    if     ( lmode == "K"  ){ return isKaon(iPFO);   }
    else if( lmode == "Pi" ){ return isPion(iPFO);   }
    else if( lmode == "p"  ){ return isProton(iPFO); }
    else                    { return false; }
  }

  Bool_t PFOTools::is_jet_mult_non0()
  {
    for ( auto iPFO_jet : PFO_jet ) {
        if( !iPFO_jet.size() ) return false;
    }
    return true;
  }

  Bool_t PFOTools::is_charge_config( ChargeConfig cc, Int_t charge0 , Int_t charge1 )
  {

    Int_t charge_config = charge0 * charge1;

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

  Bool_t PFOTools::is_momentum( PFO_Info iPFO, Float_t MINP_CUT, Float_t MAXP_CUT )
  {
    return ( MINP_CUT < iPFO.p_mag && iPFO.p_mag < MAXP_CUT);
  }

  Bool_t PFOTools::is_tpc_hits( PFO_Info iPFO, Int_t MIN_TPC_HITS )
  {
    return ( MIN_TPC_HITS < iPFO.pfo_tpc_hits );
  }

  Bool_t PFOTools::is_offset_small( PFO_Info iPFO, Int_t MAX_OFFSET )
  {
    return ( iPFO.pv < MAX_OFFSET );
  }

  Bool_t PFOTools::is_dEdxdist_bad( Float_t e_dist, Float_t mu_dist, Float_t pi_dist, Float_t k_dist, Float_t p_dist )
  {
    if( !e_dist )   return 1;
    if( !mu_dist )  return 1;
    if( !pi_dist )  return 1;
    if( !k_dist )   return 1;
    if( !p_dist )   return 1;
    return 0;
  }

}