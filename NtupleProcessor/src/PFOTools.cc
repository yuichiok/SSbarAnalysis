
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
      if (data->pfo_match[ipfo] < 0 || 1 < data->pfo_match[ipfo]) continue;

      // Make suer PFO has only one reconstructed track to avoid (lambda/sigma)
      if (data->pfo_ntracks[ipfo] != 1) continue;

      // if dEdx dist is 0
      if( is_dEdxdist_bad(data->pfo_piddedx_e_dedxdist[ipfo],
                          data->pfo_piddedx_mu_dedxdist[ipfo],
                          data->pfo_piddedx_pi_dedxdist[ipfo],
                          data->pfo_piddedx_k_dedxdist[ipfo],
                          data->pfo_piddedx_p_dedxdist[ipfo]) ) continue;

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

      PFO_jet[ijet].push_back(PFO);
      if( abs(PFO.pfo_pdgcheat) == 321 ) PFO_cheat_Ks[ijet].push_back(PFO);
      if( abs(PFO.pfo_pdgcheat) == 211 ) PFO_cheat_Pis[ijet].push_back(PFO);
      
    }

    if( ValidPFO() ){
      for (int ijet=0; ijet < 2; ijet++){
        LPFO[ijet]        = GetSortedJet(ijet).at(0);
        LPFO_["K"][ijet]  = Get_Particle_LPFO(ijet,kKaon);
        LPFO_["Pi"][ijet] = Get_Particle_LPFO(ijet,kPion);

        KLPFO[ijet]       = Get_Particle_LPFO(ijet,kKaon);
        PiLPFO[ijet]      = Get_Particle_LPFO(ijet,kPion);

        // Reconstructed PFO
        if( PFO_jet[ijet].size() > 1 ){
          SPFOs[ijet] = GetSortedJet(ijet);
          // SPFOs[ijet].erase(SPFOs[ijet].begin());
          pop_front(SPFOs[ijet]); // faster algorithm wise?

          std::copy_if(SPFOs[ijet].begin(), SPFOs[ijet].end(), std::back_inserter(SPFOs_["K"][ijet]), [](PFO_Info iPFO) {
              return isKaon(iPFO);
          });

          std::copy_if(SPFOs[ijet].begin(), SPFOs[ijet].end(), std::back_inserter(SPFOs_["Pi"][ijet]), [](PFO_Info iPFO) {
              return isPion(iPFO);
          });

          // delete here afterwards

          std::copy_if(SPFOs[ijet].begin(), SPFOs[ijet].end(), std::back_inserter(SPFOs_K[ijet]), [](PFO_Info iPFO) {
              return isKaon(iPFO);
          });

          std::copy_if(SPFOs[ijet].begin(), SPFOs[ijet].end(), std::back_inserter(SPFOs_Pi[ijet]), [](PFO_Info iPFO) {
              return isPion(iPFO);
          });

          //////

        }

        // Cheated PFO
        if( PFO_cheat_Ks[0].size() && PFO_cheat_Ks[1].size() ){
          cheat_KLPFO[ijet] = SortJet(PFO_cheat_Ks[ijet]).at(0);
          if( PFO_cheat_Ks[ijet].size() > 1 ){
            SPFOs_cheat_K[ijet] = PFO_cheat_Ks[ijet];
            pop_front(SPFOs_cheat_K[ijet]);
          }
        }

        if( PFO_cheat_Pis[0].size() && PFO_cheat_Pis[1].size() ){
          cheat_PiLPFO[ijet] = SortJet(PFO_cheat_Pis[ijet]).at(0);
          if( PFO_cheat_Pis[ijet].size() > 1 ){
            SPFOs_cheat_Pi[ijet] = PFO_cheat_Pis[ijet];
            pop_front(SPFOs_cheat_Pi[ijet]);
          }
        }

      }
    }

  }

  vector<PFO_Info> PFOTools::SortJet( vector<PFO_Info> jet )
  {
      std::sort(jet.begin(), jet.end(),std::greater<PFO_Info>());
      return jet;
  }

  Bool_t PFOTools::ValidPFO()
  {
      for ( auto iPFO_jet : PFO_jet ) {
          if( !iPFO_jet.size() ) return false;
      }
      return true;
  }

  vector<PFO_Info> PFOTools::Get_Valid_PFOs()
  {
    // Combines PFO lists in 2 jets
      std::vector<PFO_Info> PFO_Collection;
      PFO_Collection.reserve( PFO_jet[0].size() + PFO_jet[1].size() );
      PFO_Collection.insert( PFO_Collection.begin(), PFO_jet[0].begin(), PFO_jet[0].end() );
      PFO_Collection.insert( PFO_Collection.end(), PFO_jet[1].begin(), PFO_jet[1].end() );
      return PFO_Collection;
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

  PFO_Info PFOTools::Get_Particle_LPFO( int ijet, ParticleID pdg )
  {
      vector<PFO_Info> sorted_jet = GetSortedJet(ijet);
      for ( auto iPFO : sorted_jet ) {

        switch ( pdg )
        {
          case kKaon:
            if ( isKaon(iPFO) ) return iPFO;
            break;
          case kPion:
            if ( isPion(iPFO) ) return iPFO;
            break;
          case kProton:
            if ( isProton(iPFO) ) return iPFO;
            break;
          
          default:
            break;
        }

      }

      return sorted_jet.at(0);
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

  Bool_t PFOTools::isKaon( PFO_Info iPFO )
  {
      return iPFO.dEdx_dist_pdg == 321;
  }

  Bool_t PFOTools::isPion( PFO_Info iPFO )
  {
    // return iPFO.dEdx_dist_pdg == 211;
    return (iPFO.dEdx_dist_pdg == 211) && (0 < iPFO.pfo_piddedx_pi_dedxdist);
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

  Bool_t PFOTools::is_PID_config( TString lmode )
  {
    return is_PID( lmode, LPFO_[lmode][0] ) && is_PID( lmode, LPFO_[lmode][1] );
  }

  Bool_t PFOTools::is_charge_config( ChargeConfig cc, Int_t charge0 , Int_t charge1 )
  {

    // Int_t charge_config = LPFO[0].pfo_charge * LPFO[1].pfo_charge;
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

  Bool_t PFOTools::LPFO_Quality_checks( PFO_Info iPFO )
  {
    vector<Bool_t> CutTrigger;

    CutTrigger.push_back( is_momentum( iPFO, _anCfg.LPFO_p_min, _anCfg.LPFO_p_max ) );  // MIN/MAX momentum check
    CutTrigger.push_back( is_tpc_hits( iPFO, _anCfg.PFO_TPCHits_max ) );                // Number of TPC hit check
    CutTrigger.push_back( is_offset_small( iPFO, _anCfg.PFO_offset_max ) );             // Offset distance check
    
    for (auto trigger : CutTrigger ){
      if (!trigger) { return false; }
    }

    return true;

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
    if( !e_dist ) return 1;
    if( !mu_dist ) return 1;
    if( !pi_dist ) return 1;
    if( !k_dist ) return 1;
    if( !p_dist ) return 1;
    return 0;
  }

  Bool_t PFOTools::is_ss()
  {
    if ( KLPFO[0].p_mag > PiLPFO[0].p_mag &&
        KLPFO[1].p_mag > PiLPFO[1].p_mag ){
          return true;
    }else{
      return false;
    }
    
  }

  Bool_t PFOTools::is_uu_dd()
  {
    if ( PiLPFO[0].p_mag > KLPFO[0].p_mag &&
        PiLPFO[1].p_mag > KLPFO[1].p_mag ){
          return true;
    }else{
      return false;
    }
  }
}