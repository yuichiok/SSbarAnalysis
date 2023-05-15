#ifndef GUARD_HistManager_h
#define GUARD_HistManager_h

#include <iostream>
#include <map>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TFile.h> 

const int NSLABS = 15;

class HistManager
{
  public:
    HistManager();
    ~HistManager(){}

  // Methods
    virtual void InitializeHists();
    virtual void Hist2List();
    virtual void WriteLists( TFile * output );

  // Declear histograms
  // h1 hist
    enum h1_index {
      gen_q_cos,
      gen_q_qcos,

      gen_K_cos,
      gen_K_qcos,

      cheat_K_cos,
      cheat_K_qcos,

      cheat_Pi_cos,
      cheat_Pi_qcos,

      reco_K_cos,
      reco_K_qcos,
      reco_K_scos,
      reco_K_mom,
      reco_K_pdgcheat,
      gen_reco_K_sep_cos,
      jet_reco_K_sep_cos,

      reco_Pi_cos,
      reco_Pi_qcos,
      reco_Pi_scos,
      reco_Pi_mom,
      reco_Pi_pdgcheat,
      gen_reco_Pi_sep_cos,
      jet_reco_Pi_sep_cos,

      good_reco_Pi_endpt,
      good_reco_Pi_tpchits,
      good_reco_Pi_pidedx_dist,
      good_reco_Pi_kdedx_dist,
      bad_reco_Pi_endpt,
      bad_reco_Pi_tpchits,
      bad_reco_Pi_pidedx_dist,
      bad_reco_Pi_kdedx_dist,

      reco_sum_jetE,
      reco_jet_sep,

      lpfo_gen_K_mom,
      lpfo_reco_K_mom,

      gen_N_K_cos,
      reco_N_K_cos,
      N_K_corr_cos,

      gen_N_Pi_cos,
      reco_N_Pi_cos,
      N_Pi_corr_cos,

      dummy_h1,
      Last_h1 = dummy_h1
    };
    TH1F * h1[Last_h1];

    enum h1_pq_index {
      // KID
      acc_KK,
      rej_KK,

      // PiID
      acc_PiPi,
      rej_PiPi,

      dummy_pq,
      Last_h1_pq = dummy_pq
    };
    TH1F * h1_pq[Last_h1_pq];

    enum h1_particle_ratio_index {
      K_rate_gen,
      pi_rate_gen,
      p_rate_gen,

      K_rate_reco,
      pi_rate_reco,
      p_rate_reco,

      dummy_particle_ratio,
      Last_h1_dummy_particle_ratio = dummy_particle_ratio
    };
    TH1F * h1_particle_ratio[Last_h1_dummy_particle_ratio];

  // h2 hist
    enum h2_index {

      gen_K_p_cos,
      reco_K_p_cos,

      nK_gen_reco,
      npi_gen_reco,
      np_gen_reco,

      stable_K_cos,
      purity_K_cos,

      stable_Pi_cos,
      purity_Pi_cos,

      dummy_h2,
      Last_h2 = dummy_h2
    };
    TH2F * h2[Last_h2];

    enum h2_jet_index {

      jet_mult_cos,
      jet_mult_cos_noISR,

      dummy_h2_jet,
      Last_dummy_h2_jet = dummy_h2_jet
    };
    TH2F * h2_jet[Last_dummy_h2_jet];

    enum h2_particle_ratio_cos_index {
      K_rate_cos_gen,
      pi_rate_cos_gen,
      p_rate_cos_gen,

      K_rate_cos_reco,
      pi_rate_cos_reco,
      p_rate_cos_reco,

      dummy_particle_ratio_cos,
      Last_h2_dummy_particle_ratio_cos = dummy_particle_ratio_cos
    };
    TH2F * h2_particle_ratio_cos[Last_h2_dummy_particle_ratio_cos];

    enum h2_dEdx_index {

      gen_K_dEdx_p,
      gen_pi_dEdx_p,
      gen_p_dEdx_p,

      dummy_dEdx,
      Last_h2_dEdx = dummy_dEdx
    };
    TH2F * h2_dEdx[Last_h2_dEdx];

    enum particle_dEdx_index {

      // KID
      gen_ipart_reco_K_dEdx_p,
      gen_ipart_KdEdx_dist_cos,
      gen_ipart_reco_K_KdEdx_dist_cos,

      // PiID
      gen_ipart_reco_Pi_dEdx_p,
      gen_ipart_PidEdx_dist_cos,
      gen_ipart_reco_Pi_PidEdx_dist_cos,

      dummy_particle_dEdx,
      Last_h2_particle_dEdx = dummy_particle_dEdx
    };
    TString particle_list[5] = {"K", "pi", "p", "e", "mu"};
    enum particle_List_index {
      kKaon,
      kPion,
      kProton,
      kElectron,
      kMuon,
      dummy_particle_List,
      Last_particle_List = dummy_particle_List
    };
    TH2F * h2_particle_dEdx[Last_h2_particle_dEdx][Last_particle_List];

  private:
    TList* hList1                = new TList();
    TList* hList1_pq             = new TList();
    TList* hList1_particle_ratio = new TList();
    TList* hList2                = new TList();
    TList* hList2_jet            = new TList();
    TList* hList2_dEdx           = new TList();
    TList* hList2_particle_dEdx  = new TList();

};

#endif