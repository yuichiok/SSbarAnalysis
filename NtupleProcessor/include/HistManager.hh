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

      reco_K_cos,
      reco_K_qcos,

      reco_K_scos,

      gen_reco_K_sep_cos,

      reco_sum_jetE,
      reco_jet_sep,

      lpfo_gen_K_mom,
      lpfo_reco_K_mom,

      gen_N_K_cos,
      reco_N_K_cos,
      N_K_corr_cos,

      gen_N_K_cos2,
      reco_N_K_cos2,
      N_K_corr_cos2,

      dummy_h1,
      Last_h1 = dummy_h1
    };
    TH1F * h1[Last_h1];

    enum h1_pq_index {
      acc_KK,
      rej_KK,

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

      stable_cos,
      purity_cos,

      dummy_h2,
      Last_h2 = dummy_h2
    };
    TH2F * h2[Last_h2];

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

      gen_K_reco_K_dEdx_p,
      gen_pi_reco_K_dEdx_p,
      gen_p_reco_K_dEdx_p,

      gen_K_KdEdx_dist_cos,
      gen_pi_KdEdx_dist_cos,
      gen_p_KdEdx_dist_cos,

      reco_K_KdEdx_dist_cos,
      gen_K_reco_K_KdEdx_dist_cos,
      gen_pi_reco_K_KdEdx_dist_cos,
      gen_p_reco_K_KdEdx_dist_cos,

      dummy_dEdx,
      Last_h2_dEdx = dummy_dEdx
    };
    TH2F * h2_dEdx[Last_h2_dEdx];

    enum h1_K_cheat{
      d0_K_cheat_primary,
      d0_sigma_K_cheat_primary,
      d0_sigma_d0_K_cheat_primary,

      z0_K_cheat_primary,
      z0_sigma_K_cheat_primary,
      z0_sigma_z0_K_cheat_primary,

      d0_K_cheat_secondary,
      d0_sigma_K_cheat_secondary,
      d0_sigma_d0_K_cheat_secondary,

      z0_K_cheat_secondary,
      z0_sigma_K_cheat_secondary,
      z0_sigma_z0_K_cheat_secondary,

      pmag_K_cheat,
      cos_theta_K_cheat,

      dummy_K_cheat,
      Last_h1_K_cheat = dummy_K_cheat
    };
    TH1F * h1_K_cheat[Last_h1_K_cheat];

    enum h1_K_reco{
      d0_K_reco_primary,
      d0_sigma_K_reco_primary,
      d0_sigma_d0_K_reco_primary,

      z0_K_reco_primary,
      z0_sigma_K_reco_primary,
      z0_sigma_z0_K_reco_primary,

      d0_K_reco_secondary,
      d0_sigma_K_reco_secondary,
      d0_sigma_d0_K_reco_secondary,

      z0_K_reco_secondary,
      z0_sigma_K_reco_secondary,
      z0_sigma_z0_K_reco_secondary,

      pmag_K_reco,
      cos_theta_K_reco,

      dummy_K_reco,
      Last_h1_K_reco = dummy_K_reco
    };
    TH1F * h1_K_reco[Last_h1_K_reco];

    enum h_PS{
    d0_P_single,
    d0_S_single,
    d0_P_mult,
    d0_S_mult,
    z0_P_single,
    z0_S_single,
    z0_P_mult,
    z0_S_mult,
    dummy_h_PS,
    Last_h_PS = dummy_h_PS
    };
    TH1* h_PS[Last_h_PS];

    enum h_tagging {
      p_ctag,
      s_ctag,
      p_btag,
      s_btag,
      t_ctag,
      t_btag,
      tt_ctag,
      tt_btag,
      jets_info,

      dummy_tagging,
      Last_h_tagging = dummy_tagging
    };
    TH1F * h_tagging[Last_h_tagging];


    enum h_general {
      // nvtx_ctag,
      n_K_ecal,

      dummy_general,
      Last_h_general = dummy_general
    };
    TH1F * h_general[Last_h_general];

  private:
    TList* hList1                = new TList();
    TList* hList1_pq             = new TList();
    TList* hList1_particle_ratio = new TList();
    TList* hList2                = new TList();
    TList* hList2_dEdx           = new TList();

    TList* hList_K_cheat         = new TList();
    TList* hList_K_reco          = new TList();
    TList* hList_general         = new TList();
    TList* hList_tagging         = new TList();
    TList* hList_PS              = new TList();
};

#endif
