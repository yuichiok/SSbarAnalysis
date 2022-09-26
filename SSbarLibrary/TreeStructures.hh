#ifndef GUARD_TreeStructures_h
#define GUARD_TreeStructures_h

/*------------------------------------------------------------------------------
   TreeReader
 Created : 2022-09-09  okugawa
 Main class of TreeReader program.
------------------------------------------------------------------------------*/
#include "VectorTools.hh"

const static int MAX_NPARTICLES = 1000;
const static int NTRUE_JETS     = 5;
const static int NQUARKS        = 2;
const static int NJETS          = 2;

struct MC_QQbar  {

  public:
  // MC VARIABLES
    Float_t mc_quark_E[NQUARKS]      = {0};   Int_t   mc_quark_ps_n   = 0;                        Int_t   mc_stable_n = 0;
    Float_t mc_quark_px[NQUARKS]     = {0};   Float_t mc_quark_ps_E[MAX_NPARTICLES]      = {0};   Float_t mc_stable_E[MAX_NPARTICLES]         = {0};
    Float_t mc_quark_py[NQUARKS]     = {0};   Float_t mc_quark_ps_px[MAX_NPARTICLES]     = {0};   Float_t mc_stable_px[MAX_NPARTICLES]        = {0};
    Float_t mc_quark_pz[NQUARKS]     = {0};   Float_t mc_quark_ps_py[MAX_NPARTICLES]     = {0};   Float_t mc_stable_py[MAX_NPARTICLES]        = {0};
    Float_t mc_quark_m[NQUARKS]      = {0};   Float_t mc_quark_ps_pz[MAX_NPARTICLES]     = {0};   Float_t mc_stable_pz[MAX_NPARTICLES]        = {0};
    Float_t mc_quark_pdg[NQUARKS]    = {0};   Float_t mc_quark_ps_m[MAX_NPARTICLES]      = {0};   Float_t mc_stable_m[MAX_NPARTICLES]         = {0};
    Float_t mc_quark_charge[NQUARKS] = {0};   Float_t mc_quark_ps_pdg[MAX_NPARTICLES]    = {0};   Int_t   mc_stable_pdg[MAX_NPARTICLES]       = {0};
    Float_t mc_ISR_E[NQUARKS]        = {0};   Float_t mc_quark_ps_charge[MAX_NPARTICLES] = {0};   Float_t mc_stable_charge[MAX_NPARTICLES]    = {0}; 
    Float_t mc_ISR_px[NQUARKS]       = {0};   Float_t mc_quark_ps_y12 = 0;                        Int_t   mc_stable_isoverlay[MAX_NPARTICLES] = {0};
    Float_t mc_ISR_py[NQUARKS]       = {0};   Float_t mc_quark_ps_y23 = 0;                        Int_t   mc_stable_isisr[MAX_NPARTICLES]     = {0};
    Float_t mc_ISR_pz[NQUARKS]       = {0};   Float_t mc_quark_ps_d12 = 0;                        Float_t mc_stable_y12 = 0;
    Float_t mc_ISR_m[NQUARKS]        = {0};   Float_t mc_quark_ps_d23 = 0;                        Float_t mc_stable_y23 = 0;
    Float_t mc_ISR_pdg[NQUARKS]      = {0};   Float_t mc_quark_ps_jet_E[NJETS]  = {0};            Float_t mc_stable_d12 = 0;
    Float_t mc_ISR_charge[NQUARKS]   = {0};   Float_t mc_quark_ps_jet_px[NJETS] = {0};            Float_t mc_stable_d23 = 0;
                                              Float_t mc_quark_ps_jet_py[NJETS] = {0};            Float_t mc_stable_jet_E[NJETS]  = {0};
                                              Float_t mc_quark_ps_jet_pz[NJETS] = {0};            Float_t mc_stable_jet_px[NJETS] = {0};
                                                                                                  Float_t mc_stable_jet_py[NJETS] = {0};
                                                                                                  Float_t mc_stable_jet_pz[NJETS] = {0};
                                                                                      
};

struct Jet_QQbar  {

  public:
  // JET VARIABLES
    Float_t truejet_E[NTRUE_JETS]    = {0};   Float_t jet_E[NJETS]    = {0};    Float_t y23 = 0;          Float_t major_thrust_value     = 0;   Float_t major_thrust_axis[3]     = {0};
    Float_t truejet_px[NTRUE_JETS]   = {0};   Float_t jet_px[NJETS]   = {0};    Float_t y12 = 0;          Float_t minor_thrust_value     = 0;   Float_t minor_thrust_axis[3]     = {0};
    Float_t truejet_py[NTRUE_JETS]   = {0};   Float_t jet_py[NJETS]   = {0};    Float_t d23 = 0;          Float_t principle_thrust_value = 0;   Float_t principle_thrust_axis[3] = {0};
    Float_t truejet_pz[NTRUE_JETS]   = {0};   Float_t jet_pz[NJETS]   = {0};    Float_t d12 = 0;          Float_t sphericity = 0;               Float_t sphericity_tensor[3]     = {0};
    Int_t   truejet_type[NTRUE_JETS] = {0};   Float_t jet_btag[NJETS] = {0};    Float_t oblateness = 0; 
    Int_t   truejet_pdg[NTRUE_JETS]  = {0};   Float_t jet_ctag[NJETS] = {0};    Float_t aplanarity = 0;
    
};

struct PFO_QQbar  {

  public:
  // PFO VARIABLES
    Int_t   pfo_n       = 0;
    Int_t   jet_nvtx    = 0;
    Int_t   pfo_n_j1    = 0;
    Int_t   jet_nvtx_j1 = 0;
    Int_t   pfo_n_j2    = 0;
    Int_t   jet_nvtx_j2 = 0;
    Int_t   pfo_match[MAX_NPARTICLES] = {0};
    Int_t   pfo_truejet_pdg[MAX_NPARTICLES] = {0};
    Int_t   pfo_truejet_type[MAX_NPARTICLES] = {0};
    Int_t   pfo_pdgcheat[MAX_NPARTICLES] = {0};
    Int_t   pfo_nparents[MAX_NPARTICLES] = {0};
    Int_t   pfo_pdgcheat_parent[MAX_NPARTICLES][MAX_NPARTICLES] = {{0}};
    Float_t pfo_E[MAX_NPARTICLES] = {0};
    Float_t pfo_E_calo[MAX_NPARTICLES] = {0};
    Float_t pfo_px[MAX_NPARTICLES] = {0};
    Float_t pfo_py[MAX_NPARTICLES] = {0};
    Float_t pfo_pz[MAX_NPARTICLES] = {0};
    Float_t pfo_m[MAX_NPARTICLES] = {0};
    Int_t   pfo_type[MAX_NPARTICLES] = {0};
    Int_t   pfo_isoverlay[MAX_NPARTICLES] = {0};
    Int_t   pfo_isisr[MAX_NPARTICLES] = {0};
    Int_t   pfo_vtx[MAX_NPARTICLES] = {0};
    Int_t   pfo_charge[MAX_NPARTICLES] = {0};
    Int_t   pfo_ntracks[MAX_NPARTICLES] = {0};
    Int_t   pfo_tpc_hits[MAX_NPARTICLES] = {0};
    Float_t pfo_dedx[MAX_NPARTICLES] = {0};
    Float_t pfo_dedxerror[MAX_NPARTICLES] = {0};
    Float_t pfo_d0[MAX_NPARTICLES] = {0};
    Float_t pfo_d0error[MAX_NPARTICLES] = {0};
    Float_t pfo_z0[MAX_NPARTICLES] = {0};
    Float_t pfo_z0error[MAX_NPARTICLES] = {0};
    Float_t pfo_phi[MAX_NPARTICLES] = {0};
    Float_t pfo_phierror[MAX_NPARTICLES] = {0};
    Float_t pfo_omega[MAX_NPARTICLES] = {0};
    Float_t pfo_omegaerror[MAX_NPARTICLES] = {0};
    Float_t pfo_tanlambda[MAX_NPARTICLES] = {0};
    Float_t pfo_tanlambdaerror[MAX_NPARTICLES] = {0};
    Float_t pfo_chi2[MAX_NPARTICLES] = {0};
    Float_t pfo_ndf[MAX_NPARTICLES] = {0};
    Float_t pfo_vtxpt[MAX_NPARTICLES][3] = {{0}};
    Float_t pfo_endpt[MAX_NPARTICLES][3] = {{0}};
    Int_t   pfo_pid[MAX_NPARTICLES] = {0};
    Float_t pfo_pid_likelihood[MAX_NPARTICLES] = {0};
    Float_t pfo_pid_eprob[MAX_NPARTICLES] = {0};
    Float_t pfo_pid_muprob[MAX_NPARTICLES] = {0};
    Float_t pfo_pid_piprob[MAX_NPARTICLES] = {0};
    Float_t pfo_pid_kprob[MAX_NPARTICLES] = {0};
    Float_t pfo_pid_pprob[MAX_NPARTICLES] = {0};
    Float_t pfo_pid_hprob[MAX_NPARTICLES] = {0};
    Int_t   pfo_piddedx[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_likelihood[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_eprob[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_muprob[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_piprob[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_kprob[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_pprob[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_hprob[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_e_dedxdist[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_mu_dedxdist[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_pi_dedxdist[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_k_dedxdist[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_p_dedxdist[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_e_lkhood[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_mu_lkhood[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_pi_lkhood[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_k_lkhood[MAX_NPARTICLES] = {0};
    Float_t pfo_piddedx_p_lkhood[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_p_at_calo[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_closest_beta_0ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_closest_beta_10ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_closest_beta_50ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_closest_beta_100ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_fastest_beta_0ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_fastest_beta_10ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_fastest_beta_50ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_fastest_beta_100ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_cylfit_beta_0ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_cylfit_beta_10ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_cylfit_beta_50ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_cylfit_beta_100ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_closestfit_beta_0ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_closestfit_beta_10ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_closestfit_beta_50ps[MAX_NPARTICLES] = {0};
    Float_t pfo_pidtof_closestfit_beta_100ps[MAX_NPARTICLES] = {0};

};

struct PFO_Info {

  public:
    Int_t   ipfo          = -1;
    VectorTools vt;
    Float_t p_mag         =  0;
    Float_t E             =  0;
    Int_t   charge        =  0;
    Int_t   pdgcheat      =  0;
    Int_t   tpc_hits      =  0;
    Float_t pv            = -2;
    Float_t dEdx          =  0;
    Float_t kdEdx_dist    =  0;
    Float_t pidEdx_dist   =  0;
    Float_t pdEdx_dist    =  0;
    Int_t   dEdx_dist_pdg =  0;
    Float_t cos           = -2;
    Float_t qcos          = -2;
    Float_t q_sep         = -2;
    Float_t qbar_sep      = -2;

    bool operator > (const PFO_Info& apfo) const
    {
        return (p_mag > apfo.p_mag);
    }

    bool operator < (const PFO_Info& apfo) const
    {
        return (p_mag < apfo.p_mag);
    }

};

struct Branch_QQbar  {

  public:
    // List of branches
    TBranch        *b_mc_quark_E;   //!
    TBranch        *b_mc_quark_px;   //!
    TBranch        *b_mc_quark_py;   //!
    TBranch        *b_mc_quark_pz;   //!
    TBranch        *b_mc_quark_m;   //!
    TBranch        *b_mc_quark_pdg;   //!
    TBranch        *b_mc_quark_charge;   //!
    TBranch        *b_mc_ISR_E;   //!
    TBranch        *b_mc_ISR_px;   //!
    TBranch        *b_mc_ISR_py;   //!
    TBranch        *b_mc_ISR_pz;   //!
    TBranch        *b_mc_ISR_m;   //!
    TBranch        *b_mc_ISR_pdg;   //!
    TBranch        *b_mc_ISR_charge;   //!
    TBranch        *b_mc_quark_ps_n;   //!
    TBranch        *b_mc_quark_ps_E;   //!
    TBranch        *b_mc_quark_ps_px;   //!
    TBranch        *b_mc_quark_ps_py;   //!
    TBranch        *b_mc_quark_ps_pz;   //!
    TBranch        *b_mc_quark_ps_m;   //!
    TBranch        *b_mc_quark_ps_pdg;   //!
    TBranch        *b_mc_quark_ps_charge;   //!
    TBranch        *b_mc_quark_ps_y12;   //!
    TBranch        *b_mc_quark_ps_y23;   //!
    TBranch        *b_mc_quark_ps_d12;   //!
    TBranch        *b_mc_quark_ps_d23;   //!
    TBranch        *b_mc_quark_ps_jet_E;   //!
    TBranch        *b_mc_quark_ps_jet_px;   //!
    TBranch        *b_mc_quark_ps_jet_py;   //!
    TBranch        *b_mc_quark_ps_jet_pz;   //!
    TBranch        *b_mc_stable_n;   //!
    TBranch        *b_mc_stable_E;   //!
    TBranch        *b_mc_stable_px;   //!
    TBranch        *b_mc_stable_py;   //!
    TBranch        *b_mc_stable_pz;   //!
    TBranch        *b_mc_stable_m;   //!
    TBranch        *b_mc_stable_pdg;   //!
    TBranch        *b_mc_stable_charge;   //!
    TBranch        *b_mc_stable_isoverlay;   //!
    TBranch        *b_mc_stable_isisr;   //!
    TBranch        *b_mc_stable_y12;   //!
    TBranch        *b_mc_stable_y23;   //!
    TBranch        *b_mc_stable_d12;   //!
    TBranch        *b_mc_stable_d23;   //!
    TBranch        *b_mc_stable_jet_E;   //!
    TBranch        *b_mc_stable_jet_px;   //!
    TBranch        *b_mc_stable_jet_py;   //!
    TBranch        *b_mc_stable_jet_pz;   //!
    TBranch        *b_truejet_E;   //!
    TBranch        *b_truejet_px;   //!
    TBranch        *b_truejet_py;   //!
    TBranch        *b_truejet_pz;   //!
    TBranch        *b_truejet_type;   //!
    TBranch        *b_truejet_pdg;   //!
    TBranch        *b_jet_E;   //!
    TBranch        *b_jet_px;   //!
    TBranch        *b_jet_py;   //!
    TBranch        *b_jet_pz;   //!
    TBranch        *b_jet_btag;   //!
    TBranch        *b_jet_ctag;   //!
    TBranch        *b_y23;   //!
    TBranch        *b_y12;   //!
    TBranch        *b_d23;   //!
    TBranch        *b_d12;   //!
    TBranch        *b_oblateness;   //!
    TBranch        *b_aplanarity;   //!
    TBranch        *b_major_thrust_value;   //!
    TBranch        *b_major_thrust_axis;   //!
    TBranch        *b_minor_thrust_value;   //!
    TBranch        *b_minor_thrust_axis;   //!
    TBranch        *b_principle_thrust_value;   //!
    TBranch        *b_principle_thrust_axis;   //!
    TBranch        *b_sphericity;   //!
    TBranch        *b_sphericity_tensor;   //!
    TBranch        *b_pfo_n;   //!
    TBranch        *b_jet_nvtx;   //!
    TBranch        *b_pfo_n_j1;   //!
    TBranch        *b_jet_nvtx_j1;   //!
    TBranch        *b_pfo_n_j2;   //!
    TBranch        *b_jet_nvtx_j2;   //!
    TBranch        *b_pfo_match;   //!
    TBranch        *b_pfo_truejet_pdg;   //!
    TBranch        *b_pfo_truejet_type;   //!
    TBranch        *b_pfo_pdgcheat;   //!
    TBranch        *b_pfo_nparents;   //!
    TBranch        *b_pfo_pdgcheat_parent;   //!
    TBranch        *b_pfo_E;   //!
    TBranch        *b_pfo_E_calo;   //!
    TBranch        *b_pfo_px;   //!
    TBranch        *b_pfo_py;   //!
    TBranch        *b_pfo_pz;   //!
    TBranch        *b_pfo_m;   //!
    TBranch        *b_pfo_type;   //!
    TBranch        *b_pfo_isoverlay;   //!
    TBranch        *b_pfo_isisr;   //!
    TBranch        *b_pfo_vtx;   //!
    TBranch        *b_pfo_charge;   //!
    TBranch        *b_pfo_ntracks;   //!
    TBranch        *b_pfo_tpc_hits;   //!
    TBranch        *b_pfo_dedx;   //!
    TBranch        *b_pfo_dedxerror;   //!
    TBranch        *b_pfo_d0;   //!
    TBranch        *b_pfo_d0error;   //!
    TBranch        *b_pfo_z0;   //!
    TBranch        *b_pfo_z0error;   //!
    TBranch        *b_pfo_phi;   //!
    TBranch        *b_pfo_phierror;   //!
    TBranch        *b_pfo_omega;   //!
    TBranch        *b_pfo_omegaerror;   //!
    TBranch        *b_pfo_tanlambda;   //!
    TBranch        *b_pfo_tanlambdaerror;   //!
    TBranch        *b_pfo_chi2;   //!
    TBranch        *b_pfo_ndf;   //!
    TBranch        *b_pfo_vtxpt;   //!
    TBranch        *b_pfo_endpt;   //!
    TBranch        *b_pfo_pid;   //!
    TBranch        *b_pfo_pid_likelihood;   //!
    TBranch        *b_pfo_pid_eprob;   //!
    TBranch        *b_pfo_pid_muprob;   //!
    TBranch        *b_pfo_pid_piprob;   //!
    TBranch        *b_pfo_pid_kprob;   //!
    TBranch        *b_pfo_pid_pprob;   //!
    TBranch        *b_pfo_pid_hprob;   //!
    TBranch        *b_pfo_piddedx;   //!
    TBranch        *b_pfo_piddedx_likelihood;   //!
    TBranch        *b_pfo_piddedx_eprob;   //!
    TBranch        *b_pfo_piddedx_muprob;   //!
    TBranch        *b_pfo_piddedx_piprob;   //!
    TBranch        *b_pfo_piddedx_kprob;   //!
    TBranch        *b_pfo_piddedx_pprob;   //!
    TBranch        *b_pfo_piddedx_hprob;   //!
    TBranch        *b_pfo_piddedx_e_dedxdist;   //!
    TBranch        *b_pfo_piddedx_mu_dedxdist;   //!
    TBranch        *b_pfo_piddedx_pi_dedxdist;   //!
    TBranch        *b_pfo_piddedx_k_dedxdist;   //!
    TBranch        *b_pfo_piddedx_p_dedxdist;   //!
    TBranch        *b_pfo_piddedx_e_lkhood;   //!
    TBranch        *b_pfo_piddedx_mu_lkhood;   //!
    TBranch        *b_pfo_piddedx_pi_lkhood;   //!
    TBranch        *b_pfo_piddedx_k_lkhood;   //!
    TBranch        *b_pfo_piddedx_p_lkhood;   //!
    TBranch        *b_pfo_pidtof_p_at_calo;   //!
    TBranch        *b_pfo_pidtof_closest_beta_0ps;   //!
    TBranch        *b_pfo_pidtof_closest_beta_10ps;   //!
    TBranch        *b_pfo_pidtof_closest_beta_50ps;   //!
    TBranch        *b_pfo_pidtof_closest_beta_100ps;   //!
    TBranch        *b_pfo_pidtof_fastest_beta_0ps;   //!
    TBranch        *b_pfo_pidtof_fastest_beta_10ps;   //!
    TBranch        *b_pfo_pidtof_fastest_beta_50ps;   //!
    TBranch        *b_pfo_pidtof_fastest_beta_100ps;   //!
    TBranch        *b_pfo_pidtof_cylfit_beta_0ps;   //!
    TBranch        *b_pfo_pidtof_cylfit_beta_10ps;   //!
    TBranch        *b_pfo_pidtof_cylfit_beta_50ps;   //!
    TBranch        *b_pfo_pidtof_cylfit_beta_100ps;   //!
    TBranch        *b_pfo_pidtof_closestfit_beta_0ps;   //!
    TBranch        *b_pfo_pidtof_closestfit_beta_10ps;   //!
    TBranch        *b_pfo_pidtof_closestfit_beta_50ps;   //!
    TBranch        *b_pfo_pidtof_closestfit_beta_100ps;   //!

};

struct Tree_SSbar {
  public:
    Int_t   lpfo_match[NJETS] = {0};
    Int_t   lpfo_truejet_pdg[NJETS] = {0};
    Int_t   lpfo_truejet_type[NJETS] = {0};
    Int_t   lpfo_pdgcheat[NJETS] = {0};
    Int_t   lpfo_nparents[NJETS] = {0};
    Int_t   lpfo_pdgcheat_parent[NJETS][MAX_NPARTICLES] = {{0}};
    Float_t lpfo_E[NJETS] = {0};
    Float_t lpfo_px[NJETS] = {0};
    Float_t lpfo_py[NJETS] = {0};
    Float_t lpfo_pz[NJETS] = {0};
    Float_t lpfo_m[NJETS] = {0};
    Int_t   lpfo_type[NJETS] = {0};
    Int_t   lpfo_isoverlay[NJETS] = {0};
    Int_t   lpfo_isisr[NJETS] = {0};
    Int_t   lpfo_vtx[NJETS] = {0};
    Int_t   lpfo_charge[NJETS] = {0};
    Int_t   lpfo_ntracks[NJETS] = {0};
    Int_t   lpfo_tpc_hits[NJETS] = {0};
    Float_t lpfo_dedx[NJETS] = {0};
    Float_t lpfo_dedxerror[NJETS] = {0};
    Float_t lpfo_d0[NJETS] = {0};
    Float_t lpfo_d0error[NJETS] = {0};
    Float_t lpfo_z0[NJETS] = {0};
    Float_t lpfo_z0error[NJETS] = {0};
    Float_t lpfo_phi[NJETS] = {0};
    Float_t lpfo_phierror[NJETS] = {0};
    Float_t lpfo_omega[NJETS] = {0};
    Float_t lpfo_omegaerror[NJETS] = {0};
    Float_t lpfo_tanlambda[NJETS] = {0};
    Float_t lpfo_tanlambdaerror[NJETS] = {0};
    Float_t lpfo_chi2[NJETS] = {0};
    Float_t lpfo_ndf[NJETS] = {0};
    Float_t lpfo_vtxpt[NJETS][3] = {{0}};
    Float_t lpfo_endpt[NJETS][3] = {{0}};
    Int_t   lpfo_pid[NJETS] = {0};
    Float_t lpfo_pid_likelihood[NJETS] = {0};
    Float_t lpfo_pid_eprob[NJETS] = {0};
    Float_t lpfo_pid_muprob[NJETS] = {0};
    Float_t lpfo_pid_piprob[NJETS] = {0};
    Float_t lpfo_pid_kprob[NJETS] = {0};
    Float_t lpfo_pid_pprob[NJETS] = {0};
    Float_t lpfo_pid_hprob[NJETS] = {0};
    Int_t   lpfo_piddedx[NJETS] = {0};
    Float_t lpfo_piddedx_likelihood[NJETS] = {0};
    Float_t lpfo_piddedx_eprob[NJETS] = {0};
    Float_t lpfo_piddedx_muprob[NJETS] = {0};
    Float_t lpfo_piddedx_piprob[NJETS] = {0};
    Float_t lpfo_piddedx_kprob[NJETS] = {0};
    Float_t lpfo_piddedx_pprob[NJETS] = {0};
    Float_t lpfo_piddedx_hprob[NJETS] = {0};
    Float_t lpfo_piddedx_e_dedxdist[NJETS] = {0};
    Float_t lpfo_piddedx_mu_dedxdist[NJETS] = {0};
    Float_t lpfo_piddedx_pi_dedxdist[NJETS] = {0};
    Float_t lpfo_piddedx_k_dedxdist[NJETS] = {0};
    Float_t lpfo_piddedx_p_dedxdist[NJETS] = {0};
    Float_t lpfo_piddedx_e_lkhood[NJETS] = {0};
    Float_t lpfo_piddedx_mu_lkhood[NJETS] = {0};
    Float_t lpfo_piddedx_pi_lkhood[NJETS] = {0};
    Float_t lpfo_piddedx_k_lkhood[NJETS] = {0};
    Float_t lpfo_piddedx_p_lkhood[NJETS] = {0};
    Float_t lpfo_pidtof_p_at_calo[NJETS] = {0};
    Float_t lpfo_pidtof_closest_beta_0ps[NJETS] = {0};
    Float_t lpfo_pidtof_closest_beta_10ps[NJETS] = {0};
    Float_t lpfo_pidtof_closest_beta_50ps[NJETS] = {0};
    Float_t lpfo_pidtof_closest_beta_100ps[NJETS] = {0};
    Float_t lpfo_pidtof_fastest_beta_0ps[NJETS] = {0};
    Float_t lpfo_pidtof_fastest_beta_10ps[NJETS] = {0};
    Float_t lpfo_pidtof_fastest_beta_50ps[NJETS] = {0};
    Float_t lpfo_pidtof_fastest_beta_100ps[NJETS] = {0};
    Float_t lpfo_pidtof_cylfit_beta_0ps[NJETS] = {0};
    Float_t lpfo_pidtof_cylfit_beta_10ps[NJETS] = {0};
    Float_t lpfo_pidtof_cylfit_beta_50ps[NJETS] = {0};
    Float_t lpfo_pidtof_cylfit_beta_100ps[NJETS] = {0};
    Float_t lpfo_pidtof_closestfit_beta_0ps[NJETS] = {0};
    Float_t lpfo_pidtof_closestfit_beta_10ps[NJETS] = {0};
    Float_t lpfo_pidtof_closestfit_beta_50ps[NJETS] = {0};
    Float_t lpfo_pidtof_closestfit_beta_100ps[NJETS] = {0};

};

#endif