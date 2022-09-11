#ifndef GUARD_TreeStructures_h
#define GUARD_TreeStructures_h

/*------------------------------------------------------------------------------
   TreeReader
 Created : 2022-09-09  okugawa
 Main class of TreeReader program.
------------------------------------------------------------------------------*/

const static int MAX_NPARTICLES = 1000;
const static int NTRUE_JETS     = 5;
const static int NQUARKS        = 2;
const static int NJETS          = 2;

struct MC_QQbar  {

  public:
  // MC VARIABLES
    Float_t mc_quark_E[NQUARKS];        Int_t   mc_quark_ps_n;                        Int_t   mc_stable_n;
    Float_t mc_quark_px[NQUARKS];       Float_t mc_quark_ps_E[MAX_NPARTICLES];        Float_t mc_stable_E[MAX_NPARTICLES];
    Float_t mc_quark_py[NQUARKS];       Float_t mc_quark_ps_px[MAX_NPARTICLES];       Float_t mc_stable_px[MAX_NPARTICLES];
    Float_t mc_quark_pz[NQUARKS];       Float_t mc_quark_ps_py[MAX_NPARTICLES];       Float_t mc_stable_py[MAX_NPARTICLES];
    Float_t mc_quark_m[NQUARKS];        Float_t mc_quark_ps_pz[MAX_NPARTICLES];       Float_t mc_stable_pz[MAX_NPARTICLES];
    Float_t mc_quark_pdg[NQUARKS];      Float_t mc_quark_ps_m[MAX_NPARTICLES];        Float_t mc_stable_m[MAX_NPARTICLES];
    Float_t mc_quark_charge[NQUARKS];   Float_t mc_quark_ps_pdg[MAX_NPARTICLES];      Int_t   mc_stable_pdg[MAX_NPARTICLES];
    Float_t mc_ISR_E[NQUARKS];          Float_t mc_quark_ps_charge[MAX_NPARTICLES];   Float_t mc_stable_charge[MAX_NPARTICLES]; 
    Float_t mc_ISR_px[NQUARKS];         Float_t mc_quark_ps_y12;                      Int_t   mc_stable_isoverlay[MAX_NPARTICLES];
    Float_t mc_ISR_py[NQUARKS];         Float_t mc_quark_ps_y23;                      Int_t   mc_stable_isisr[MAX_NPARTICLES];
    Float_t mc_ISR_pz[NQUARKS];         Float_t mc_quark_ps_d12;                      Float_t mc_stable_y12;
    Float_t mc_ISR_m[NQUARKS];          Float_t mc_quark_ps_d23;                      Float_t mc_stable_y23;
    Float_t mc_ISR_pdg[NQUARKS];        Float_t mc_quark_ps_jet_E[NJETS];             Float_t mc_stable_d12;
    Float_t mc_ISR_charge[NQUARKS];     Float_t mc_quark_ps_jet_px[NJETS];            Float_t mc_stable_d23;
                                        Float_t mc_quark_ps_jet_py[NJETS];            Float_t mc_stable_jet_E[NJETS];
                                        Float_t mc_quark_ps_jet_pz[NJETS];            Float_t mc_stable_jet_px[NJETS];
                                                                                      Float_t mc_stable_jet_py[NJETS];
                                                                                      Float_t mc_stable_jet_pz[NJETS];
                                                                                      
};

struct Jet_QQbar  {

  public:
  // JET VARIABLES
    Float_t truejet_E[NTRUE_JETS];      Float_t jet_E[NJETS];       Float_t y23;          Float_t major_thrust_value;       Float_t major_thrust_axis[3];
    Float_t truejet_px[NTRUE_JETS];     Float_t jet_px[NJETS];      Float_t y12;          Float_t minor_thrust_value;       Float_t minor_thrust_axis[3];
    Float_t truejet_py[NTRUE_JETS];     Float_t jet_py[NJETS];      Float_t d23;          Float_t principle_thrust_value;   Float_t principle_thrust_axis[3];
    Float_t truejet_pz[NTRUE_JETS];     Float_t jet_pz[NJETS];      Float_t d12;          Float_t sphericity;               Float_t sphericity_tensor[3];
    Int_t   truejet_type[NTRUE_JETS];   Float_t jet_btag[NJETS];    Float_t oblateness; 
    Int_t   truejet_pdg[NTRUE_JETS];    Float_t jet_ctag[NJETS];    Float_t aplanarity;
    
};

struct PFO_QQbar  {

  public:
  // PFO VARIABLES
    Int_t pfo_n;
    Int_t jet_nvtx;
    Int_t pfo_n_j1;
    Int_t jet_nvtx_j1;
    Int_t pfo_n_j2;
    Int_t jet_nvtx_j2;
    Int_t pfo_match[MAX_NPARTICLES];
    Int_t pfo_truejet_pdg[MAX_NPARTICLES];
    Int_t pfo_truejet_type[MAX_NPARTICLES];
    Int_t pfo_pdgcheat[MAX_NPARTICLES];
    Int_t pfo_nparents[MAX_NPARTICLES];
    Int_t pfo_pdgcheat_parent[MAX_NPARTICLES][MAX_NPARTICLES];
    Float_t pfo_E[MAX_NPARTICLES];
    Float_t pfo_E_calo[MAX_NPARTICLES];
    Float_t pfo_px[MAX_NPARTICLES];
    Float_t pfo_py[MAX_NPARTICLES];
    Float_t pfo_pz[MAX_NPARTICLES];
    Float_t pfo_m[MAX_NPARTICLES];
    Int_t pfo_type[MAX_NPARTICLES];
    Int_t pfo_isoverlay[MAX_NPARTICLES];
    Int_t pfo_isisr[MAX_NPARTICLES];
    Int_t pfo_vtx[MAX_NPARTICLES];
    Int_t pfo_charge[MAX_NPARTICLES];
    Int_t pfo_ntracks[MAX_NPARTICLES];
    Int_t pfo_tpc_hits[MAX_NPARTICLES];
    Float_t pfo_dedx[MAX_NPARTICLES];
    Float_t pfo_dedxerror[MAX_NPARTICLES];
    Float_t pfo_d0[MAX_NPARTICLES];
    Float_t pfo_d0error[MAX_NPARTICLES];
    Float_t pfo_z0[MAX_NPARTICLES];
    Float_t pfo_z0error[MAX_NPARTICLES];
    Float_t pfo_phi[MAX_NPARTICLES];
    Float_t pfo_phierror[MAX_NPARTICLES];
    Float_t pfo_omega[MAX_NPARTICLES];
    Float_t pfo_omegaerror[MAX_NPARTICLES];
    Float_t pfo_tanlambda[MAX_NPARTICLES];
    Float_t pfo_tanlambdaerror[MAX_NPARTICLES];
    Float_t pfo_chi2[MAX_NPARTICLES];
    Float_t pfo_ndf[MAX_NPARTICLES];
    Float_t pfo_vtxpt[MAX_NPARTICLES][3];
    Float_t pfo_endpt[MAX_NPARTICLES][3];
    Int_t pfo_pid[MAX_NPARTICLES];
    Float_t pfo_pid_likelihood[MAX_NPARTICLES];
    Float_t pfo_pid_eprob[MAX_NPARTICLES];
    Float_t pfo_pid_muprob[MAX_NPARTICLES];
    Float_t pfo_pid_piprob[MAX_NPARTICLES];
    Float_t pfo_pid_kprob[MAX_NPARTICLES];
    Float_t pfo_pid_pprob[MAX_NPARTICLES];
    Float_t pfo_pid_hprob[MAX_NPARTICLES];
    Int_t pfo_piddedx[MAX_NPARTICLES];
    Float_t pfo_piddedx_likelihood[MAX_NPARTICLES];
    Float_t pfo_piddedx_eprob[MAX_NPARTICLES];
    Float_t pfo_piddedx_muprob[MAX_NPARTICLES];
    Float_t pfo_piddedx_piprob[MAX_NPARTICLES];
    Float_t pfo_piddedx_kprob[MAX_NPARTICLES];
    Float_t pfo_piddedx_pprob[MAX_NPARTICLES];
    Float_t pfo_piddedx_hprob[MAX_NPARTICLES];
    Float_t pfo_piddedx_e_dedxdist[MAX_NPARTICLES];
    Float_t pfo_piddedx_mu_dedxdist[MAX_NPARTICLES];
    Float_t pfo_piddedx_pi_dedxdist[MAX_NPARTICLES];
    Float_t pfo_piddedx_k_dedxdist[MAX_NPARTICLES];
    Float_t pfo_piddedx_p_dedxdist[MAX_NPARTICLES];
    Float_t pfo_piddedx_e_lkhood[MAX_NPARTICLES];
    Float_t pfo_piddedx_mu_lkhood[MAX_NPARTICLES];
    Float_t pfo_piddedx_pi_lkhood[MAX_NPARTICLES];
    Float_t pfo_piddedx_k_lkhood[MAX_NPARTICLES];
    Float_t pfo_piddedx_p_lkhood[MAX_NPARTICLES];
    Float_t pfo_pidtof_p_at_calo[MAX_NPARTICLES];
    Float_t pfo_pidtof_closest_beta_0ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_closest_beta_10ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_closest_beta_50ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_closest_beta_100ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_fastest_beta_0ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_fastest_beta_10ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_fastest_beta_50ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_fastest_beta_100ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_cylfit_beta_0ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_cylfit_beta_10ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_cylfit_beta_50ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_cylfit_beta_100ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_closestfit_beta_0ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_closestfit_beta_10ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_closestfit_beta_50ps[MAX_NPARTICLES];
    Float_t pfo_pidtof_closestfit_beta_100ps[MAX_NPARTICLES];

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

#endif