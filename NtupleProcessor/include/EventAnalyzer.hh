#ifndef GUARD_EventAnalyzer_h
#define GUARD_EventAnalyzer_h

/*------------------------------------------------------------------------------
   NtupleProcessor
 Created : 2022-09-02  okugawa
 Main class of NtupleProcessor program. Created to handle input and running of
 ntuple processing.
 Takes input variables (datasets, configurations, etc) and sets up the
 appropriate classes to handle each portion of the analysis process.
------------------------------------------------------------------------------*/
#include <TString.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <fstream>

class EventAnalyzer
{
  public:
    EventAnalyzer(TString o="");
    virtual ~EventAnalyzer(){}

  // methods
    bool mapTree(TTree*);                         // Maps class variables to an input TTree.
    void evalCriteria();                          // Evaluates the class' list of event selection criteria
    virtual Bool_t   Notify();


  // Running Variables
    TString options;         // Options input with TreeIterator.

    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
    Float_t         mc_quark_E[2];
    Float_t         mc_quark_px[2];
    Float_t         mc_quark_py[2];
    Float_t         mc_quark_pz[2];
    Float_t         mc_quark_m[2];
    Float_t         mc_quark_pdg[2];
    Float_t         mc_quark_charge[2];
    Float_t         mc_ISR_E[2];
    Float_t         mc_ISR_px[2];
    Float_t         mc_ISR_py[2];
    Float_t         mc_ISR_pz[2];
    Float_t         mc_ISR_m[2];
    Float_t         mc_ISR_pdg[2];
    Float_t         mc_ISR_charge[2];
    Int_t           mc_quark_ps_n;
    Float_t         mc_quark_ps_E[1000];   //[mc_quark_ps_n]
    Float_t         mc_quark_ps_px[1000];   //[mc_quark_ps_n]
    Float_t         mc_quark_ps_py[1000];   //[mc_quark_ps_n]
    Float_t         mc_quark_ps_pz[1000];   //[mc_quark_ps_n]
    Float_t         mc_quark_ps_m[1000];   //[mc_quark_ps_n]
    Float_t         mc_quark_ps_pdg[1000];   //[mc_quark_ps_n]
    Float_t         mc_quark_ps_charge[1000];   //[mc_quark_ps_n]
    Float_t         mc_quark_ps_y12;
    Float_t         mc_quark_ps_y23;
    Float_t         mc_quark_ps_d12;
    Float_t         mc_quark_ps_d23;
    Float_t         mc_quark_ps_jet_E[2];
    Float_t         mc_quark_ps_jet_px[2];
    Float_t         mc_quark_ps_jet_py[2];
    Float_t         mc_quark_ps_jet_pz[2];
    Int_t           mc_stable_n;
    Float_t         mc_stable_E[1000];   //[mc_stable_n]
    Float_t         mc_stable_px[1000];   //[mc_stable_n]
    Float_t         mc_stable_py[1000];   //[mc_stable_n]
    Float_t         mc_stable_pz[1000];   //[mc_stable_n]
    Float_t         mc_stable_m[1000];   //[mc_stable_n]
    Int_t           mc_stable_pdg[1000];   //[mc_stable_n]
    Float_t         mc_stable_charge[1000];   //[mc_stable_n]
    Int_t           mc_stable_isoverlay[1000];   //[mc_stable_n]
    Int_t           mc_stable_isisr[1000];   //[mc_stable_n]
    Float_t         mc_stable_y12;
    Float_t         mc_stable_y23;
    Float_t         mc_stable_d12;
    Float_t         mc_stable_d23;
    Float_t         mc_stable_jet_E[2];
    Float_t         mc_stable_jet_px[2];
    Float_t         mc_stable_jet_py[2];
    Float_t         mc_stable_jet_pz[2];
    Float_t         truejet_E[5];
    Float_t         truejet_px[5];
    Float_t         truejet_py[5];
    Float_t         truejet_pz[5];
    Int_t           truejet_type[5];
    Int_t           truejet_pdg[5];
    Float_t         jet_E[2];
    Float_t         jet_px[2];
    Float_t         jet_py[2];
    Float_t         jet_pz[2];
    Float_t         jet_btag[2];
    Float_t         jet_ctag[2];
    Float_t         y23;
    Float_t         y12;
    Float_t         d23;
    Float_t         d12;
    Float_t         oblateness;
    Float_t         aplanarity;
    Float_t         major_thrust_value;
    Float_t         major_thrust_axis[3];
    Float_t         minor_thrust_value;
    Float_t         minor_thrust_axis[3];
    Float_t         principle_thrust_value;
    Float_t         principle_thrust_axis[3];
    Float_t         sphericity;
    Float_t         sphericity_tensor[3];
    Int_t           pfo_n;
    Int_t           jet_nvtx;
    Int_t           pfo_n_j1;
    Int_t           jet_nvtx_j1;
    Int_t           pfo_n_j2;
    Int_t           jet_nvtx_j2;
    Int_t           pfo_match[1000];   //[pfo_n]
    Int_t           pfo_truejet_pdg[1000];   //[pfo_n]
    Int_t           pfo_truejet_type[1000];   //[pfo_n]
    Int_t           pfo_pdgcheat[1000];   //[pfo_n]
    Int_t           pfo_nparents[1000];   //[pfo_n]
    Int_t           pfo_pdgcheat_parent[1000][1000];   //[pfo_n]
    Float_t         pfo_E[1000];   //[pfo_n]
    Float_t         pfo_E_calo[1000];   //[pfo_n]
    Float_t         pfo_px[1000];   //[pfo_n]
    Float_t         pfo_py[1000];   //[pfo_n]
    Float_t         pfo_pz[1000];   //[pfo_n]
    Float_t         pfo_m[1000];   //[pfo_n]
    Int_t           pfo_type[1000];   //[pfo_n]
    Int_t           pfo_isoverlay[1000];   //[pfo_n]
    Int_t           pfo_isisr[1000];   //[pfo_n]
    Int_t           pfo_vtx[1000];   //[pfo_n]
    Int_t           pfo_charge[1000];   //[pfo_n]
    Int_t           pfo_ntracks[1000];   //[pfo_n]
    Int_t           pfo_tpc_hits[1000];   //[pfo_n]
    Float_t         pfo_dedx[1000];   //[pfo_n]
    Float_t         pfo_dedxerror[1000];   //[pfo_n]
    Float_t         pfo_d0[1000];   //[pfo_n]
    Float_t         pfo_d0error[1000];   //[pfo_n]
    Float_t         pfo_z0[1000];   //[pfo_n]
    Float_t         pfo_z0error[1000];   //[pfo_n]
    Float_t         pfo_phi[1000];   //[pfo_n]
    Float_t         pfo_phierror[1000];   //[pfo_n]
    Float_t         pfo_omega[1000];   //[pfo_n]
    Float_t         pfo_omegaerror[1000];   //[pfo_n]
    Float_t         pfo_tanlambda[1000];   //[pfo_n]
    Float_t         pfo_tanlambdaerror[1000];   //[pfo_n]
    Float_t         pfo_chi2[1000];   //[pfo_n]
    Float_t         pfo_ndf[1000];   //[pfo_n]
    Float_t         pfo_vtxpt[1000][3];   //[pfo_n]
    Float_t         pfo_endpt[1000][3];   //[pfo_n]
    Int_t           pfo_pid[1000];   //[pfo_n]
    Float_t         pfo_pid_likelihood[1000];   //[pfo_n]
    Float_t         pfo_pid_eprob[1000];   //[pfo_n]
    Float_t         pfo_pid_muprob[1000];   //[pfo_n]
    Float_t         pfo_pid_piprob[1000];   //[pfo_n]
    Float_t         pfo_pid_kprob[1000];   //[pfo_n]
    Float_t         pfo_pid_pprob[1000];   //[pfo_n]
    Float_t         pfo_pid_hprob[1000];   //[pfo_n]
    Int_t           pfo_piddedx[1000];   //[pfo_n]
    Float_t         pfo_piddedx_likelihood[1000];   //[pfo_n]
    Float_t         pfo_piddedx_eprob[1000];   //[pfo_n]
    Float_t         pfo_piddedx_muprob[1000];   //[pfo_n]
    Float_t         pfo_piddedx_piprob[1000];   //[pfo_n]
    Float_t         pfo_piddedx_kprob[1000];   //[pfo_n]
    Float_t         pfo_piddedx_pprob[1000];   //[pfo_n]
    Float_t         pfo_piddedx_hprob[1000];   //[pfo_n]
    Float_t         pfo_piddedx_e_dedxdist[1000];   //[pfo_n]
    Float_t         pfo_piddedx_mu_dedxdist[1000];   //[pfo_n]
    Float_t         pfo_piddedx_pi_dedxdist[1000];   //[pfo_n]
    Float_t         pfo_piddedx_k_dedxdist[1000];   //[pfo_n]
    Float_t         pfo_piddedx_p_dedxdist[1000];   //[pfo_n]
    Float_t         pfo_piddedx_e_lkhood[1000];   //[pfo_n]
    Float_t         pfo_piddedx_mu_lkhood[1000];   //[pfo_n]
    Float_t         pfo_piddedx_pi_lkhood[1000];   //[pfo_n]
    Float_t         pfo_piddedx_k_lkhood[1000];   //[pfo_n]
    Float_t         pfo_piddedx_p_lkhood[1000];   //[pfo_n]
    Float_t         pfo_pidtof_p_at_calo[1000];   //[pfo_n]
    Float_t         pfo_pidtof_closest_beta_0ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_closest_beta_10ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_closest_beta_50ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_closest_beta_100ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_fastest_beta_0ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_fastest_beta_10ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_fastest_beta_50ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_fastest_beta_100ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_cylfit_beta_0ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_cylfit_beta_10ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_cylfit_beta_50ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_cylfit_beta_100ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_closestfit_beta_0ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_closestfit_beta_10ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_closestfit_beta_50ps[1000];   //[pfo_n]
    Float_t         pfo_pidtof_closestfit_beta_100ps[1000];   //[pfo_n]

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

  // Extra data
    long patEventsAnalyzed;     // Number of events that were processed to make the Ntuple.
    long entriesInNtuple  ;     // Number of events that were processed to make the Ntuple.

  private: 

};

#endif