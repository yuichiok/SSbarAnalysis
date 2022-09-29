#ifndef GUARD_TreeVariables_h
#define GUARD_TreeVariables_h

/*------------------------------------------------------------------------------
   TreeVariables
 Created : 2022-09-14  okugawa
 Main class of VectorTool program.
------------------------------------------------------------------------------*/
#include <TTree.h>

struct TreeVariables
{
  
const static int MAX_NPARTICLES = 1000;
const static int NTRUE_JETS     = 5;
const static int NQUARKS        = 2;
const static int NJETS          = 2;


  public:
    Int_t   lpfo_config = 0;
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

    ClassDef(TreeVariables,1)

};

#endif