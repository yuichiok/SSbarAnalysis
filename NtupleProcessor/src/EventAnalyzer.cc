
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
#include "../include/EventAnalyzer.hh"
#include "../include/TreeReader.hh"

using std::cout;   using std::endl;

EventAnalyzer::EventAnalyzer(TString o)
: options(o)
{
  // TEST output
    cout << "    NtupleProcessor: Created." << endl;

    patEventsAnalyzed = 0;
    entriesInNtuple   = 0;

}

bool EventAnalyzer::mapTree(TTree* tree)
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

    fChain->SetBranchAddress("mc_quark_E", _stats.mc_quark_E, &_branch.b_mc_quark_E);
    fChain->SetBranchAddress("mc_quark_px", _stats.mc_quark_px, &_branch.b_mc_quark_px);
    fChain->SetBranchAddress("mc_quark_py", _stats.mc_quark_py, &_branch.b_mc_quark_py);
    fChain->SetBranchAddress("mc_quark_pz", _stats.mc_quark_pz, &_branch.b_mc_quark_pz);
    fChain->SetBranchAddress("mc_quark_m", _stats.mc_quark_m, &_branch.b_mc_quark_m);
    fChain->SetBranchAddress("mc_quark_pdg", _stats.mc_quark_pdg, &_branch.b_mc_quark_pdg);
    fChain->SetBranchAddress("mc_quark_charge", _stats.mc_quark_charge, &_branch.b_mc_quark_charge);
    fChain->SetBranchAddress("mc_ISR_E", _stats.mc_ISR_E, &_branch.b_mc_ISR_E);
    fChain->SetBranchAddress("mc_ISR_px", _stats.mc_ISR_px, &_branch.b_mc_ISR_px);
    fChain->SetBranchAddress("mc_ISR_py", _stats.mc_ISR_py, &_branch.b_mc_ISR_py);
    fChain->SetBranchAddress("mc_ISR_pz", _stats.mc_ISR_pz, &_branch.b_mc_ISR_pz);
    fChain->SetBranchAddress("mc_ISR_m", _stats.mc_ISR_m, &_branch.b_mc_ISR_m);
    fChain->SetBranchAddress("mc_ISR_pdg", _stats.mc_ISR_pdg, &_branch.b_mc_ISR_pdg);
    fChain->SetBranchAddress("mc_ISR_charge", _stats.mc_ISR_charge, &_branch.b_mc_ISR_charge);
    fChain->SetBranchAddress("mc_quark_ps_n", &_stats.mc_quark_ps_n, &_branch.b_mc_quark_ps_n);
    fChain->SetBranchAddress("mc_quark_ps_E", _stats.mc_quark_ps_E, &_branch.b_mc_quark_ps_E);
    fChain->SetBranchAddress("mc_quark_ps_px", _stats.mc_quark_ps_px, &_branch.b_mc_quark_ps_px);
    fChain->SetBranchAddress("mc_quark_ps_py", _stats.mc_quark_ps_py, &_branch.b_mc_quark_ps_py);
    fChain->SetBranchAddress("mc_quark_ps_pz", _stats.mc_quark_ps_pz, &_branch.b_mc_quark_ps_pz);
    fChain->SetBranchAddress("mc_quark_ps_m", _stats.mc_quark_ps_m, &_branch.b_mc_quark_ps_m);
    fChain->SetBranchAddress("mc_quark_ps_pdg", _stats.mc_quark_ps_pdg, &_branch.b_mc_quark_ps_pdg);
    fChain->SetBranchAddress("mc_quark_ps_charge", _stats.mc_quark_ps_charge, &_branch.b_mc_quark_ps_charge);
    fChain->SetBranchAddress("mc_quark_ps_y12", &_stats.mc_quark_ps_y12, &_branch.b_mc_quark_ps_y12);
    fChain->SetBranchAddress("mc_quark_ps_y23", &_stats.mc_quark_ps_y23, &_branch.b_mc_quark_ps_y23);
    fChain->SetBranchAddress("mc_quark_ps_d12", &_stats.mc_quark_ps_d12, &_branch.b_mc_quark_ps_d12);
    fChain->SetBranchAddress("mc_quark_ps_d23", &_stats.mc_quark_ps_d23, &_branch.b_mc_quark_ps_d23);
    fChain->SetBranchAddress("mc_quark_ps_jet_E", _stats.mc_quark_ps_jet_E, &_branch.b_mc_quark_ps_jet_E);
    fChain->SetBranchAddress("mc_quark_ps_jet_px", _stats.mc_quark_ps_jet_px, &_branch.b_mc_quark_ps_jet_px);
    fChain->SetBranchAddress("mc_quark_ps_jet_py", _stats.mc_quark_ps_jet_py, &_branch.b_mc_quark_ps_jet_py);
    fChain->SetBranchAddress("mc_quark_ps_jet_pz", _stats.mc_quark_ps_jet_pz, &_branch.b_mc_quark_ps_jet_pz);
    fChain->SetBranchAddress("mc_stable_n", &_stats.mc_stable_n, &_branch.b_mc_stable_n);
    fChain->SetBranchAddress("mc_stable_E", _stats.mc_stable_E, &_branch.b_mc_stable_E);
    fChain->SetBranchAddress("mc_stable_px", _stats.mc_stable_px, &_branch.b_mc_stable_px);
    fChain->SetBranchAddress("mc_stable_py", _stats.mc_stable_py, &_branch.b_mc_stable_py);
    fChain->SetBranchAddress("mc_stable_pz", _stats.mc_stable_pz, &_branch.b_mc_stable_pz);
    fChain->SetBranchAddress("mc_stable_m", _stats.mc_stable_m, &_branch.b_mc_stable_m);
    fChain->SetBranchAddress("mc_stable_pdg", _stats.mc_stable_pdg, &_branch.b_mc_stable_pdg);
    fChain->SetBranchAddress("mc_stable_charge", _stats.mc_stable_charge, &_branch.b_mc_stable_charge);
    fChain->SetBranchAddress("mc_stable_isoverlay", _stats.mc_stable_isoverlay, &_branch.b_mc_stable_isoverlay);
    fChain->SetBranchAddress("mc_stable_isisr", _stats.mc_stable_isisr, &_branch.b_mc_stable_isisr);
    fChain->SetBranchAddress("mc_stable_y12", &_stats.mc_stable_y12, &_branch.b_mc_stable_y12);
    fChain->SetBranchAddress("mc_stable_y23", &_stats.mc_stable_y23, &_branch.b_mc_stable_y23);
    fChain->SetBranchAddress("mc_stable_d12", &_stats.mc_stable_d12, &_branch.b_mc_stable_d12);
    fChain->SetBranchAddress("mc_stable_d23", &_stats.mc_stable_d23, &_branch.b_mc_stable_d23);
    fChain->SetBranchAddress("mc_stable_jet_E", _stats.mc_stable_jet_E, &_branch.b_mc_stable_jet_E);
    fChain->SetBranchAddress("mc_stable_jet_px", _stats.mc_stable_jet_px, &_branch.b_mc_stable_jet_px);
    fChain->SetBranchAddress("mc_stable_jet_py", _stats.mc_stable_jet_py, &_branch.b_mc_stable_jet_py);
    fChain->SetBranchAddress("mc_stable_jet_pz", _stats.mc_stable_jet_pz, &_branch.b_mc_stable_jet_pz);
    fChain->SetBranchAddress("truejet_E", _stats.truejet_E, &_branch.b_truejet_E);
    fChain->SetBranchAddress("truejet_px", _stats.truejet_px, &_branch.b_truejet_px);
    fChain->SetBranchAddress("truejet_py", _stats.truejet_py, &_branch.b_truejet_py);
    fChain->SetBranchAddress("truejet_pz", _stats.truejet_pz, &_branch.b_truejet_pz);
    fChain->SetBranchAddress("truejet_type", _stats.truejet_type, &_branch.b_truejet_type);
    fChain->SetBranchAddress("truejet_pdg", _stats.truejet_pdg, &_branch.b_truejet_pdg);
    fChain->SetBranchAddress("jet_E", _stats.jet_E, &_branch.b_jet_E);
    fChain->SetBranchAddress("jet_px", _stats.jet_px, &_branch.b_jet_px);
    fChain->SetBranchAddress("jet_py", _stats.jet_py, &_branch.b_jet_py);
    fChain->SetBranchAddress("jet_pz", _stats.jet_pz, &_branch.b_jet_pz);
    fChain->SetBranchAddress("jet_btag", _stats.jet_btag, &_branch.b_jet_btag);
    fChain->SetBranchAddress("jet_ctag", _stats.jet_ctag, &_branch.b_jet_ctag);
    fChain->SetBranchAddress("y23", &_stats.y23, &_branch.b_y23);
    fChain->SetBranchAddress("y12", &_stats.y12, &_branch.b_y12);
    fChain->SetBranchAddress("d23", &_stats.d23, &_branch.b_d23);
    fChain->SetBranchAddress("d12", &_stats.d12, &_branch.b_d12);
    fChain->SetBranchAddress("oblateness", &_stats.oblateness, &_branch.b_oblateness);
    fChain->SetBranchAddress("aplanarity", &_stats.aplanarity, &_branch.b_aplanarity);
    fChain->SetBranchAddress("major_thrust_value", &_stats.major_thrust_value, &_branch.b_major_thrust_value);
    fChain->SetBranchAddress("major_thrust_axis", _stats.major_thrust_axis, &_branch.b_major_thrust_axis);
    fChain->SetBranchAddress("minor_thrust_value", &_stats.minor_thrust_value, &_branch.b_minor_thrust_value);
    fChain->SetBranchAddress("minor_thrust_axis", _stats.minor_thrust_axis, &_branch.b_minor_thrust_axis);
    fChain->SetBranchAddress("principle_thrust_value", &_stats.principle_thrust_value, &_branch.b_principle_thrust_value);
    fChain->SetBranchAddress("principle_thrust_axis", _stats.principle_thrust_axis, &_branch.b_principle_thrust_axis);
    fChain->SetBranchAddress("sphericity", &_stats.sphericity, &_branch.b_sphericity);
    fChain->SetBranchAddress("sphericity_tensor", _stats.sphericity_tensor, &_branch.b_sphericity_tensor);
    fChain->SetBranchAddress("pfo_n", &_stats.pfo_n, &_branch.b_pfo_n);
    fChain->SetBranchAddress("jet_nvtx", &_stats.jet_nvtx, &_branch.b_jet_nvtx);
    fChain->SetBranchAddress("pfo_n_j1", &_stats.pfo_n_j1, &_branch.b_pfo_n_j1);
    fChain->SetBranchAddress("jet_nvtx_j1", &_stats.jet_nvtx_j1, &_branch.b_jet_nvtx_j1);
    fChain->SetBranchAddress("pfo_n_j2", &_stats.pfo_n_j2, &_branch.b_pfo_n_j2);
    fChain->SetBranchAddress("jet_nvtx_j2", &_stats.jet_nvtx_j2, &_branch.b_jet_nvtx_j2);
    fChain->SetBranchAddress("pfo_match", _stats.pfo_match, &_branch.b_pfo_match);
    fChain->SetBranchAddress("pfo_truejet_pdg", _stats.pfo_truejet_pdg, &_branch.b_pfo_truejet_pdg);
    fChain->SetBranchAddress("pfo_truejet_type", _stats.pfo_truejet_type, &_branch.b_pfo_truejet_type);
    fChain->SetBranchAddress("pfo_pdgcheat", _stats.pfo_pdgcheat, &_branch.b_pfo_pdgcheat);
    fChain->SetBranchAddress("pfo_nparents", _stats.pfo_nparents, &_branch.b_pfo_nparents);
    fChain->SetBranchAddress("pfo_pdgcheat_parent", _stats.pfo_pdgcheat_parent, &_branch.b_pfo_pdgcheat_parent);
    fChain->SetBranchAddress("pfo_E", _stats.pfo_E, &_branch.b_pfo_E);
    fChain->SetBranchAddress("pfo_E_calo", _stats.pfo_E_calo, &_branch.b_pfo_E_calo);
    fChain->SetBranchAddress("pfo_px", _stats.pfo_px, &_branch.b_pfo_px);
    fChain->SetBranchAddress("pfo_py", _stats.pfo_py, &_branch.b_pfo_py);
    fChain->SetBranchAddress("pfo_pz", _stats.pfo_pz, &_branch.b_pfo_pz);
    fChain->SetBranchAddress("pfo_m", _stats.pfo_m, &_branch.b_pfo_m);
    fChain->SetBranchAddress("pfo_type", _stats.pfo_type, &_branch.b_pfo_type);
    fChain->SetBranchAddress("pfo_isoverlay", _stats.pfo_isoverlay, &_branch.b_pfo_isoverlay);
    fChain->SetBranchAddress("pfo_isisr", _stats.pfo_isisr, &_branch.b_pfo_isisr);
    fChain->SetBranchAddress("pfo_vtx", _stats.pfo_vtx, &_branch.b_pfo_vtx);
    fChain->SetBranchAddress("pfo_charge", _stats.pfo_charge, &_branch.b_pfo_charge);
    fChain->SetBranchAddress("pfo_ntracks", _stats.pfo_ntracks, &_branch.b_pfo_ntracks);
    fChain->SetBranchAddress("pfo_tpc_hits", _stats.pfo_tpc_hits, &_branch.b_pfo_tpc_hits);
    fChain->SetBranchAddress("pfo_dedx", _stats.pfo_dedx, &_branch.b_pfo_dedx);
    fChain->SetBranchAddress("pfo_dedxerror", _stats.pfo_dedxerror, &_branch.b_pfo_dedxerror);
    fChain->SetBranchAddress("pfo_d0", _stats.pfo_d0, &_branch.b_pfo_d0);
    fChain->SetBranchAddress("pfo_d0error", _stats.pfo_d0error, &_branch.b_pfo_d0error);
    fChain->SetBranchAddress("pfo_z0", _stats.pfo_z0, &_branch.b_pfo_z0);
    fChain->SetBranchAddress("pfo_z0error", _stats.pfo_z0error, &_branch.b_pfo_z0error);
    fChain->SetBranchAddress("pfo_phi", _stats.pfo_phi, &_branch.b_pfo_phi);
    fChain->SetBranchAddress("pfo_phierror", _stats.pfo_phierror, &_branch.b_pfo_phierror);
    fChain->SetBranchAddress("pfo_omega", _stats.pfo_omega, &_branch.b_pfo_omega);
    fChain->SetBranchAddress("pfo_omegaerror", _stats.pfo_omegaerror, &_branch.b_pfo_omegaerror);
    fChain->SetBranchAddress("pfo_tanlambda", _stats.pfo_tanlambda, &_branch.b_pfo_tanlambda);
    fChain->SetBranchAddress("pfo_tanlambdaerror", _stats.pfo_tanlambdaerror, &_branch.b_pfo_tanlambdaerror);
    fChain->SetBranchAddress("pfo_chi2", _stats.pfo_chi2, &_branch.b_pfo_chi2);
    fChain->SetBranchAddress("pfo_ndf", _stats.pfo_ndf, &_branch.b_pfo_ndf);
    fChain->SetBranchAddress("pfo_vtxpt", _stats.pfo_vtxpt, &_branch.b_pfo_vtxpt);
    fChain->SetBranchAddress("pfo_endpt", _stats.pfo_endpt, &_branch.b_pfo_endpt);
    fChain->SetBranchAddress("pfo_pid", _stats.pfo_pid, &_branch.b_pfo_pid);
    fChain->SetBranchAddress("pfo_pid_likelihood", _stats.pfo_pid_likelihood, &_branch.b_pfo_pid_likelihood);
    fChain->SetBranchAddress("pfo_pid_eprob", _stats.pfo_pid_eprob, &_branch.b_pfo_pid_eprob);
    fChain->SetBranchAddress("pfo_pid_muprob", _stats.pfo_pid_muprob, &_branch.b_pfo_pid_muprob);
    fChain->SetBranchAddress("pfo_pid_piprob", _stats.pfo_pid_piprob, &_branch.b_pfo_pid_piprob);
    fChain->SetBranchAddress("pfo_pid_kprob", _stats.pfo_pid_kprob, &_branch.b_pfo_pid_kprob);
    fChain->SetBranchAddress("pfo_pid_pprob", _stats.pfo_pid_pprob, &_branch.b_pfo_pid_pprob);
    fChain->SetBranchAddress("pfo_pid_hprob", _stats.pfo_pid_hprob, &_branch.b_pfo_pid_hprob);
    fChain->SetBranchAddress("pfo_piddedx", _stats.pfo_piddedx, &_branch.b_pfo_piddedx);
    fChain->SetBranchAddress("pfo_piddedx_likelihood", _stats.pfo_piddedx_likelihood, &_branch.b_pfo_piddedx_likelihood);
    fChain->SetBranchAddress("pfo_piddedx_eprob", _stats.pfo_piddedx_eprob, &_branch.b_pfo_piddedx_eprob);
    fChain->SetBranchAddress("pfo_piddedx_muprob", _stats.pfo_piddedx_muprob, &_branch.b_pfo_piddedx_muprob);
    fChain->SetBranchAddress("pfo_piddedx_piprob", _stats.pfo_piddedx_piprob, &_branch.b_pfo_piddedx_piprob);
    fChain->SetBranchAddress("pfo_piddedx_kprob", _stats.pfo_piddedx_kprob, &_branch.b_pfo_piddedx_kprob);
    fChain->SetBranchAddress("pfo_piddedx_pprob", _stats.pfo_piddedx_pprob, &_branch.b_pfo_piddedx_pprob);
    fChain->SetBranchAddress("pfo_piddedx_hprob", _stats.pfo_piddedx_hprob, &_branch.b_pfo_piddedx_hprob);
    fChain->SetBranchAddress("pfo_piddedx_e_dedxdist", _stats.pfo_piddedx_e_dedxdist, &_branch.b_pfo_piddedx_e_dedxdist);
    fChain->SetBranchAddress("pfo_piddedx_mu_dedxdist", _stats.pfo_piddedx_mu_dedxdist, &_branch.b_pfo_piddedx_mu_dedxdist);
    fChain->SetBranchAddress("pfo_piddedx_pi_dedxdist", _stats.pfo_piddedx_pi_dedxdist, &_branch.b_pfo_piddedx_pi_dedxdist);
    fChain->SetBranchAddress("pfo_piddedx_k_dedxdist", _stats.pfo_piddedx_k_dedxdist, &_branch.b_pfo_piddedx_k_dedxdist);
    fChain->SetBranchAddress("pfo_piddedx_p_dedxdist", _stats.pfo_piddedx_p_dedxdist, &_branch.b_pfo_piddedx_p_dedxdist);
    fChain->SetBranchAddress("pfo_piddedx_e_lkhood", _stats.pfo_piddedx_e_lkhood, &_branch.b_pfo_piddedx_e_lkhood);
    fChain->SetBranchAddress("pfo_piddedx_mu_lkhood", _stats.pfo_piddedx_mu_lkhood, &_branch.b_pfo_piddedx_mu_lkhood);
    fChain->SetBranchAddress("pfo_piddedx_pi_lkhood", _stats.pfo_piddedx_pi_lkhood, &_branch.b_pfo_piddedx_pi_lkhood);
    fChain->SetBranchAddress("pfo_piddedx_k_lkhood", _stats.pfo_piddedx_k_lkhood, &_branch.b_pfo_piddedx_k_lkhood);
    fChain->SetBranchAddress("pfo_piddedx_p_lkhood", _stats.pfo_piddedx_p_lkhood, &_branch.b_pfo_piddedx_p_lkhood);
    fChain->SetBranchAddress("pfo_pidtof_p_at_calo", _stats.pfo_pidtof_p_at_calo, &_branch.b_pfo_pidtof_p_at_calo);
    fChain->SetBranchAddress("pfo_pidtof_closest_beta_0ps", _stats.pfo_pidtof_closest_beta_0ps, &_branch.b_pfo_pidtof_closest_beta_0ps);
    fChain->SetBranchAddress("pfo_pidtof_closest_beta_10ps", _stats.pfo_pidtof_closest_beta_10ps, &_branch.b_pfo_pidtof_closest_beta_10ps);
    fChain->SetBranchAddress("pfo_pidtof_closest_beta_50ps", _stats.pfo_pidtof_closest_beta_50ps, &_branch.b_pfo_pidtof_closest_beta_50ps);
    fChain->SetBranchAddress("pfo_pidtof_closest_beta_100ps", _stats.pfo_pidtof_closest_beta_100ps, &_branch.b_pfo_pidtof_closest_beta_100ps);
    fChain->SetBranchAddress("pfo_pidtof_fastest_beta_0ps", _stats.pfo_pidtof_fastest_beta_0ps, &_branch.b_pfo_pidtof_fastest_beta_0ps);
    fChain->SetBranchAddress("pfo_pidtof_fastest_beta_10ps", _stats.pfo_pidtof_fastest_beta_10ps, &_branch.b_pfo_pidtof_fastest_beta_10ps);
    fChain->SetBranchAddress("pfo_pidtof_fastest_beta_50ps", _stats.pfo_pidtof_fastest_beta_50ps, &_branch.b_pfo_pidtof_fastest_beta_50ps);
    fChain->SetBranchAddress("pfo_pidtof_fastest_beta_100ps", _stats.pfo_pidtof_fastest_beta_100ps, &_branch.b_pfo_pidtof_fastest_beta_100ps);
    fChain->SetBranchAddress("pfo_pidtof_cylfit_beta_0ps", _stats.pfo_pidtof_cylfit_beta_0ps, &_branch.b_pfo_pidtof_cylfit_beta_0ps);
    fChain->SetBranchAddress("pfo_pidtof_cylfit_beta_10ps", _stats.pfo_pidtof_cylfit_beta_10ps, &_branch.b_pfo_pidtof_cylfit_beta_10ps);
    fChain->SetBranchAddress("pfo_pidtof_cylfit_beta_50ps", _stats.pfo_pidtof_cylfit_beta_50ps, &_branch.b_pfo_pidtof_cylfit_beta_50ps);
    fChain->SetBranchAddress("pfo_pidtof_cylfit_beta_100ps", _stats.pfo_pidtof_cylfit_beta_100ps, &_branch.b_pfo_pidtof_cylfit_beta_100ps);
    fChain->SetBranchAddress("pfo_pidtof_closestfit_beta_0ps", _stats.pfo_pidtof_closestfit_beta_0ps, &_branch.b_pfo_pidtof_closestfit_beta_0ps);
    fChain->SetBranchAddress("pfo_pidtof_closestfit_beta_10ps", _stats.pfo_pidtof_closestfit_beta_10ps, &_branch.b_pfo_pidtof_closestfit_beta_10ps);
    fChain->SetBranchAddress("pfo_pidtof_closestfit_beta_50ps", _stats.pfo_pidtof_closestfit_beta_50ps, &_branch.b_pfo_pidtof_closestfit_beta_50ps);
    fChain->SetBranchAddress("pfo_pidtof_closestfit_beta_100ps", _stats.pfo_pidtof_closestfit_beta_100ps, &_branch.b_pfo_pidtof_closestfit_beta_100ps);

    Notify();

    return true;

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

void EventAnalyzer::evalCriteria()
{ // Evaluates the class' list of event selection criteria

  cout << "mc_quark_E = " << _stats.mc_quark_E[0] << endl;

}