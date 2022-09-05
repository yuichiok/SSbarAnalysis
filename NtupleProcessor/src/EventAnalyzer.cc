
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

    fChain->SetBranchAddress("mc_quark_E", mc_quark_E, &b_mc_quark_E);
    fChain->SetBranchAddress("mc_quark_px", mc_quark_px, &b_mc_quark_px);
    fChain->SetBranchAddress("mc_quark_py", mc_quark_py, &b_mc_quark_py);
    fChain->SetBranchAddress("mc_quark_pz", mc_quark_pz, &b_mc_quark_pz);
    fChain->SetBranchAddress("mc_quark_m", mc_quark_m, &b_mc_quark_m);
    fChain->SetBranchAddress("mc_quark_pdg", mc_quark_pdg, &b_mc_quark_pdg);
    fChain->SetBranchAddress("mc_quark_charge", mc_quark_charge, &b_mc_quark_charge);
    fChain->SetBranchAddress("mc_ISR_E", mc_ISR_E, &b_mc_ISR_E);
    fChain->SetBranchAddress("mc_ISR_px", mc_ISR_px, &b_mc_ISR_px);
    fChain->SetBranchAddress("mc_ISR_py", mc_ISR_py, &b_mc_ISR_py);
    fChain->SetBranchAddress("mc_ISR_pz", mc_ISR_pz, &b_mc_ISR_pz);
    fChain->SetBranchAddress("mc_ISR_m", mc_ISR_m, &b_mc_ISR_m);
    fChain->SetBranchAddress("mc_ISR_pdg", mc_ISR_pdg, &b_mc_ISR_pdg);
    fChain->SetBranchAddress("mc_ISR_charge", mc_ISR_charge, &b_mc_ISR_charge);
    fChain->SetBranchAddress("mc_quark_ps_n", &mc_quark_ps_n, &b_mc_quark_ps_n);
    fChain->SetBranchAddress("mc_quark_ps_E", mc_quark_ps_E, &b_mc_quark_ps_E);
    fChain->SetBranchAddress("mc_quark_ps_px", mc_quark_ps_px, &b_mc_quark_ps_px);
    fChain->SetBranchAddress("mc_quark_ps_py", mc_quark_ps_py, &b_mc_quark_ps_py);
    fChain->SetBranchAddress("mc_quark_ps_pz", mc_quark_ps_pz, &b_mc_quark_ps_pz);
    fChain->SetBranchAddress("mc_quark_ps_m", mc_quark_ps_m, &b_mc_quark_ps_m);
    fChain->SetBranchAddress("mc_quark_ps_pdg", mc_quark_ps_pdg, &b_mc_quark_ps_pdg);
    fChain->SetBranchAddress("mc_quark_ps_charge", mc_quark_ps_charge, &b_mc_quark_ps_charge);
    fChain->SetBranchAddress("mc_quark_ps_y12", &mc_quark_ps_y12, &b_mc_quark_ps_y12);
    fChain->SetBranchAddress("mc_quark_ps_y23", &mc_quark_ps_y23, &b_mc_quark_ps_y23);
    fChain->SetBranchAddress("mc_quark_ps_d12", &mc_quark_ps_d12, &b_mc_quark_ps_d12);
    fChain->SetBranchAddress("mc_quark_ps_d23", &mc_quark_ps_d23, &b_mc_quark_ps_d23);
    fChain->SetBranchAddress("mc_quark_ps_jet_E", mc_quark_ps_jet_E, &b_mc_quark_ps_jet_E);
    fChain->SetBranchAddress("mc_quark_ps_jet_px", mc_quark_ps_jet_px, &b_mc_quark_ps_jet_px);
    fChain->SetBranchAddress("mc_quark_ps_jet_py", mc_quark_ps_jet_py, &b_mc_quark_ps_jet_py);
    fChain->SetBranchAddress("mc_quark_ps_jet_pz", mc_quark_ps_jet_pz, &b_mc_quark_ps_jet_pz);
    fChain->SetBranchAddress("mc_stable_n", &mc_stable_n, &b_mc_stable_n);
    fChain->SetBranchAddress("mc_stable_E", mc_stable_E, &b_mc_stable_E);
    fChain->SetBranchAddress("mc_stable_px", mc_stable_px, &b_mc_stable_px);
    fChain->SetBranchAddress("mc_stable_py", mc_stable_py, &b_mc_stable_py);
    fChain->SetBranchAddress("mc_stable_pz", mc_stable_pz, &b_mc_stable_pz);
    fChain->SetBranchAddress("mc_stable_m", mc_stable_m, &b_mc_stable_m);
    fChain->SetBranchAddress("mc_stable_pdg", mc_stable_pdg, &b_mc_stable_pdg);
    fChain->SetBranchAddress("mc_stable_charge", mc_stable_charge, &b_mc_stable_charge);
    fChain->SetBranchAddress("mc_stable_isoverlay", mc_stable_isoverlay, &b_mc_stable_isoverlay);
    fChain->SetBranchAddress("mc_stable_isisr", mc_stable_isisr, &b_mc_stable_isisr);
    fChain->SetBranchAddress("mc_stable_y12", &mc_stable_y12, &b_mc_stable_y12);
    fChain->SetBranchAddress("mc_stable_y23", &mc_stable_y23, &b_mc_stable_y23);
    fChain->SetBranchAddress("mc_stable_d12", &mc_stable_d12, &b_mc_stable_d12);
    fChain->SetBranchAddress("mc_stable_d23", &mc_stable_d23, &b_mc_stable_d23);
    fChain->SetBranchAddress("mc_stable_jet_E", mc_stable_jet_E, &b_mc_stable_jet_E);
    fChain->SetBranchAddress("mc_stable_jet_px", mc_stable_jet_px, &b_mc_stable_jet_px);
    fChain->SetBranchAddress("mc_stable_jet_py", mc_stable_jet_py, &b_mc_stable_jet_py);
    fChain->SetBranchAddress("mc_stable_jet_pz", mc_stable_jet_pz, &b_mc_stable_jet_pz);
    fChain->SetBranchAddress("truejet_E", truejet_E, &b_truejet_E);
    fChain->SetBranchAddress("truejet_px", truejet_px, &b_truejet_px);
    fChain->SetBranchAddress("truejet_py", truejet_py, &b_truejet_py);
    fChain->SetBranchAddress("truejet_pz", truejet_pz, &b_truejet_pz);
    fChain->SetBranchAddress("truejet_type", truejet_type, &b_truejet_type);
    fChain->SetBranchAddress("truejet_pdg", truejet_pdg, &b_truejet_pdg);
    fChain->SetBranchAddress("jet_E", jet_E, &b_jet_E);
    fChain->SetBranchAddress("jet_px", jet_px, &b_jet_px);
    fChain->SetBranchAddress("jet_py", jet_py, &b_jet_py);
    fChain->SetBranchAddress("jet_pz", jet_pz, &b_jet_pz);
    fChain->SetBranchAddress("jet_btag", jet_btag, &b_jet_btag);
    fChain->SetBranchAddress("jet_ctag", jet_ctag, &b_jet_ctag);
    fChain->SetBranchAddress("y23", &y23, &b_y23);
    fChain->SetBranchAddress("y12", &y12, &b_y12);
    fChain->SetBranchAddress("d23", &d23, &b_d23);
    fChain->SetBranchAddress("d12", &d12, &b_d12);
    fChain->SetBranchAddress("oblateness", &oblateness, &b_oblateness);
    fChain->SetBranchAddress("aplanarity", &aplanarity, &b_aplanarity);
    fChain->SetBranchAddress("major_thrust_value", &major_thrust_value, &b_major_thrust_value);
    fChain->SetBranchAddress("major_thrust_axis", major_thrust_axis, &b_major_thrust_axis);
    fChain->SetBranchAddress("minor_thrust_value", &minor_thrust_value, &b_minor_thrust_value);
    fChain->SetBranchAddress("minor_thrust_axis", minor_thrust_axis, &b_minor_thrust_axis);
    fChain->SetBranchAddress("principle_thrust_value", &principle_thrust_value, &b_principle_thrust_value);
    fChain->SetBranchAddress("principle_thrust_axis", principle_thrust_axis, &b_principle_thrust_axis);
    fChain->SetBranchAddress("sphericity", &sphericity, &b_sphericity);
    fChain->SetBranchAddress("sphericity_tensor", sphericity_tensor, &b_sphericity_tensor);
    fChain->SetBranchAddress("pfo_n", &pfo_n, &b_pfo_n);
    fChain->SetBranchAddress("jet_nvtx", &jet_nvtx, &b_jet_nvtx);
    fChain->SetBranchAddress("pfo_n_j1", &pfo_n_j1, &b_pfo_n_j1);
    fChain->SetBranchAddress("jet_nvtx_j1", &jet_nvtx_j1, &b_jet_nvtx_j1);
    fChain->SetBranchAddress("pfo_n_j2", &pfo_n_j2, &b_pfo_n_j2);
    fChain->SetBranchAddress("jet_nvtx_j2", &jet_nvtx_j2, &b_jet_nvtx_j2);
    fChain->SetBranchAddress("pfo_match", pfo_match, &b_pfo_match);
    fChain->SetBranchAddress("pfo_truejet_pdg", pfo_truejet_pdg, &b_pfo_truejet_pdg);
    fChain->SetBranchAddress("pfo_truejet_type", pfo_truejet_type, &b_pfo_truejet_type);
    fChain->SetBranchAddress("pfo_pdgcheat", pfo_pdgcheat, &b_pfo_pdgcheat);
    fChain->SetBranchAddress("pfo_nparents", pfo_nparents, &b_pfo_nparents);
    fChain->SetBranchAddress("pfo_pdgcheat_parent", pfo_pdgcheat_parent, &b_pfo_pdgcheat_parent);
    fChain->SetBranchAddress("pfo_E", pfo_E, &b_pfo_E);
    fChain->SetBranchAddress("pfo_E_calo", pfo_E_calo, &b_pfo_E_calo);
    fChain->SetBranchAddress("pfo_px", pfo_px, &b_pfo_px);
    fChain->SetBranchAddress("pfo_py", pfo_py, &b_pfo_py);
    fChain->SetBranchAddress("pfo_pz", pfo_pz, &b_pfo_pz);
    fChain->SetBranchAddress("pfo_m", pfo_m, &b_pfo_m);
    fChain->SetBranchAddress("pfo_type", pfo_type, &b_pfo_type);
    fChain->SetBranchAddress("pfo_isoverlay", pfo_isoverlay, &b_pfo_isoverlay);
    fChain->SetBranchAddress("pfo_isisr", pfo_isisr, &b_pfo_isisr);
    fChain->SetBranchAddress("pfo_vtx", pfo_vtx, &b_pfo_vtx);
    fChain->SetBranchAddress("pfo_charge", pfo_charge, &b_pfo_charge);
    fChain->SetBranchAddress("pfo_ntracks", pfo_ntracks, &b_pfo_ntracks);
    fChain->SetBranchAddress("pfo_tpc_hits", pfo_tpc_hits, &b_pfo_tpc_hits);
    fChain->SetBranchAddress("pfo_dedx", pfo_dedx, &b_pfo_dedx);
    fChain->SetBranchAddress("pfo_dedxerror", pfo_dedxerror, &b_pfo_dedxerror);
    fChain->SetBranchAddress("pfo_d0", pfo_d0, &b_pfo_d0);
    fChain->SetBranchAddress("pfo_d0error", pfo_d0error, &b_pfo_d0error);
    fChain->SetBranchAddress("pfo_z0", pfo_z0, &b_pfo_z0);
    fChain->SetBranchAddress("pfo_z0error", pfo_z0error, &b_pfo_z0error);
    fChain->SetBranchAddress("pfo_phi", pfo_phi, &b_pfo_phi);
    fChain->SetBranchAddress("pfo_phierror", pfo_phierror, &b_pfo_phierror);
    fChain->SetBranchAddress("pfo_omega", pfo_omega, &b_pfo_omega);
    fChain->SetBranchAddress("pfo_omegaerror", pfo_omegaerror, &b_pfo_omegaerror);
    fChain->SetBranchAddress("pfo_tanlambda", pfo_tanlambda, &b_pfo_tanlambda);
    fChain->SetBranchAddress("pfo_tanlambdaerror", pfo_tanlambdaerror, &b_pfo_tanlambdaerror);
    fChain->SetBranchAddress("pfo_chi2", pfo_chi2, &b_pfo_chi2);
    fChain->SetBranchAddress("pfo_ndf", pfo_ndf, &b_pfo_ndf);
    fChain->SetBranchAddress("pfo_vtxpt", pfo_vtxpt, &b_pfo_vtxpt);
    fChain->SetBranchAddress("pfo_endpt", pfo_endpt, &b_pfo_endpt);
    fChain->SetBranchAddress("pfo_pid", pfo_pid, &b_pfo_pid);
    fChain->SetBranchAddress("pfo_pid_likelihood", pfo_pid_likelihood, &b_pfo_pid_likelihood);
    fChain->SetBranchAddress("pfo_pid_eprob", pfo_pid_eprob, &b_pfo_pid_eprob);
    fChain->SetBranchAddress("pfo_pid_muprob", pfo_pid_muprob, &b_pfo_pid_muprob);
    fChain->SetBranchAddress("pfo_pid_piprob", pfo_pid_piprob, &b_pfo_pid_piprob);
    fChain->SetBranchAddress("pfo_pid_kprob", pfo_pid_kprob, &b_pfo_pid_kprob);
    fChain->SetBranchAddress("pfo_pid_pprob", pfo_pid_pprob, &b_pfo_pid_pprob);
    fChain->SetBranchAddress("pfo_pid_hprob", pfo_pid_hprob, &b_pfo_pid_hprob);
    fChain->SetBranchAddress("pfo_piddedx", pfo_piddedx, &b_pfo_piddedx);
    fChain->SetBranchAddress("pfo_piddedx_likelihood", pfo_piddedx_likelihood, &b_pfo_piddedx_likelihood);
    fChain->SetBranchAddress("pfo_piddedx_eprob", pfo_piddedx_eprob, &b_pfo_piddedx_eprob);
    fChain->SetBranchAddress("pfo_piddedx_muprob", pfo_piddedx_muprob, &b_pfo_piddedx_muprob);
    fChain->SetBranchAddress("pfo_piddedx_piprob", pfo_piddedx_piprob, &b_pfo_piddedx_piprob);
    fChain->SetBranchAddress("pfo_piddedx_kprob", pfo_piddedx_kprob, &b_pfo_piddedx_kprob);
    fChain->SetBranchAddress("pfo_piddedx_pprob", pfo_piddedx_pprob, &b_pfo_piddedx_pprob);
    fChain->SetBranchAddress("pfo_piddedx_hprob", pfo_piddedx_hprob, &b_pfo_piddedx_hprob);
    fChain->SetBranchAddress("pfo_piddedx_e_dedxdist", pfo_piddedx_e_dedxdist, &b_pfo_piddedx_e_dedxdist);
    fChain->SetBranchAddress("pfo_piddedx_mu_dedxdist", pfo_piddedx_mu_dedxdist, &b_pfo_piddedx_mu_dedxdist);
    fChain->SetBranchAddress("pfo_piddedx_pi_dedxdist", pfo_piddedx_pi_dedxdist, &b_pfo_piddedx_pi_dedxdist);
    fChain->SetBranchAddress("pfo_piddedx_k_dedxdist", pfo_piddedx_k_dedxdist, &b_pfo_piddedx_k_dedxdist);
    fChain->SetBranchAddress("pfo_piddedx_p_dedxdist", pfo_piddedx_p_dedxdist, &b_pfo_piddedx_p_dedxdist);
    fChain->SetBranchAddress("pfo_piddedx_e_lkhood", pfo_piddedx_e_lkhood, &b_pfo_piddedx_e_lkhood);
    fChain->SetBranchAddress("pfo_piddedx_mu_lkhood", pfo_piddedx_mu_lkhood, &b_pfo_piddedx_mu_lkhood);
    fChain->SetBranchAddress("pfo_piddedx_pi_lkhood", pfo_piddedx_pi_lkhood, &b_pfo_piddedx_pi_lkhood);
    fChain->SetBranchAddress("pfo_piddedx_k_lkhood", pfo_piddedx_k_lkhood, &b_pfo_piddedx_k_lkhood);
    fChain->SetBranchAddress("pfo_piddedx_p_lkhood", pfo_piddedx_p_lkhood, &b_pfo_piddedx_p_lkhood);
    fChain->SetBranchAddress("pfo_pidtof_p_at_calo", pfo_pidtof_p_at_calo, &b_pfo_pidtof_p_at_calo);
    fChain->SetBranchAddress("pfo_pidtof_closest_beta_0ps", pfo_pidtof_closest_beta_0ps, &b_pfo_pidtof_closest_beta_0ps);
    fChain->SetBranchAddress("pfo_pidtof_closest_beta_10ps", pfo_pidtof_closest_beta_10ps, &b_pfo_pidtof_closest_beta_10ps);
    fChain->SetBranchAddress("pfo_pidtof_closest_beta_50ps", pfo_pidtof_closest_beta_50ps, &b_pfo_pidtof_closest_beta_50ps);
    fChain->SetBranchAddress("pfo_pidtof_closest_beta_100ps", pfo_pidtof_closest_beta_100ps, &b_pfo_pidtof_closest_beta_100ps);
    fChain->SetBranchAddress("pfo_pidtof_fastest_beta_0ps", pfo_pidtof_fastest_beta_0ps, &b_pfo_pidtof_fastest_beta_0ps);
    fChain->SetBranchAddress("pfo_pidtof_fastest_beta_10ps", pfo_pidtof_fastest_beta_10ps, &b_pfo_pidtof_fastest_beta_10ps);
    fChain->SetBranchAddress("pfo_pidtof_fastest_beta_50ps", pfo_pidtof_fastest_beta_50ps, &b_pfo_pidtof_fastest_beta_50ps);
    fChain->SetBranchAddress("pfo_pidtof_fastest_beta_100ps", pfo_pidtof_fastest_beta_100ps, &b_pfo_pidtof_fastest_beta_100ps);
    fChain->SetBranchAddress("pfo_pidtof_cylfit_beta_0ps", pfo_pidtof_cylfit_beta_0ps, &b_pfo_pidtof_cylfit_beta_0ps);
    fChain->SetBranchAddress("pfo_pidtof_cylfit_beta_10ps", pfo_pidtof_cylfit_beta_10ps, &b_pfo_pidtof_cylfit_beta_10ps);
    fChain->SetBranchAddress("pfo_pidtof_cylfit_beta_50ps", pfo_pidtof_cylfit_beta_50ps, &b_pfo_pidtof_cylfit_beta_50ps);
    fChain->SetBranchAddress("pfo_pidtof_cylfit_beta_100ps", pfo_pidtof_cylfit_beta_100ps, &b_pfo_pidtof_cylfit_beta_100ps);
    fChain->SetBranchAddress("pfo_pidtof_closestfit_beta_0ps", pfo_pidtof_closestfit_beta_0ps, &b_pfo_pidtof_closestfit_beta_0ps);
    fChain->SetBranchAddress("pfo_pidtof_closestfit_beta_10ps", pfo_pidtof_closestfit_beta_10ps, &b_pfo_pidtof_closestfit_beta_10ps);
    fChain->SetBranchAddress("pfo_pidtof_closestfit_beta_50ps", pfo_pidtof_closestfit_beta_50ps, &b_pfo_pidtof_closestfit_beta_50ps);
    fChain->SetBranchAddress("pfo_pidtof_closestfit_beta_100ps", pfo_pidtof_closestfit_beta_100ps, &b_pfo_pidtof_closestfit_beta_100ps);

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
}