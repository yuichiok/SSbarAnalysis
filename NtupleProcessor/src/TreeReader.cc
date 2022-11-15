
/*------------------------------------------------------------------------------
TreeReader.cpp
 Created : 2022-09-08  okugawa
------------------------------------------------------------------------------*/

#include <iostream>
#include <TString.h>
#include <TFile.h> 
#include "TreeReader.hh"

using std::cout;   using std::endl;

TreeReader::TreeReader() {}

void TreeReader::InitializeMCReadTree(TTree *_hTree, MC_QQbar & _data, Branch_QQbar & _branch)
{
    _hTree->SetBranchAddress("mc_quark_E", _data.mc_quark_E, &_branch.b_mc_quark_E);
    _hTree->SetBranchAddress("mc_quark_px", _data.mc_quark_px, &_branch.b_mc_quark_px);
    _hTree->SetBranchAddress("mc_quark_py", _data.mc_quark_py, &_branch.b_mc_quark_py);
    _hTree->SetBranchAddress("mc_quark_pz", _data.mc_quark_pz, &_branch.b_mc_quark_pz);
    _hTree->SetBranchAddress("mc_quark_m", _data.mc_quark_m, &_branch.b_mc_quark_m);
    _hTree->SetBranchAddress("mc_quark_pdg", _data.mc_quark_pdg, &_branch.b_mc_quark_pdg);
    _hTree->SetBranchAddress("mc_quark_charge", _data.mc_quark_charge, &_branch.b_mc_quark_charge);
    _hTree->SetBranchAddress("mc_ISR_E", _data.mc_ISR_E, &_branch.b_mc_ISR_E);
    _hTree->SetBranchAddress("mc_ISR_px", _data.mc_ISR_px, &_branch.b_mc_ISR_px);
    _hTree->SetBranchAddress("mc_ISR_py", _data.mc_ISR_py, &_branch.b_mc_ISR_py);
    _hTree->SetBranchAddress("mc_ISR_pz", _data.mc_ISR_pz, &_branch.b_mc_ISR_pz);
    _hTree->SetBranchAddress("mc_ISR_m", _data.mc_ISR_m, &_branch.b_mc_ISR_m);
    _hTree->SetBranchAddress("mc_ISR_pdg", _data.mc_ISR_pdg, &_branch.b_mc_ISR_pdg);
    _hTree->SetBranchAddress("mc_ISR_charge", _data.mc_ISR_charge, &_branch.b_mc_ISR_charge);
    _hTree->SetBranchAddress("mc_quark_ps_n", &_data.mc_quark_ps_n, &_branch.b_mc_quark_ps_n);
    _hTree->SetBranchAddress("mc_quark_ps_E", _data.mc_quark_ps_E, &_branch.b_mc_quark_ps_E);
    _hTree->SetBranchAddress("mc_quark_ps_px", _data.mc_quark_ps_px, &_branch.b_mc_quark_ps_px);
    _hTree->SetBranchAddress("mc_quark_ps_py", _data.mc_quark_ps_py, &_branch.b_mc_quark_ps_py);
    _hTree->SetBranchAddress("mc_quark_ps_pz", _data.mc_quark_ps_pz, &_branch.b_mc_quark_ps_pz);
    _hTree->SetBranchAddress("mc_quark_ps_m", _data.mc_quark_ps_m, &_branch.b_mc_quark_ps_m);
    _hTree->SetBranchAddress("mc_quark_ps_pdg", _data.mc_quark_ps_pdg, &_branch.b_mc_quark_ps_pdg);
    _hTree->SetBranchAddress("mc_quark_ps_charge", _data.mc_quark_ps_charge, &_branch.b_mc_quark_ps_charge);
    _hTree->SetBranchAddress("mc_quark_ps_y12", &_data.mc_quark_ps_y12, &_branch.b_mc_quark_ps_y12);
    _hTree->SetBranchAddress("mc_quark_ps_y23", &_data.mc_quark_ps_y23, &_branch.b_mc_quark_ps_y23);
    _hTree->SetBranchAddress("mc_quark_ps_d12", &_data.mc_quark_ps_d12, &_branch.b_mc_quark_ps_d12);
    _hTree->SetBranchAddress("mc_quark_ps_d23", &_data.mc_quark_ps_d23, &_branch.b_mc_quark_ps_d23);
    _hTree->SetBranchAddress("mc_quark_ps_jet_E", _data.mc_quark_ps_jet_E, &_branch.b_mc_quark_ps_jet_E);
    _hTree->SetBranchAddress("mc_quark_ps_jet_px", _data.mc_quark_ps_jet_px, &_branch.b_mc_quark_ps_jet_px);
    _hTree->SetBranchAddress("mc_quark_ps_jet_py", _data.mc_quark_ps_jet_py, &_branch.b_mc_quark_ps_jet_py);
    _hTree->SetBranchAddress("mc_quark_ps_jet_pz", _data.mc_quark_ps_jet_pz, &_branch.b_mc_quark_ps_jet_pz);
    _hTree->SetBranchAddress("mc_stable_n", &_data.mc_stable_n, &_branch.b_mc_stable_n);
    _hTree->SetBranchAddress("mc_stable_E", _data.mc_stable_E, &_branch.b_mc_stable_E);
    _hTree->SetBranchAddress("mc_stable_px", _data.mc_stable_px, &_branch.b_mc_stable_px);
    _hTree->SetBranchAddress("mc_stable_py", _data.mc_stable_py, &_branch.b_mc_stable_py);
    _hTree->SetBranchAddress("mc_stable_pz", _data.mc_stable_pz, &_branch.b_mc_stable_pz);
    _hTree->SetBranchAddress("mc_stable_m", _data.mc_stable_m, &_branch.b_mc_stable_m);
    _hTree->SetBranchAddress("mc_stable_pdg", _data.mc_stable_pdg, &_branch.b_mc_stable_pdg);
    _hTree->SetBranchAddress("mc_stable_charge", _data.mc_stable_charge, &_branch.b_mc_stable_charge);
    _hTree->SetBranchAddress("mc_stable_isoverlay", _data.mc_stable_isoverlay, &_branch.b_mc_stable_isoverlay);
    _hTree->SetBranchAddress("mc_stable_isisr", _data.mc_stable_isisr, &_branch.b_mc_stable_isisr);
    _hTree->SetBranchAddress("mc_stable_y12", &_data.mc_stable_y12, &_branch.b_mc_stable_y12);
    _hTree->SetBranchAddress("mc_stable_y23", &_data.mc_stable_y23, &_branch.b_mc_stable_y23);
    _hTree->SetBranchAddress("mc_stable_d12", &_data.mc_stable_d12, &_branch.b_mc_stable_d12);
    _hTree->SetBranchAddress("mc_stable_d23", &_data.mc_stable_d23, &_branch.b_mc_stable_d23);
    _hTree->SetBranchAddress("mc_stable_jet_E", _data.mc_stable_jet_E, &_branch.b_mc_stable_jet_E);
    _hTree->SetBranchAddress("mc_stable_jet_px", _data.mc_stable_jet_px, &_branch.b_mc_stable_jet_px);
    _hTree->SetBranchAddress("mc_stable_jet_py", _data.mc_stable_jet_py, &_branch.b_mc_stable_jet_py);
    _hTree->SetBranchAddress("mc_stable_jet_pz", _data.mc_stable_jet_pz, &_branch.b_mc_stable_jet_pz);
}

void TreeReader::InitializeJetReadTree(TTree *_hTree, Jet_QQbar & _data, Branch_QQbar & _branch)
{
    _hTree->SetBranchAddress("truejet_E", _data.truejet_E, &_branch.b_truejet_E);
    _hTree->SetBranchAddress("truejet_px", _data.truejet_px, &_branch.b_truejet_px);
    _hTree->SetBranchAddress("truejet_py", _data.truejet_py, &_branch.b_truejet_py);
    _hTree->SetBranchAddress("truejet_pz", _data.truejet_pz, &_branch.b_truejet_pz);
    _hTree->SetBranchAddress("truejet_type", _data.truejet_type, &_branch.b_truejet_type);
    _hTree->SetBranchAddress("truejet_pdg", _data.truejet_pdg, &_branch.b_truejet_pdg);
    _hTree->SetBranchAddress("jet_E", _data.jet_E, &_branch.b_jet_E);
    _hTree->SetBranchAddress("jet_px", _data.jet_px, &_branch.b_jet_px);
    _hTree->SetBranchAddress("jet_py", _data.jet_py, &_branch.b_jet_py);
    _hTree->SetBranchAddress("jet_pz", _data.jet_pz, &_branch.b_jet_pz);
    _hTree->SetBranchAddress("jet_btag", _data.jet_btag, &_branch.b_jet_btag);
    _hTree->SetBranchAddress("jet_ctag", _data.jet_ctag, &_branch.b_jet_ctag);
    _hTree->SetBranchAddress("y23", &_data.y23, &_branch.b_y23);
    _hTree->SetBranchAddress("y12", &_data.y12, &_branch.b_y12);
    _hTree->SetBranchAddress("d23", &_data.d23, &_branch.b_d23);
    _hTree->SetBranchAddress("d12", &_data.d12, &_branch.b_d12);
    _hTree->SetBranchAddress("oblateness", &_data.oblateness, &_branch.b_oblateness);
    _hTree->SetBranchAddress("aplanarity", &_data.aplanarity, &_branch.b_aplanarity);
    _hTree->SetBranchAddress("major_thrust_value", &_data.major_thrust_value, &_branch.b_major_thrust_value);
    _hTree->SetBranchAddress("major_thrust_axis", _data.major_thrust_axis, &_branch.b_major_thrust_axis);
    _hTree->SetBranchAddress("minor_thrust_value", &_data.minor_thrust_value, &_branch.b_minor_thrust_value);
    _hTree->SetBranchAddress("minor_thrust_axis", _data.minor_thrust_axis, &_branch.b_minor_thrust_axis);
    _hTree->SetBranchAddress("principle_thrust_value", &_data.principle_thrust_value, &_branch.b_principle_thrust_value);
    _hTree->SetBranchAddress("principle_thrust_axis", _data.principle_thrust_axis, &_branch.b_principle_thrust_axis);
    _hTree->SetBranchAddress("sphericity", &_data.sphericity, &_branch.b_sphericity);
    _hTree->SetBranchAddress("sphericity_tensor", _data.sphericity_tensor, &_branch.b_sphericity_tensor);
}

void TreeReader::InitializePFOReadTree(TTree *_hTree, PFO_QQbar & _data, Branch_QQbar & _branch)
{
    _hTree->SetBranchAddress("pfo_n", &_data.pfo_n, &_branch.b_pfo_n);
    _hTree->SetBranchAddress("jet_nvtx", &_data.jet_nvtx, &_branch.b_jet_nvtx);
    _hTree->SetBranchAddress("pfo_n_j1", &_data.pfo_n_j1, &_branch.b_pfo_n_j1);
    _hTree->SetBranchAddress("jet_nvtx_j1", &_data.jet_nvtx_j1, &_branch.b_jet_nvtx_j1);
    _hTree->SetBranchAddress("pfo_n_j2", &_data.pfo_n_j2, &_branch.b_pfo_n_j2);
    _hTree->SetBranchAddress("jet_nvtx_j2", &_data.jet_nvtx_j2, &_branch.b_jet_nvtx_j2);
    _hTree->SetBranchAddress("pfo_match", _data.pfo_match, &_branch.b_pfo_match);
    _hTree->SetBranchAddress("pfo_truejet_pdg", _data.pfo_truejet_pdg, &_branch.b_pfo_truejet_pdg);
    _hTree->SetBranchAddress("pfo_truejet_type", _data.pfo_truejet_type, &_branch.b_pfo_truejet_type);
    _hTree->SetBranchAddress("pfo_pdgcheat", _data.pfo_pdgcheat, &_branch.b_pfo_pdgcheat);
    _hTree->SetBranchAddress("pfo_pdgcheat_id", _data.pfo_pdgcheat_id, &_branch.b_pfo_pdgcheat_id);
    _hTree->SetBranchAddress("pfo_nparents", _data.pfo_nparents, &_branch.b_pfo_nparents);
    _hTree->SetBranchAddress("pfo_pdgcheat_parent", _data.pfo_pdgcheat_parent, &_branch.b_pfo_pdgcheat_parent);
    _hTree->SetBranchAddress("pfo_E", _data.pfo_E, &_branch.b_pfo_E);
    _hTree->SetBranchAddress("pfo_px", _data.pfo_px, &_branch.b_pfo_px);
    _hTree->SetBranchAddress("pfo_py", _data.pfo_py, &_branch.b_pfo_py);
    _hTree->SetBranchAddress("pfo_pz", _data.pfo_pz, &_branch.b_pfo_pz);
    _hTree->SetBranchAddress("pfo_m", _data.pfo_m, &_branch.b_pfo_m);
    _hTree->SetBranchAddress("pfo_type", _data.pfo_type, &_branch.b_pfo_type);
    _hTree->SetBranchAddress("pfo_isoverlay", _data.pfo_isoverlay, &_branch.b_pfo_isoverlay);
    _hTree->SetBranchAddress("pfo_isisr", _data.pfo_isisr, &_branch.b_pfo_isisr);
    _hTree->SetBranchAddress("pfo_vtx", _data.pfo_vtx, &_branch.b_pfo_vtx);
    _hTree->SetBranchAddress("pfo_charge", _data.pfo_charge, &_branch.b_pfo_charge);
    _hTree->SetBranchAddress("pfo_ntracks", _data.pfo_ntracks, &_branch.b_pfo_ntracks);
    _hTree->SetBranchAddress("pfo_tpc_hits", _data.pfo_tpc_hits, &_branch.b_pfo_tpc_hits);
    _hTree->SetBranchAddress("pfo_dedx", _data.pfo_dedx, &_branch.b_pfo_dedx);
    _hTree->SetBranchAddress("pfo_dedxerror", _data.pfo_dedxerror, &_branch.b_pfo_dedxerror);
    _hTree->SetBranchAddress("pfo_d0", _data.pfo_d0, &_branch.b_pfo_d0);
    _hTree->SetBranchAddress("pfo_d0error", _data.pfo_d0error, &_branch.b_pfo_d0error);
    _hTree->SetBranchAddress("pfo_z0", _data.pfo_z0, &_branch.b_pfo_z0);
    _hTree->SetBranchAddress("pfo_z0error", _data.pfo_z0error, &_branch.b_pfo_z0error);
    _hTree->SetBranchAddress("pfo_phi", _data.pfo_phi, &_branch.b_pfo_phi);
    _hTree->SetBranchAddress("pfo_phierror", _data.pfo_phierror, &_branch.b_pfo_phierror);
    _hTree->SetBranchAddress("pfo_omega", _data.pfo_omega, &_branch.b_pfo_omega);
    _hTree->SetBranchAddress("pfo_omegaerror", _data.pfo_omegaerror, &_branch.b_pfo_omegaerror);
    _hTree->SetBranchAddress("pfo_tanlambda", _data.pfo_tanlambda, &_branch.b_pfo_tanlambda);
    _hTree->SetBranchAddress("pfo_tanlambdaerror", _data.pfo_tanlambdaerror, &_branch.b_pfo_tanlambdaerror);
    _hTree->SetBranchAddress("pfo_chi2", _data.pfo_chi2, &_branch.b_pfo_chi2);
    _hTree->SetBranchAddress("pfo_ndf", _data.pfo_ndf, &_branch.b_pfo_ndf);
    _hTree->SetBranchAddress("pfo_vtxpt", _data.pfo_vtxpt, &_branch.b_pfo_vtxpt);
    _hTree->SetBranchAddress("pfo_endpt", _data.pfo_endpt, &_branch.b_pfo_endpt);
    _hTree->SetBranchAddress("pfo_pid", _data.pfo_pid, &_branch.b_pfo_pid);
    _hTree->SetBranchAddress("pfo_pid_likelihood", _data.pfo_pid_likelihood, &_branch.b_pfo_pid_likelihood);
    _hTree->SetBranchAddress("pfo_pid_eprob", _data.pfo_pid_eprob, &_branch.b_pfo_pid_eprob);
    _hTree->SetBranchAddress("pfo_pid_muprob", _data.pfo_pid_muprob, &_branch.b_pfo_pid_muprob);
    _hTree->SetBranchAddress("pfo_pid_piprob", _data.pfo_pid_piprob, &_branch.b_pfo_pid_piprob);
    _hTree->SetBranchAddress("pfo_pid_kprob", _data.pfo_pid_kprob, &_branch.b_pfo_pid_kprob);
    _hTree->SetBranchAddress("pfo_pid_pprob", _data.pfo_pid_pprob, &_branch.b_pfo_pid_pprob);
    _hTree->SetBranchAddress("pfo_pid_hprob", _data.pfo_pid_hprob, &_branch.b_pfo_pid_hprob);
    _hTree->SetBranchAddress("pfo_piddedx", _data.pfo_piddedx, &_branch.b_pfo_piddedx);
    _hTree->SetBranchAddress("pfo_piddedx_likelihood", _data.pfo_piddedx_likelihood, &_branch.b_pfo_piddedx_likelihood);
    _hTree->SetBranchAddress("pfo_piddedx_eprob", _data.pfo_piddedx_eprob, &_branch.b_pfo_piddedx_eprob);
    _hTree->SetBranchAddress("pfo_piddedx_muprob", _data.pfo_piddedx_muprob, &_branch.b_pfo_piddedx_muprob);
    _hTree->SetBranchAddress("pfo_piddedx_piprob", _data.pfo_piddedx_piprob, &_branch.b_pfo_piddedx_piprob);
    _hTree->SetBranchAddress("pfo_piddedx_kprob", _data.pfo_piddedx_kprob, &_branch.b_pfo_piddedx_kprob);
    _hTree->SetBranchAddress("pfo_piddedx_pprob", _data.pfo_piddedx_pprob, &_branch.b_pfo_piddedx_pprob);
    _hTree->SetBranchAddress("pfo_piddedx_hprob", _data.pfo_piddedx_hprob, &_branch.b_pfo_piddedx_hprob);
    _hTree->SetBranchAddress("pfo_piddedx_e_dedxdist", _data.pfo_piddedx_e_dedxdist, &_branch.b_pfo_piddedx_e_dedxdist);
    _hTree->SetBranchAddress("pfo_piddedx_mu_dedxdist", _data.pfo_piddedx_mu_dedxdist, &_branch.b_pfo_piddedx_mu_dedxdist);
    _hTree->SetBranchAddress("pfo_piddedx_pi_dedxdist", _data.pfo_piddedx_pi_dedxdist, &_branch.b_pfo_piddedx_pi_dedxdist);
    _hTree->SetBranchAddress("pfo_piddedx_k_dedxdist", _data.pfo_piddedx_k_dedxdist, &_branch.b_pfo_piddedx_k_dedxdist);
    _hTree->SetBranchAddress("pfo_piddedx_p_dedxdist", _data.pfo_piddedx_p_dedxdist, &_branch.b_pfo_piddedx_p_dedxdist);
    _hTree->SetBranchAddress("pfo_piddedx_e_lkhood", _data.pfo_piddedx_e_lkhood, &_branch.b_pfo_piddedx_e_lkhood);
    _hTree->SetBranchAddress("pfo_piddedx_mu_lkhood", _data.pfo_piddedx_mu_lkhood, &_branch.b_pfo_piddedx_mu_lkhood);
    _hTree->SetBranchAddress("pfo_piddedx_pi_lkhood", _data.pfo_piddedx_pi_lkhood, &_branch.b_pfo_piddedx_pi_lkhood);
    _hTree->SetBranchAddress("pfo_piddedx_k_lkhood", _data.pfo_piddedx_k_lkhood, &_branch.b_pfo_piddedx_k_lkhood);
    _hTree->SetBranchAddress("pfo_piddedx_p_lkhood", _data.pfo_piddedx_p_lkhood, &_branch.b_pfo_piddedx_p_lkhood);
    _hTree->SetBranchAddress("pfo_pidtof_p_at_calo", _data.pfo_pidtof_p_at_calo, &_branch.b_pfo_pidtof_p_at_calo);
    _hTree->SetBranchAddress("pfo_pidtof_closest_beta_0ps", _data.pfo_pidtof_closest_beta_0ps, &_branch.b_pfo_pidtof_closest_beta_0ps);
    _hTree->SetBranchAddress("pfo_pidtof_closest_beta_10ps", _data.pfo_pidtof_closest_beta_10ps, &_branch.b_pfo_pidtof_closest_beta_10ps);
    _hTree->SetBranchAddress("pfo_pidtof_closest_beta_50ps", _data.pfo_pidtof_closest_beta_50ps, &_branch.b_pfo_pidtof_closest_beta_50ps);
    _hTree->SetBranchAddress("pfo_pidtof_closest_beta_100ps", _data.pfo_pidtof_closest_beta_100ps, &_branch.b_pfo_pidtof_closest_beta_100ps);
    _hTree->SetBranchAddress("pfo_pidtof_fastest_beta_0ps", _data.pfo_pidtof_fastest_beta_0ps, &_branch.b_pfo_pidtof_fastest_beta_0ps);
    _hTree->SetBranchAddress("pfo_pidtof_fastest_beta_10ps", _data.pfo_pidtof_fastest_beta_10ps, &_branch.b_pfo_pidtof_fastest_beta_10ps);
    _hTree->SetBranchAddress("pfo_pidtof_fastest_beta_50ps", _data.pfo_pidtof_fastest_beta_50ps, &_branch.b_pfo_pidtof_fastest_beta_50ps);
    _hTree->SetBranchAddress("pfo_pidtof_fastest_beta_100ps", _data.pfo_pidtof_fastest_beta_100ps, &_branch.b_pfo_pidtof_fastest_beta_100ps);
    _hTree->SetBranchAddress("pfo_pidtof_cylfit_beta_0ps", _data.pfo_pidtof_cylfit_beta_0ps, &_branch.b_pfo_pidtof_cylfit_beta_0ps);
    _hTree->SetBranchAddress("pfo_pidtof_cylfit_beta_10ps", _data.pfo_pidtof_cylfit_beta_10ps, &_branch.b_pfo_pidtof_cylfit_beta_10ps);
    _hTree->SetBranchAddress("pfo_pidtof_cylfit_beta_50ps", _data.pfo_pidtof_cylfit_beta_50ps, &_branch.b_pfo_pidtof_cylfit_beta_50ps);
    _hTree->SetBranchAddress("pfo_pidtof_cylfit_beta_100ps", _data.pfo_pidtof_cylfit_beta_100ps, &_branch.b_pfo_pidtof_cylfit_beta_100ps);
    _hTree->SetBranchAddress("pfo_pidtof_closestfit_beta_0ps", _data.pfo_pidtof_closestfit_beta_0ps, &_branch.b_pfo_pidtof_closestfit_beta_0ps);
    _hTree->SetBranchAddress("pfo_pidtof_closestfit_beta_10ps", _data.pfo_pidtof_closestfit_beta_10ps, &_branch.b_pfo_pidtof_closestfit_beta_10ps);
    _hTree->SetBranchAddress("pfo_pidtof_closestfit_beta_50ps", _data.pfo_pidtof_closestfit_beta_50ps, &_branch.b_pfo_pidtof_closestfit_beta_50ps);
    _hTree->SetBranchAddress("pfo_pidtof_closestfit_beta_100ps", _data.pfo_pidtof_closestfit_beta_100ps, &_branch.b_pfo_pidtof_closestfit_beta_100ps);

}