#include "TreeWriter.hh"
using std::string;
using std::vector;

TreeWriter:: TreeWriter() {}
void TreeWriter::InitializeLPFOTree(TTree * _hTree, Tree_SSbar& data)
{
  _hTree->Branch("lpfo_match", data.lpfo_match, "lpfo_match[2]/I");
  _hTree->Branch("lpfo_truejet_pdg", data.lpfo_truejet_pdg, "lpfo_truejet_pdg[2]/I");
  _hTree->Branch("lpfo_truejet_type", data.lpfo_truejet_type, "lpfo_truejet_type[2]/I");
  _hTree->Branch("lpfo_pdgcheat", data.lpfo_pdgcheat , "lpfo_pdgcheat[2]/I");
  _hTree->Branch("lpfo_nparents", data.lpfo_nparents , "lpfo_nparents[2]/I");
  _hTree->Branch("lpfo_pdgcheat_parent", data.lpfo_pdgcheat_parent, "lpfo_pdgcheat_parent[2][1000]/I");
  _hTree->Branch("lpfo_E", data.lpfo_E, "lpfo_E[2]/F");
  _hTree->Branch("lpfo_px", data.lpfo_px, "lpfo_px[2]/F");
  _hTree->Branch("lpfo_py", data.lpfo_py, "lpfo_py[2]/F");
  _hTree->Branch("lpfo_pz", data.lpfo_pz, "lpfo_pz[2]/F");
  _hTree->Branch("lpfo_m", data.lpfo_m, "lpfo_m[2]/F");
  _hTree->Branch("lpfo_type", data.lpfo_type, "lpfo_type[2]/I");
  _hTree->Branch("lpfo_isoverlay", data.lpfo_isoverlay, "lpfo_isoverlay[2]/I");
  _hTree->Branch("lpfo_isisr", data.lpfo_isisr, "lpfo_isisr[2]/I");
  _hTree->Branch("lpfo_vtx", data.lpfo_vtx , "lpfo_vtx[2]/I");
  _hTree->Branch("lpfo_charge", data.lpfo_charge, "lpfo_charge[2]/I");
  _hTree->Branch("lpfo_ntracks", data.lpfo_ntracks, "lpfo_ntracks[2]/I");
  _hTree->Branch("lpfo_tpc_hits", data.lpfo_tpc_hits , "lpfo_tpc_hits[2]/I");
  _hTree->Branch("lpfo_dedx", data.lpfo_dedx, "lpfo_dedx[2]/F");
  _hTree->Branch("lpfo_dedxerror", data.lpfo_dedxerror, "lpfo_dedxerror[2]/F");
  _hTree->Branch("lpfo_d0", data.lpfo_d0, "lpfo_d0[2]/F");
  _hTree->Branch("lpfo_d0error", data.lpfo_d0error, "lpfo_d0error[2]/F");
  _hTree->Branch("lpfo_z0", data.lpfo_z0, "lpfo_z0[2]/F");
  _hTree->Branch("lpfo_z0error", data.lpfo_z0error, "lpfo_z0error[2]/F");
  _hTree->Branch("lpfo_phi", data.lpfo_phi , "lpfo_phi[2]/F");
  _hTree->Branch("lpfo_phierror", data.lpfo_phierror , "lpfo_phierror[2]/F");
  _hTree->Branch("lpfo_omega", data.lpfo_omega, "lpfo_omega[2]/F");
  _hTree->Branch("lpfo_omegaerror", data.lpfo_omegaerror, "lpfo_omegaerror[2]/F");
  _hTree->Branch("lpfo_tanlambda", data.lpfo_tanlambda, "lpfo_tanlambda[2]/F");
  _hTree->Branch("lpfo_tanlambdaerror", data.lpfo_tanlambdaerror, "lpfo_tanlambdaerror[2]/F");
  _hTree->Branch("lpfo_chi2", data.lpfo_chi2, "lpfo_chi2[2]/F");
  _hTree->Branch("lpfo_ndf", data.lpfo_ndf , "lpfo_ndf[2]/F");
  _hTree->Branch("lpfo_vtxpt", data.lpfo_vtxpt, "lpfo_vtxpt[2][3]/F");
  _hTree->Branch("lpfo_endpt", data.lpfo_endpt, "lpfo_endpt[2][3]/F");
  _hTree->Branch("lpfo_pid", data.lpfo_pid , "lpfo_pid[2]/I");
  _hTree->Branch("lpfo_pid_likelihood", data.lpfo_pid_likelihood, "lpfo_pid_likelihood[2]/F");
  _hTree->Branch("lpfo_pid_eprob", data.lpfo_pid_eprob, "lpfo_pid_eprob[2]/F");
  _hTree->Branch("lpfo_pid_muprob", data.lpfo_pid_muprob, "lpfo_pid_muprob[2]/F");
  _hTree->Branch("lpfo_pid_piprob", data.lpfo_pid_piprob, "lpfo_pid_piprob[2]/F");
  _hTree->Branch("lpfo_pid_kprob", data.lpfo_pid_kprob, "lpfo_pid_kprob[2]/F");
  _hTree->Branch("lpfo_pid_pprob", data.lpfo_pid_pprob, "lpfo_pid_pprob[2]/F");
  _hTree->Branch("lpfo_pid_hprob", data.lpfo_pid_hprob, "lpfo_pid_hprob[2]/F");
  _hTree->Branch("lpfo_piddedx", data.lpfo_piddedx, "lpfo_piddedx[2]/I");
  _hTree->Branch("lpfo_piddedx_likelihood", data.lpfo_piddedx_likelihood , "lpfo_piddedx_likelihood[2]/F");
  _hTree->Branch("lpfo_piddedx_eprob", data.lpfo_piddedx_eprob , "lpfo_piddedx_eprob[2]/F");
  _hTree->Branch("lpfo_piddedx_muprob", data.lpfo_piddedx_muprob, "lpfo_piddedx_muprob[2]/F");
  _hTree->Branch("lpfo_piddedx_piprob", data.lpfo_piddedx_piprob, "lpfo_piddedx_piprob[2]/F");
  _hTree->Branch("lpfo_piddedx_kprob", data.lpfo_piddedx_kprob , "lpfo_piddedx_kprob[2]/F");
  _hTree->Branch("lpfo_piddedx_pprob", data.lpfo_piddedx_pprob , "lpfo_piddedx_pprob[2]/F");
  _hTree->Branch("lpfo_piddedx_hprob", data.lpfo_piddedx_hprob , "lpfo_piddedx_hprob[2]/F");
  _hTree->Branch("lpfo_piddedx_e_dedxdist", data.lpfo_piddedx_e_dedxdist , "lpfo_piddedx_e_dedxdist[2]/F");
  _hTree->Branch("lpfo_piddedx_mu_dedxdist", data.lpfo_piddedx_mu_dedxdist, "lpfo_piddedx_mu_dedxdist[2]/F");
  _hTree->Branch("lpfo_piddedx_pi_dedxdist", data.lpfo_piddedx_pi_dedxdist, "lpfo_piddedx_pi_dedxdist[2]/F");
  _hTree->Branch("lpfo_piddedx_k_dedxdist", data.lpfo_piddedx_k_dedxdist , "lpfo_piddedx_k_dedxdist[2]/F");
  _hTree->Branch("lpfo_piddedx_p_dedxdist", data.lpfo_piddedx_p_dedxdist , "lpfo_piddedx_p_dedxdist[2]/F");
  _hTree->Branch("lpfo_piddedx_e_lkhood", data.lpfo_piddedx_e_lkhood, "lpfo_piddedx_e_lkhood[2]/F");
  _hTree->Branch("lpfo_piddedx_mu_lkhood", data.lpfo_piddedx_mu_lkhood, "lpfo_piddedx_mu_lkhood[2]/F");
  _hTree->Branch("lpfo_piddedx_pi_lkhood", data.lpfo_piddedx_pi_lkhood, "lpfo_piddedx_pi_lkhood[2]/F");
  _hTree->Branch("lpfo_piddedx_k_lkhood", data.lpfo_piddedx_k_lkhood, "lpfo_piddedx_k_lkhood[2]/F");
  _hTree->Branch("lpfo_piddedx_p_lkhood", data.lpfo_piddedx_p_lkhood, "lpfo_piddedx_p_lkhood[2]/F");
  _hTree->Branch("lpfo_pidtof_p_at_calo", data.lpfo_pidtof_p_at_calo, "lpfo_pidtof_p_at_calo[2]/F");
  _hTree->Branch("lpfo_pidtof_closest_beta_0ps", data.lpfo_pidtof_closest_beta_0ps , "lpfo_pidtof_closest_beta_0ps[2]/F");
  _hTree->Branch("lpfo_pidtof_closest_beta_10ps", data.lpfo_pidtof_closest_beta_10ps, "lpfo_pidtof_closest_beta_10ps[2]/F");
  _hTree->Branch("lpfo_pidtof_closest_beta_50ps", data.lpfo_pidtof_closest_beta_50ps, "lpfo_pidtof_closest_beta_50ps[2]/F");
  _hTree->Branch("lpfo_pidtof_closest_beta_100ps", data.lpfo_pidtof_closest_beta_100ps, "lpfo_pidtof_closest_beta_100ps[2]/F");
  _hTree->Branch("lpfo_pidtof_fastest_beta_0ps", data.lpfo_pidtof_fastest_beta_0ps , "lpfo_pidtof_fastest_beta_0ps[2]/F");
  _hTree->Branch("lpfo_pidtof_fastest_beta_10ps", data.lpfo_pidtof_fastest_beta_10ps, "lpfo_pidtof_fastest_beta_10ps[2]/F");
  _hTree->Branch("lpfo_pidtof_fastest_beta_50ps", data.lpfo_pidtof_fastest_beta_50ps, "lpfo_pidtof_fastest_beta_50ps[2]/F");
  _hTree->Branch("lpfo_pidtof_fastest_beta_100ps", data.lpfo_pidtof_fastest_beta_100ps, "lpfo_pidtof_fastest_beta_100ps[2]/F");
  _hTree->Branch("lpfo_pidtof_cylfit_beta_0ps", data.lpfo_pidtof_cylfit_beta_0ps, "lpfo_pidtof_cylfit_beta_0ps[2]/F");
  _hTree->Branch("lpfo_pidtof_cylfit_beta_10ps", data.lpfo_pidtof_cylfit_beta_10ps , "lpfo_pidtof_cylfit_beta_10ps[2]/F");
  _hTree->Branch("lpfo_pidtof_cylfit_beta_50ps", data.lpfo_pidtof_cylfit_beta_50ps , "lpfo_pidtof_cylfit_beta_50ps[2]/F");
  _hTree->Branch("lpfo_pidtof_cylfit_beta_100ps", data.lpfo_pidtof_cylfit_beta_100ps, "lpfo_pidtof_cylfit_beta_100ps[2]/F");
  _hTree->Branch("lpfo_pidtof_closestfit_beta_0ps", data.lpfo_pidtof_closestfit_beta_0ps, "lpfo_pidtof_closestfit_beta_0ps[2]/F");
  _hTree->Branch("lpfo_pidtof_closestfit_beta_10ps", data.lpfo_pidtof_closestfit_beta_10ps, "lpfo_pidtof_closestfit_beta_10ps[2]/F");
  _hTree->Branch("lpfo_pidtof_closestfit_beta_50ps", data.lpfo_pidtof_closestfit_beta_50ps, "lpfo_pidtof_closestfit_beta_50ps[2]/F");
  _hTree->Branch("lpfo_pidtof_closestfit_beta_100ps", data.lpfo_pidtof_closestfit_beta_100ps , "lpfo_pidtof_closestfit_beta_100ps[2]/F");

}

void TreeWriter::WriteLPFOVariables(PFOTools pt, PFO_QQbar *pqq, TreeVariables *data)
{
  int iLPFOs[2] = {pt.LPFO[0].ipfo, pt.LPFO[1].ipfo};
  int i = 0;
  for (auto iLPFO : iLPFOs){
    data->lpfo_match                       [i] = pqq->pfo_match                       [iLPFO];
    data->lpfo_truejet_pdg                 [i] = pqq->pfo_truejet_pdg                 [iLPFO];
    data->lpfo_truejet_type                [i] = pqq->pfo_truejet_type                [iLPFO];
    data->lpfo_pdgcheat                    [i] = pqq->pfo_pdgcheat                    [iLPFO];
    data->lpfo_nparents                    [i] = pqq->pfo_nparents                    [iLPFO];
    for (int j=0; j<pqq->pfo_nparents[iLPFO]; j++)
   {data->lpfo_pdgcheat_parent          [i][j] = pqq->pfo_pdgcheat_parent             [iLPFO][j];}
    data->lpfo_E                           [i] = pqq->pfo_E                           [iLPFO];
    data->lpfo_px                          [i] = pqq->pfo_px                          [iLPFO];
    data->lpfo_py                          [i] = pqq->pfo_py                          [iLPFO];
    data->lpfo_pz                          [i] = pqq->pfo_pz                          [iLPFO];
    data->lpfo_m                           [i] = pqq->pfo_m                           [iLPFO];
    data->lpfo_type                        [i] = pqq->pfo_type                        [iLPFO];
    data->lpfo_isoverlay                   [i] = pqq->pfo_isoverlay                   [iLPFO];
    data->lpfo_isisr                       [i] = pqq->pfo_isisr                       [iLPFO];
    data->lpfo_vtx                         [i] = pqq->pfo_vtx                         [iLPFO];
    data->lpfo_charge                      [i] = pqq->pfo_charge                      [iLPFO];
    data->lpfo_ntracks                     [i] = pqq->pfo_ntracks                     [iLPFO];
    data->lpfo_tpc_hits                    [i] = pqq->pfo_tpc_hits                    [iLPFO];
    data->lpfo_dedx                        [i] = pqq->pfo_dedx                        [iLPFO];
    data->lpfo_dedxerror                   [i] = pqq->pfo_dedxerror                   [iLPFO];
    data->lpfo_d0                          [i] = pqq->pfo_d0                          [iLPFO];
    data->lpfo_d0error                     [i] = pqq->pfo_d0error                     [iLPFO];
    data->lpfo_z0                          [i] = pqq->pfo_z0                          [iLPFO];
    data->lpfo_z0error                     [i] = pqq->pfo_z0error                     [iLPFO];
    data->lpfo_phi                         [i] = pqq->pfo_phi                         [iLPFO];
    data->lpfo_phierror                    [i] = pqq->pfo_phierror                    [iLPFO];
    data->lpfo_omega                       [i] = pqq->pfo_omega                       [iLPFO];
    data->lpfo_omegaerror                  [i] = pqq->pfo_omegaerror                  [iLPFO];
    data->lpfo_tanlambda                   [i] = pqq->pfo_tanlambda                   [iLPFO];
    data->lpfo_tanlambdaerror              [i] = pqq->pfo_tanlambdaerror              [iLPFO];
    data->lpfo_chi2                        [i] = pqq->pfo_chi2                        [iLPFO];
    data->lpfo_ndf                         [i] = pqq->pfo_ndf                         [iLPFO];
    for (int j=0; j<3; j++)
   {data->lpfo_vtxpt                    [i][j] = pqq->pfo_vtxpt                    [iLPFO][j];
    data->lpfo_endpt                    [i][j] = pqq->pfo_endpt                    [iLPFO][j];}
    data->lpfo_pid                         [i] = pqq->pfo_pid                         [iLPFO];
    data->lpfo_pid_likelihood              [i] = pqq->pfo_pid_likelihood              [iLPFO];
    data->lpfo_pid_eprob                   [i] = pqq->pfo_pid_eprob                   [iLPFO];
    data->lpfo_pid_muprob                  [i] = pqq->pfo_pid_muprob                  [iLPFO];
    data->lpfo_pid_piprob                  [i] = pqq->pfo_pid_piprob                  [iLPFO];
    data->lpfo_pid_kprob                   [i] = pqq->pfo_pid_kprob                   [iLPFO];
    data->lpfo_pid_pprob                   [i] = pqq->pfo_pid_pprob                   [iLPFO];
    data->lpfo_pid_hprob                   [i] = pqq->pfo_pid_hprob                   [iLPFO];
    data->lpfo_piddedx                     [i] = pqq->pfo_piddedx                     [iLPFO];
    data->lpfo_piddedx_likelihood          [i] = pqq->pfo_piddedx_likelihood          [iLPFO];
    data->lpfo_piddedx_eprob               [i] = pqq->pfo_piddedx_eprob               [iLPFO];
    data->lpfo_piddedx_muprob              [i] = pqq->pfo_piddedx_muprob              [iLPFO];
    data->lpfo_piddedx_piprob              [i] = pqq->pfo_piddedx_piprob              [iLPFO];
    data->lpfo_piddedx_kprob               [i] = pqq->pfo_piddedx_kprob               [iLPFO];
    data->lpfo_piddedx_pprob               [i] = pqq->pfo_piddedx_pprob               [iLPFO];
    data->lpfo_piddedx_hprob               [i] = pqq->pfo_piddedx_hprob               [iLPFO];
    data->lpfo_piddedx_e_dedxdist          [i] = pqq->pfo_piddedx_e_dedxdist          [iLPFO];
    data->lpfo_piddedx_mu_dedxdist         [i] = pqq->pfo_piddedx_mu_dedxdist         [iLPFO];
    data->lpfo_piddedx_pi_dedxdist         [i] = pqq->pfo_piddedx_pi_dedxdist         [iLPFO];
    data->lpfo_piddedx_k_dedxdist          [i] = pqq->pfo_piddedx_k_dedxdist          [iLPFO];
    data->lpfo_piddedx_p_dedxdist          [i] = pqq->pfo_piddedx_p_dedxdist          [iLPFO];
    data->lpfo_piddedx_e_lkhood            [i] = pqq->pfo_piddedx_e_lkhood            [iLPFO];
    data->lpfo_piddedx_mu_lkhood           [i] = pqq->pfo_piddedx_mu_lkhood           [iLPFO];
    data->lpfo_piddedx_pi_lkhood           [i] = pqq->pfo_piddedx_pi_lkhood           [iLPFO];
    data->lpfo_piddedx_k_lkhood            [i] = pqq->pfo_piddedx_k_lkhood            [iLPFO];
    data->lpfo_piddedx_p_lkhood            [i] = pqq->pfo_piddedx_p_lkhood            [iLPFO];
    data->lpfo_pidtof_p_at_calo            [i] = pqq->pfo_pidtof_p_at_calo            [iLPFO];
    data->lpfo_pidtof_closest_beta_0ps     [i] = pqq->pfo_pidtof_closest_beta_0ps     [iLPFO];
    data->lpfo_pidtof_closest_beta_10ps    [i] = pqq->pfo_pidtof_closest_beta_10ps    [iLPFO];
    data->lpfo_pidtof_closest_beta_50ps    [i] = pqq->pfo_pidtof_closest_beta_50ps    [iLPFO];
    data->lpfo_pidtof_closest_beta_100ps   [i] = pqq->pfo_pidtof_closest_beta_100ps   [iLPFO];
    data->lpfo_pidtof_fastest_beta_0ps     [i] = pqq->pfo_pidtof_fastest_beta_0ps     [iLPFO];
    data->lpfo_pidtof_fastest_beta_10ps    [i] = pqq->pfo_pidtof_fastest_beta_10ps    [iLPFO];
    data->lpfo_pidtof_fastest_beta_50ps    [i] = pqq->pfo_pidtof_fastest_beta_50ps    [iLPFO];
    data->lpfo_pidtof_fastest_beta_100ps   [i] = pqq->pfo_pidtof_fastest_beta_100ps   [iLPFO];
    data->lpfo_pidtof_cylfit_beta_0ps      [i] = pqq->pfo_pidtof_cylfit_beta_0ps      [iLPFO];
    data->lpfo_pidtof_cylfit_beta_10ps     [i] = pqq->pfo_pidtof_cylfit_beta_10ps     [iLPFO];
    data->lpfo_pidtof_cylfit_beta_50ps     [i] = pqq->pfo_pidtof_cylfit_beta_50ps     [iLPFO];
    data->lpfo_pidtof_cylfit_beta_100ps    [i] = pqq->pfo_pidtof_cylfit_beta_100ps    [iLPFO];
    data->lpfo_pidtof_closestfit_beta_0ps  [i] = pqq->pfo_pidtof_closestfit_beta_0ps  [iLPFO];
    data->lpfo_pidtof_closestfit_beta_10ps [i] = pqq->pfo_pidtof_closestfit_beta_10ps [iLPFO];
    data->lpfo_pidtof_closestfit_beta_50ps [i] = pqq->pfo_pidtof_closestfit_beta_50ps [iLPFO];
    data->lpfo_pidtof_closestfit_beta_100ps[i] = pqq->pfo_pidtof_closestfit_beta_100ps[iLPFO];
    i++;
  }

}