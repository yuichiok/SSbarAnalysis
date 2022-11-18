#include "TreeWriter.hh"
using std::string;
using std::vector;

TreeWriter:: TreeWriter() {}
void TreeWriter::InitializeDataTree(TTree * _hTree, Tree_Data& data)
{
    _hTree->Branch("sum_jet_E",&data.sum_jet_E,"sum_jet_E/F");
    _hTree->Branch("jet_acol",&data.jet_acol,"jet_acol/F");

    _hTree->Branch("dEdx_pdg_match",   &data.dEdx_pdg_match,   "dEdx_pdg_match/I");

    _hTree->Branch("N_K_Gen",   &data.N_K_Gen,   "N_K_Gen/I");
    _hTree->Branch("N_K_PFO",   &data.N_K_PFO,   "N_K_PFO/I");
    _hTree->Branch("N_K_corr",  &data.N_K_corr,  "N_K_corr/I");
    _hTree->Branch("stability", &data.stability, "stability/F");
    _hTree->Branch("purity",    &data.purity,    "purity/F");

  // Valid PFO Collection
    _hTree->Branch("n_valid_pfo", &data.n_valid_pfo, "n_valid_pfo/I");
    _hTree->Branch("vpfo_E", data.vpfo_E, "vpfo_E[n_valid_pfo]/F");
    _hTree->Branch("vpfo_p", data.vpfo_p, "vpfo_p[n_valid_pfo]/F");
    _hTree->Branch("vpfo_cos", data.vpfo_cos, "vpfo_cos[n_valid_pfo]/F");
    _hTree->Branch("vpfo_dedx", data.vpfo_dedx, "vpfo_dedx[n_valid_pfo]/F");
    _hTree->Branch("vpfo_pdgcheat", data.vpfo_pdgcheat, "vpfo_pdgcheat[n_valid_pfo]/I");
    _hTree->Branch("vpfo_piddedx_k_dedxdist", data.vpfo_piddedx_k_dedxdist, "vpfo_piddedx_k_dedxdist[n_valid_pfo]/F");
    _hTree->Branch("vpfo_piddedx_pi_dedxdist", data.vpfo_piddedx_pi_dedxdist, "vpfo_piddedx_pi_dedxdist[n_valid_pfo]/F");
    _hTree->Branch("vpfo_piddedx_p_dedxdist", data.vpfo_piddedx_p_dedxdist, "vpfo_piddedx_p_dedxdist[n_valid_pfo]/F");

    _hTree->Branch("LPFO_cos",   data.LPFO_cos,    "LPFO_cos[2]/F");
    _hTree->Branch("LPFO_qcos",  data.LPFO_qcos,  "LPFO_qcos[2]/F");
}

void TreeWriter::WriteLPFO_Info(PFOTools pt, PFO_QQbar *pqq, TreeVariables *data)
{
  int iLPFOs[2] = {pt.LPFO[0].ipfo, pt.LPFO[1].ipfo};
  int i = 0;
  for (auto iLPFO : iLPFOs){
    data->lpfo_match                       [i] = pqq->pfo_match                       [iLPFO];
    data->lpfo_truejet_pdg                 [i] = pqq->pfo_truejet_pdg                 [iLPFO];
    data->lpfo_truejet_type                [i] = pqq->pfo_truejet_type                [iLPFO];
    data->lpfo_pdgcheat                    [i] = pqq->pfo_pdgcheat                    [iLPFO];
    data->lpfo_pdgcheat_id                 [i] = pqq->pfo_pdgcheat_id                 [iLPFO];
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