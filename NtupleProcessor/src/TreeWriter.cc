#include "TreeWriter.hh"
using std::string;
using std::vector;

namespace QQbarAnalysis
{
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
    /*
      _hTree->Branch("n_valid_pfo", &data.n_valid_pfo, "n_valid_pfo/I");
      _hTree->Branch("vpfo_E", data.vpfo_E, "vpfo_E[n_valid_pfo]/F");
      _hTree->Branch("vpfo_p", data.vpfo_p, "vpfo_p[n_valid_pfo]/F");
      _hTree->Branch("vpfo_cos", data.vpfo_cos, "vpfo_cos[n_valid_pfo]/F");
      _hTree->Branch("vpfo_dedx", data.vpfo_dedx, "vpfo_dedx[n_valid_pfo]/F");
      _hTree->Branch("vpfo_pdgcheat", data.vpfo_pdgcheat, "vpfo_pdgcheat[n_valid_pfo]/I");
      _hTree->Branch("vpfo_piddedx_k_dedxdist", data.vpfo_piddedx_k_dedxdist, "vpfo_piddedx_k_dedxdist[n_valid_pfo]/F");
      _hTree->Branch("vpfo_piddedx_pi_dedxdist", data.vpfo_piddedx_pi_dedxdist, "vpfo_piddedx_pi_dedxdist[n_valid_pfo]/F");
      _hTree->Branch("vpfo_piddedx_p_dedxdist", data.vpfo_piddedx_p_dedxdist, "vpfo_piddedx_p_dedxdist[n_valid_pfo]/F");
    */

      _hTree->Branch("LPFO_cos",   data.LPFO_cos,    "LPFO_cos[2]/F");
      _hTree->Branch("LPFO_qcos",  data.LPFO_qcos,  "LPFO_qcos[2]/F");
  }

}