#include <iostream>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TFile.h> 

#include "../include/HistManager.hh"

using std::cout;   using std::endl;

HistManager::HistManager() {}

void HistManager::InitializeHists()
{
  Int_t cos_bin = 100;

  // TH1F
    h1[gen_q_cos]       = new TH1F("h_gen_q_cos","; Generated q#bar{q} cos#theta; Entries",cos_bin,-1,1);
    h1[gen_q_qcos]      = new TH1F("h_gen_q_qcos","; Generated q#bar{q} qcos#theta; Entries",cos_bin,-1,1);

    h1[gen_K_cos]       = new TH1F("h_gen_K_cos","; Generated Kaon cos#theta; Entries",cos_bin,-1,1);
    h1[gen_K_qcos]      = new TH1F("h_gen_K_qcos","; Generated Kaon qcos#theta; Entries",cos_bin,-1,1);
    
    h1[reco_K_cos]      = new TH1F("h_reco_K_cos","; LPFO Kaon cos#theta (no filter); Entries",cos_bin,-1,1);
    h1[reco_K_qcos]     = new TH1F("h_reco_K_qcos","; LPFO Kaon cos#theta (no filter); Entries",cos_bin,-1,1);

  // ISR parameters
    h1[reco_sum_jetE]   = new TH1F("h_reco_sum_jetE", ";Visible Energy (GeV);", 100, 0, 300);
    h1[reco_jet_sep]    = new TH1F("h_reco_jet_sep", ";Jet sep |cos#theta|;", 100, 0, 1);

    h1[lpfo_gen_K_mom]  = new TH1F("h_lpfo_gen_K_mom","; Leading Gen Kaon momentum (GeV); Entries",100,0,100);
    h1[lpfo_reco_K_mom] = new TH1F("h_lpfo_reco_K_mom","; Leading Reco Kaon momentum (GeV); Entries",100,0,100);

  // Number of Gen Reco Kaons
    h1[gen_N_K_cos]     = new TH1F("h_gen_N_K_cos",";cos#theta;Generated N Kaons",cos_bin,-1,1);
    h1[reco_N_K_cos]    = new TH1F("h_reco_N_K_cos",";cos#theta;Reconstructed N Kaons",   cos_bin,-1,1);
    h1[N_K_corr_cos]    = new TH1F("h_N_K_corr_cos",";cos#theta;Correctly Reconstructed N Kaons",   cos_bin,-1,1);

  // pq method
    h1_pq[acc_KK]      = new TH1F("h_acc_KK_cos",";Accepted K^{+}K^{-} cos#theta;N_{acc}",cos_bin,-1,1);
    h1_pq[rej_KK]      = new TH1F("h_rej_KK_cos",";Rejected K^{+}K^{-} cos#theta;N_{rej}",cos_bin,-1,1);
    h1_pq[p_KK]        = new TH1F("h_p_KK_cos",";K^{+}K^{-} cos#theta;p value",cos_bin,-1,1);

  // TH2F
    h2[stable_cos]      = new TH2F("h2_stable_cos",";cos#theta;Stability",cos_bin,-1,1,50,0,1);
    h2[purity_cos]      = new TH2F("h2_purity_cos",";cos#theta;Purity",   cos_bin,-1,1,50,0,1);

    Hist2List();

}

void HistManager::Hist2List()
{
  for (auto ih : h1) {
    hList1->Add(ih);
  }

  for (auto ih : h2) {
    hList2->Add(ih);
  }
}

void HistManager::WriteLists( TFile * output)
{
  // Focus to this file
  output->cd();

  // Write histogram lists
  hList1->Write();
  hList2->Write();

}