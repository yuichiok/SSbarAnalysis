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
  Int_t   cos_bin = 100;
  Float_t bins_cos_fine[] = {-1.0,-0.98,-0.96,-0.94,-0.92,-0.90,-0.88,-0.86,-0.84,-0.82,-0.80,-0.75,-0.70,-0.60,-0.50,-0.40,-0.30,-0.20,-0.10,
                            0.0,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.75,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.0};
  Int_t   nbins_cos = sizeof(bins_cos_fine) / sizeof(Float_t) - 1;

  // TH1F
    h1[gen_q_cos]       = new TH1F("h_gen_q_cos","; Generated q#bar{q} cos#theta; Entries",nbins_cos,bins_cos_fine);
    h1[gen_q_qcos]      = new TH1F("h_gen_q_qcos","; Generated q#bar{q} qcos#theta; Entries",nbins_cos,bins_cos_fine);

    h1[gen_K_cos]       = new TH1F("h_gen_K_cos","; Generated Kaon cos#theta; Entries",nbins_cos,bins_cos_fine);
    h1[gen_K_qcos]      = new TH1F("h_gen_K_qcos","; Generated Kaon qcos#theta; Entries",nbins_cos,bins_cos_fine);
    
    h1[reco_K_cos]      = new TH1F("h_reco_K_cos",";LPFO Kaon cos#theta; Entries",nbins_cos,bins_cos_fine);
    h1[reco_K_qcos]     = new TH1F("h_reco_K_qcos",";LPFO Kaon cos#theta; Entries",nbins_cos,bins_cos_fine);
    
    h1[reco_K_scos]     = new TH1F("h_reco_K_scos",";LPFO Kaon cos#theta; Entries",nbins_cos,bins_cos_fine);

    h1[gen_reco_K_sep_cos] = new TH1F("h_gen_reco_K_sep_cos",";Kaon cos#theta_{gen-reco}; Entries",nbins_cos,bins_cos_fine);

  // ISR parameters
    h1[reco_sum_jetE]   = new TH1F("h_reco_sum_jetE", ";Visible Energy (GeV);", 100, 0, 300);
    h1[reco_jet_sep]    = new TH1F("h_reco_jet_sep", ";Jet sep |cos#theta|;", 100, 0, 1);

    h1[lpfo_gen_K_mom]  = new TH1F("h_lpfo_gen_K_mom","; Leading Gen Kaon momentum (GeV); Entries",100,0,100);
    h1[lpfo_reco_K_mom] = new TH1F("h_lpfo_reco_K_mom","; Leading Reco Kaon momentum (GeV); Entries",100,0,100);

  // Number of Gen Reco Kaons
    h1[gen_N_K_cos]     = new TH1F("h_gen_N_K_cos",";cos#theta;Generated N Kaons", nbins_cos,bins_cos_fine);
    h1[reco_N_K_cos]    = new TH1F("h_reco_N_K_cos",";cos#theta;Reconstructed N Kaons", nbins_cos,bins_cos_fine);
    h1[N_K_corr_cos]    = new TH1F("h_N_K_corr_cos",";cos#theta;Correctly Reconstructed N Kaons", nbins_cos,bins_cos_fine);

  // pq method
    h1_pq[acc_KK]      = new TH1F("h_acc_KK_cos",";Accepted K^{+}K^{-} cos#theta;N_{acc}",nbins_cos,bins_cos_fine);
    h1_pq[rej_KK]      = new TH1F("h_rej_KK_cos",";Rejected K^{+}K^{-} cos#theta;N_{rej}",nbins_cos,bins_cos_fine);

  // TH2F
    h2[gen_K_p_cos]  = new TH2F("h2_gen_K_p_cos" ,";cos#theta;p (GeV)",cos_bin,-1,1,100,0,60);
    h2[reco_K_p_cos] = new TH2F("h2_reco_K_p_cos",";cos#theta;p (GeV)",cos_bin,-1,1,100,0,60);
    h2[stable_cos]   = new TH2F("h2_stable_cos",";cos#theta;Stability",cos_bin,-1,1,50,0,1);
    h2[purity_cos]   = new TH2F("h2_purity_cos",";cos#theta;Purity",   cos_bin,-1,1,50,0,1);

    Hist2List();

}

void HistManager::Hist2List()
{
  for (auto ih : h1) {
    hList1->Add(ih);
  }

  for (auto ih : h1_pq) {
    hList1_pq->Add(ih);
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

  TDirectory * d_pq = output->mkdir("pq");
    d_pq->cd();
    hList1_pq->Write();
    output->cd();

}