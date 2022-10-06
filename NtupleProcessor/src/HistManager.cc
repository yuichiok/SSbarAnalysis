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
  // TH1F
    h1[gen_q_cos]       = new TH1F("h_gen_q_cos","; Generated q#bar{q} cos#theta; Entries",100,-1,1);
    h1[gen_q_qcos]      = new TH1F("h_gen_q_qcos","; Generated q#bar{q} qcos#theta; Entries",100,-1,1);

    h1[gen_K_cos]       = new TH1F("h_gen_K_cos","; Generated Kaon cos#theta; Entries",100,-1,1);
    h1[gen_K_qcos]      = new TH1F("h_gen_K_qcos","; Generated Kaon qcos#theta; Entries",100,-1,1);
    
    h1[reco_K_cos]      = new TH1F("h_reco_K_cos","; LPFO Kaon cos#theta (no filter); Entries",100,-1,1);
    h1[reco_K_qcos]     = new TH1F("h_reco_K_qcos","; LPFO Kaon cos#theta (no filter); Entries",100,-1,1);

    // ISR parameters
    h1[reco_sum_jetE]   = new TH1F("h_reco_sum_jetE", ";Visible Energy (GeV);", 100, 0, 300);
    h1[reco_jet_sep]    = new TH1F("h_reco_jet_sep", ";Jet sep |cos#theta|;", 100, 0, 1);

    h1[lpfo_gen_K_mom]  = new TH1F("h_lpfo_gen_K_mom","; Leading Gen Kaon momentum (GeV); Entries",100,0,100);
    h1[lpfo_reco_K_mom] = new TH1F("h_lpfo_reco_K_mom","; Leading Reco Kaon momentum (GeV); Entries",100,0,100);

  // TH2F
    Int_t cos_bin = 100;
    h2[gen_N_K_cos]     = new TH2F("h2_gen_N_K_cos",";cos#theta;Generated N Kaons",cos_bin,-1,1,10,0,10);
    h2[reco_N_K_cos]    = new TH2F("h2_reco_N_K_cos",";cos#theta;Reconstructed N Kaons",   cos_bin,-1,1,10,0,10);
    h2[N_K_corr_cos]    = new TH2F("h2_N_K_corr_cos",";cos#theta;Correctly Reconstructed N Kaons",   cos_bin,-1,1,10,0,10);

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