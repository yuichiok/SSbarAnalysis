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
    h[gen_q_cos]   = new TH1F("h_gen_q_cos","; Generated q#bar{q} cos#theta; Entries",100,-1,1);
    h[gen_q_qcos]  = new TH1F("h_gen_q_qcos","; Generated q#bar{q} qcos#theta; Entries",100,-1,1);

    h[gen_K_cos]   = new TH1F("h_gen_K_cos","; Generated Kaon cos#theta; Entries",100,-1,1);
    h[gen_K_qcos]  = new TH1F("h_gen_K_qcos","; Generated Kaon qcos#theta; Entries",100,-1,1);
    
    h[reco_K_cos]  = new TH1F("h_reco_K_cos","; LPFO Kaon cos#theta (no filter); Entries",100,-1,1);
    h[reco_K_qcos] = new TH1F("h_reco_K_qcos","; LPFO Kaon cos#theta (no filter); Entries",100,-1,1);

    h[lpfo_gen_K_mom]  = new TH1F("h_lpfo_gen_K_mom","; Leading Gen Kaon momentum (GeV); Entries",100,0,100);
    h[lpfo_reco_K_mom] = new TH1F("h_lpfo_reco_K_mom","; Leading Reco Kaon momentum (GeV); Entries",100,0,100);

    Hist2List();

}

void HistManager::Hist2List()
{
  for (auto ih : h) {
    hList->Add(ih);
  }

}

void HistManager::WriteLists( TFile * output)
{
  // Focus to this file
  output->cd();

  // Write histogram lists
  hList->Write();

}