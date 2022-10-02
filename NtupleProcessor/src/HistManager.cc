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
    h_gen_K_cos  = new TH1F("h_gen_K_cos","; Generated Kaon cos#theta; Entries",100,-1,1);
    h_reco_K_cos = new TH1F("h_lpfo_reco_K_cos","; LPFO Kaon cos#theta (no filter); Entries",100,-1,1);

    h_lpfo_gen_K_mom  = new TH1F("h_lpfo_gen_K_mom","; Leading Gen Kaon momentum (GeV); Entries",100,0,100);
    h_lpfo_reco_K_mom = new TH1F("h_lpfo_reco_K_mom","; Leading Reco Kaon momentum (GeV); Entries",100,0,100);

    Hist2List();

}

void HistManager::Hist2List()
{
    hList->Add(h_gen_K_cos);
    hList->Add(h_reco_K_cos);

    hList->Add(h_lpfo_gen_K_mom);
    hList->Add(h_lpfo_reco_K_mom);
}

void HistManager::WriteLists( TFile * output)
{
  // Focus to this file
  output->cd();

  // Write histogram lists
  hList->Write();

}