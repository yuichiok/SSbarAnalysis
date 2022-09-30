#include <iostream>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TFile.h> 

#include "../include/HistManager.hh"

using std::cout;   using std::endl;

HistManager::HistManager()
{
  // Initialization
    // InitializeHists();
    // Hist2List();
}

void HistManager::InitializeHists()
{
    h_lpfo_K_mom = new TH1F("h_lpfo_K_mom","; Leading Kaon momentum (GeV); Entries",100,0,100);

    Hist2List();

}

void HistManager::Hist2List()
{
    hList->Add(h_lpfo_K_mom);

}

void HistManager::WriteLists( TFile * output)
{
  // Focus to this file
  output->cd();

  // Write histogram lists
  hList->Write();

}