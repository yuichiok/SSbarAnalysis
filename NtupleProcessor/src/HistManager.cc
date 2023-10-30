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
  h_primary[static_cast<int>(h_primary::d0)] = new TH2F("d0_primary"  ,"Reconstructed primary K;d_0, [mm];ctag",1000,0,0.2,100,0,1);
  h_primary[static_cast<int>(h_primary::d0_sigma)] = new TH2F("d0_sigma_primary"  ,"Reconstructed primary K;#sigma(d_0), [mm];ctag",100,0,0.01,100,0,1);
  h_primary[static_cast<int>(h_primary::d0_sigma_d0)] = new TH2F("d0_sigma_d0_primary"  ,"Reconstructed primary K;#frac{d0}{#sigma(d0)};ctag",1000,0,100,100,0,1);

  h_primary[static_cast<int>(h_primary::z0)] = new TH2F("z0_primary"  ,"Reconstructed primary K;z_0, [mm];ctag",10000,0,2,100,0,1);
  h_primary[static_cast<int>(h_primary::z0_sigma)] = new TH2F("z0_sigma_primary"  ,"Reconstructed primary K;#sigma(z_0), [mm];ctag",1000,0,0.1,100,0,1);
  h_primary[static_cast<int>(h_primary::z0_sigma_z0)] = new TH2F("z0_sigma_z0_primary"  ,"Reconstructed primary K;#frac{z0}{#sigma(z0)};ctag",1000,0,100,100,0,1);

  h_secondary[static_cast<int>(h_secondary::d0)] = new TH2F("d0_secondary"  ,"Reconstructed secondary K;d_0, [mm];ctag",1000,0,0.2,100,0,1);
  h_secondary[static_cast<int>(h_secondary::d0_sigma)] = new TH2F("d0_sigma_secondary"  ,"Reconstructed secondary K;#sigma(d_0), [mm];ctag",100,0,0.01,100,0,1);
  h_secondary[static_cast<int>(h_secondary::d0_sigma_d0)] = new TH2F("d0_sigma_d0_secondary"  ,"Reconstructed secondary K;#frac{d0}{#sigma(d0)};ctag",1000,0,100,100,0,1);

  h_secondary[static_cast<int>(h_secondary::z0)] = new TH2F("z0_secondary"  ,"Reconstructed secondary K;z_0, [mm];ctag",10000,0,2,100,0,1);
  h_secondary[static_cast<int>(h_secondary::z0_sigma)] = new TH2F("z0_sigma_secondary"  ,"Reconstructed secondary K;#sigma(z_0), [mm];ctag",100,0,0.01,100,0,1);
  h_secondary[static_cast<int>(h_secondary::z0_sigma_z0)] = new TH2F("z0_sigma_z0_secondary"  ,"Reconstructed secondary K;#frac{z0}{#sigma(z0)};ctag",1000,0,100,100,0,1);

  h_garbage[static_cast<int>(h_garbage::d0)] = new TH2F("d0_garbage"  ,"Reconstructed garbage K;d_0, [mm];ctag",1000,0,0.2,100,0,1);
  h_garbage[static_cast<int>(h_garbage::d0_sigma)] = new TH2F("d0_sigma_garbage"  ,"Reconstructed garbage K;#sigma(d_0), [mm];ctag",100,0,0.1,100,0,1);
  h_garbage[static_cast<int>(h_garbage::d0_sigma_d0)] = new TH2F("d0_sigma_d0_garbage"  ,"Reconstructed garbage K;#frac{d0}{#sigma(d0)};ctag",1000,0,100,100,0,1);

  h_garbage[static_cast<int>(h_garbage::z0)] = new TH2F("z0_garbage"  ,"Reconstructed garbage K;z_0, [mm];ctag",10000,0,2,100,0,1);
  h_garbage[static_cast<int>(h_garbage::z0_sigma)] = new TH2F("z0_sigma_garbage"  ,"Reconstructed garbage K;#sigma(z_0), [mm];ctag",1000,0,0.1,100,0,1);
  h_garbage[static_cast<int>(h_garbage::z0_sigma_z0)] = new TH2F("z0_sigma_z0_garbage"  ,"Reconstructed garbage K;#frac{z0}{#sigma(z0)};ctag",1000,0,100,100,0,1);

  h_total[static_cast<int>(h_total::d0)] = new TH2F("d0_total"  ,"Reconstructed garbage K;d_0, [mm];ctag",1000,0,0.2,100,0,1);
  h_total[static_cast<int>(h_total::d0_sigma)] = new TH2F("d0_sigma_total"  ,"Reconstructed garbage K;#sigma(d_0), [mm];ctag",1000,0,0.1,100,0,1);
  h_total[static_cast<int>(h_total::d0_sigma_d0)] = new TH2F("d0_sigma_d0_total"  ,"Reconstructed garbage K;#frac{d0}{#sigma(d0)};ctag",1000,0,100,100,0,1);

  h_total[static_cast<int>(h_total::z0)] = new TH2F("z0_total"  ,"Reconstructed all K;z_0, [mm];ctag",10000,0,2,100,0,1);
  h_total[static_cast<int>(h_total::z0_sigma)] = new TH2F("z0_sigma__total"  ,"Reconstructed all K;#sigma(z_0), [mm];ctag",1000,0,0.1,100,0,1);
  h_total[static_cast<int>(h_total::z0_sigma_z0)] = new TH2F("z0_sigma_z0_total"  ,"Reconstructed all K;#frac{z0}{#sigma(z0)};ctag",1000,0,100,100,0,1);


  h_cos[static_cast<int>(h_cos::cos_theta)] = new TH2F("cos_theta"  ,"Reconstructed all K;cos(#theta);ctag",200,-1,1,100,0,1);
  h_cos[static_cast<int>(h_cos::cos_theta_acc)] = new TH2F("cos_theta_acc"  ,"Reconstructed all acc K;cos(#theta);ctag",200,-1,1,100,0,1);
  h_cos[static_cast<int>(h_cos::cos_theta_rej)] = new TH2F("cos_theta_rej"  ,"Reconstructed all rej K;cos(#theta);ctag",200,-1,1,100,0,1);
  h_cos[static_cast<int>(h_cos::cos_theta_primary)] = new TH2F("cos_theta_primary"  ,"Reconstructed primary K;cos(#theta);ctag",200,-1,1,100,0,1);
  h_cos[static_cast<int>(h_cos::cos_theta_acc_primary)] = new TH2F("cos_theta_acc_primary"  ,"Reconstructed primary acc K;cos(#theta);ctag",200,-1,1,100,0,1);
  h_cos[static_cast<int>(h_cos::cos_theta_rej_primary)] = new TH2F("cos_theta_rej_primary"  ,"Reconstructed primary rej K;cos(#theta);ctag",200,-1,1,100,0,1);


  h_general[static_cast<int>(h_general::cuts)] = new TH1F("cuts"  ,"Cuts stages",10,0,10);
  h_general[static_cast<int>(h_general::cuts)]->GetXaxis()->SetBinLabel(1,"LPFO_checks[Kaon]");
  h_general[static_cast<int>(h_general::cuts)]->GetXaxis()->SetBinLabel(2,"p_mag");
  h_general[static_cast<int>(h_general::cuts)]->GetXaxis()->SetBinLabel(3,"double_tag && is_ss");
  h_general[static_cast<int>(h_general::cuts)]->GetXaxis()->SetBinLabel(4,"sign_check");
  h_general[static_cast<int>(h_general::cuts)]->GetXaxis()->SetBinLabel(5,"P");
  h_general[static_cast<int>(h_general::cuts)]->GetXaxis()->SetBinLabel(6,"S");

  h_general[static_cast<int>(h_general::ctag_primary)] = new TH1F("ctag_primary"  ,"ctag of primary K",100,0,1);
  h_general[static_cast<int>(h_general::ctag_secondary)] = new TH1F("ctag_secondary"  ,"ctag of secondary K",100,0,1);
  h_general[static_cast<int>(h_general::ctag_total)] = new TH1F("ctag_total"  ,"ctag of all K",100,0,1);

  Hist2List();
}

void HistManager::Hist2List()
{  
  for (auto ih : h_primary) {
    hList_primary->Add(ih);
  }

  for (auto ih : h_secondary) {
    hList_secondary->Add(ih);
  }

  for (auto ih : h_garbage) {
    hList_garbage->Add(ih);
  }

  for (auto ih : h_total) {
    hList_total->Add(ih);
  }

  for (auto ih : h_general) {
    hList_general->Add(ih);
  }
  
  for (auto ih : h_cos) {
    hList_cos->Add(ih);
  }  
}

void HistManager::WriteLists( TFile * output)
{
  output->cd();

  hList_general->Write();
  hList_cos->Write();
  TDirectory * d_ctag = output->mkdir("ctag-analysis");
    d_ctag->cd();
    hList_primary->Write();
    hList_secondary->Write();
    hList_garbage->Write();
    hList_total->Write();
}