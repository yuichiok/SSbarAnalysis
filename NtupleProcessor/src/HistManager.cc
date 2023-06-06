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

  // Own

    h1_K_reco[d0_K_reco_primary] = new TH1F("d0_K_reco_primary"  ,"Reconstructed primary K;d_0, [mm];",2000,0,2);
    h1_K_reco[d0_sigma_K_reco_primary] = new TH1F("d0_sigma_K_reco_primary"  ,"Reconstructed primary K;#sigma(d_0), [mm];",1000,0,0.1);
    h1_K_reco[d0_sigma_d0_K_reco_primary] = new TH1F("d0_sigma_d0_K_reco_primary"  ,"Reconstructed primary K;#frac{d0}{#sigma(d0)};",1000,0,100);

    h1_K_reco[z0_K_reco_primary] = new TH1F("z0_K_reco_primary"  ,"Reconstructed primary K;z_0, [mm];",2000,0,2);
    h1_K_reco[z0_sigma_K_reco_primary] = new TH1F("z0_sigma_K_reco_primary"  ,"Reconstructed primary K;#sigma(z_0), [mm];",1000,0,0.1);
    h1_K_reco[z0_sigma_z0_K_reco_primary] = new TH1F("z0_sigma_z0_K_reco_primary"  ,"Reconstructed primary K;#frac{z0}{#sigma(z0)};",1000,0,100);

    h1_K_reco[d0_K_reco_secondary] = new TH1F("d0_K_reco_secondary"  ,"Reconstructed secondary K;d_0, [mm];",2000,0,2);
    h1_K_reco[d0_sigma_K_reco_secondary] = new TH1F("d0_sigma_K_reco_secondary"  ,"Reconstructed secondary K;#sigma(d_0), [mm];",1000,0,0.1);
    h1_K_reco[d0_sigma_d0_K_reco_secondary] = new TH1F("d0_sigma_d0_K_reco_secondary"  ,"Reconstructed secondary K;#frac{d0}{#sigma(d0)};",1000,0,100);

    h1_K_reco[z0_K_reco_secondary] = new TH1F("z0_K_reco_secondary"  ,"Reconstructed secondary K;z_0, [um];",2000,0,2);
    h1_K_reco[z0_sigma_K_reco_secondary] = new TH1F("z0_sigma_K_reco_secondary"  ,"Reconstructed secondary K;#sigma(z_0), [mm];",1000,0,0.1);
    h1_K_reco[z0_sigma_z0_K_reco_secondary] = new TH1F("z0_sigma_z0_K_reco_secondary"  ,"Reconstructed secondary K;#frac{z0}{#sigma(z0)};",1000,0,100);

    h1_K_reco[d0_K_reco_pseudo] = new TH1F("d0_K_reco_pseudo"  ,"Reconstructed pseudo K;d_0, [mm];",2000,0,2);
    h1_K_reco[d0_sigma_K_reco_pseudo] = new TH1F("d0_sigma_K_reco_pseudo"  ,"Reconstructed pseudo K;#sigma(d_0), [mm];",1000,0,0.1);
    h1_K_reco[d0_sigma_d0_K_reco_pseudo] = new TH1F("d0_sigma_d0_K_reco_pseudo"  ,"Reconstructed pseudo K;#frac{d0}{#sigma(d0)};",1000,0,100);

    h1_K_reco[z0_K_reco_pseudo] = new TH1F("z0_K_reco_pseudo"  ,"Reconstructed pseudo K;z_0, [um];",2000,0,2);
    h1_K_reco[z0_sigma_K_reco_pseudo] = new TH1F("z0_sigma_K_reco_pseudo"  ,"Reconstructed pseudo K;#sigma(z_0), [mm];",1000,0,0.1);
    h1_K_reco[z0_sigma_z0_K_reco_pseudo] = new TH1F("z0_sigma_z0_K_reco_pseudo"  ,"Reconstructed pseudo K;#frac{z0}{#sigma(z0)};",1000,0,100);

    h1_K_reco[d0_K_reco_garbage] = new TH1F("d0_K_reco_garbage"  ,"Reconstructed garbage K;d_0, [mm];",2000,0,2);
    h1_K_reco[d0_sigma_K_reco_garbage] = new TH1F("d0_sigma_K_reco_garbage"  ,"Reconstructed garbage K;#sigma(d_0), [mm];",1000,0,0.1);
    h1_K_reco[d0_sigma_d0_K_reco_garbage] = new TH1F("d0_sigma_d0_K_reco_garbage"  ,"Reconstructed garbage K;#frac{d0}{#sigma(d0)};",1000,0,100);

    h1_K_reco[z0_K_reco_garbage] = new TH1F("z0_K_reco_garbage"  ,"Reconstructed garbage K;z_0, [um];",2000,0,2);
    h1_K_reco[z0_sigma_K_reco_garbage] = new TH1F("z0_sigma_K_reco_garbage"  ,"Reconstructed garbage K;#sigma(z_0), [mm];",1000,0,0.1);
    h1_K_reco[z0_sigma_z0_K_reco_garbage] = new TH1F("z0_sigma_z0_K_reco_garbage"  ,"Reconstructed garbage K;#frac{z0}{#sigma(z0)};",1000,0,100);


    h1_K_reco[d0_K_reco_primary_ctag_cut] = new TH1F("d0_K_reco_primary_ctag_cut"  ,"Reconstructed primary K with ctag cut;d_0, [mm];",2000,0,2);
    h1_K_reco[d0_sigma_K_reco_primary_ctag_cut] = new TH1F("d0_sigma_K_reco_primary_ctag_cut"  ,"Reconstructed primary K ctag cut;#sigma(d_0), [mm];",1000,0,0.1);
    h1_K_reco[d0_sigma_d0_K_reco_primary_ctag_cut] = new TH1F("d0_sigma_d0_K_reco_primary_ctag_cut"  ,"Reconstructed primary K ctag cut;#frac{d0}{#sigma(d0)};",1000,0,100);

    h1_K_reco[z0_K_reco_primary_ctag_cut] = new TH1F("z0_K_reco_primary_ctag_cut"  ,"Reconstructed primary K; ctag cutz_0, [mm];",2000,0,2);
    h1_K_reco[z0_sigma_K_reco_primary_ctag_cut] = new TH1F("z0_sigma_K_reco_primary_ctag_cut"  ,"Reconstructed primary K ctag cut;#sigma(z_0), [mm];",1000,0,0.1);
    h1_K_reco[z0_sigma_z0_K_reco_primary_ctag_cut] = new TH1F("z0_sigma_z0_K_reco_primary_ctag_cut"  ,"Reconstructed primary K ctag cut;#frac{z0}{#sigma(z0)};",1000,0,100);

    h1_K_reco[d0_K_reco_secondary_ctag_cut] = new TH1F("d0_K_reco_secondary_ctag_cut"  ,"Reconstructed secondary K ctag cut;d_0, [mm];",2000,0,2);
    h1_K_reco[d0_sigma_K_reco_secondary_ctag_cut] = new TH1F("d0_sigma_K_reco_secondary_ctag_cut"  ,"Reconstructed secondary K ctag cut;#sigma(d_0), [mm];",1000,0,0.1);
    h1_K_reco[d0_sigma_d0_K_reco_secondary_ctag_cut] = new TH1F("d0_sigma_d0_K_reco_secondary_ctag_cut"  ,"Reconstructed secondary K ctag cut;#frac{d0}{#sigma(d0)};",1000,0,100);

    h1_K_reco[z0_K_reco_secondary_ctag_cut] = new TH1F("z0_K_reco_secondary_ctag_cut"  ,"Reconstructed secondary K ctag cut;z_0, [um];",2000,0,2);
    h1_K_reco[z0_sigma_K_reco_secondary_ctag_cut] = new TH1F("z0_sigma_K_reco_secondary_ctag_cut"  ,"Reconstructed secondary K ctag cut;#sigma(z_0), [mm];",1000,0,0.1);
    h1_K_reco[z0_sigma_z0_K_reco_secondary_ctag_cut] = new TH1F("z0_sigma_z0_K_reco_secondary_ctag_cut"  ,"Reconstructed secondary K ctag cut;#frac{z0}{#sigma(z0)};",1000,0,100);

    h1_K_reco[cuts] = new TH1F("cuts"  ,"Cuts stages",10,0,10);
    h1_K_reco[cuts]->GetXaxis()->SetBinLabel(1,"No cuts");
    h1_K_reco[cuts]->GetXaxis()->SetBinLabel(2,"LPFO_checks && double_tag && is_ss");
    h1_K_reco[cuts]->GetXaxis()->SetBinLabel(3,"sign_check");
    h1_K_reco[cuts]->GetXaxis()->SetBinLabel(4,"P");
    h1_K_reco[cuts]->GetXaxis()->SetBinLabel(5,"S");
    

    h1_K_reco[pmag_K_reco] = new TH1F("pmag_K_reco"  ,"Reconstructed K;p, [GeV];",100,0,100); 
    h1_K_reco[cos_theta_K_reco] = new TH1F("cos_theta_K_reco"  ,"Reconstructed K;cos#theta;",200,-1,1); 

    h_PS[d0_P_single] = new TH1F("d0_P_single"  ,"Primary vertex single prong;d_0, [um];",2000,0,2);
    h_PS[d0_S_single] = new TH1F("d0_S_single"  ,"Secondary vertex single prong;d_0, [mm];",2000,0,2);
    h_PS[d0_P_mult] = new TH1F("d0_P_mult"  ,"Primary vertex mult prong;d_0, [um];",2000,0,2);
    h_PS[d0_S_mult] = new TH1F("d0_S_mult"  ,"Secondary vertex mult prong;d_0, [mm];",2000,0,2);

    h_PS[z0_P_single] = new TH1F("z0_P_single"  ,"Primary vertex single prong;z_0, [um];",2000,0,2);
    h_PS[z0_S_single] = new TH1F("z0_S_single"  ,"Secondary vertex single prong;z_0, [mm];",2000,0,2);
    h_PS[z0_P_mult] = new TH1F("z0_P_mult"  ,"Primary vertex mult prong;z_0, [um];",2000,0,2);
    h_PS[z0_S_mult] = new TH1F("z0_S_mult"  ,"Secondary vertex mult prong;z_0, [mm];",2000,0,2);

    h_tagging[p_ctag] = new TH1F("p_ctag","ctag for primary vertex",100,0,1.);
    h_tagging[s_ctag] = new TH1F("s_ctag","ctag for secondary vertex",100,0,1.);
    h_tagging[t_ctag] = new TH1F("t_ctag","ctag for no vertex",100,0,1.);
    
    h_tagging[p_btag] = new TH1F("p_btag","btag for primary vertex",100,0,1.);
    h_tagging[s_btag] = new TH1F("s_btag","btag for secondary vertex",100,0,1.);
    h_tagging[t_btag] = new TH1F("t_btag","btag for no vertex",100,0,1.);

    h_tagging[ctag] = new TH1F("ctag"  ,"Ctag values;ctag;",100,0,1);

    h_tagging[jets_info] = new TH1F("jets_info","Info", 10,0,10);
    h_tagging[jets_info]->GetXaxis()->SetBinLabel(1,"Events");
    h_tagging[jets_info]->GetXaxis()->SetBinLabel(2,"Jet 1");
    h_tagging[jets_info]->GetXaxis()->SetBinLabel(3,"Jet 2");
    h_tagging[jets_info]->GetXaxis()->SetBinLabel(4,"Jets with no vtx");
    h_tagging[jets_info]->GetXaxis()->SetBinLabel(5,"Jets with S");
    h_tagging[jets_info]->GetXaxis()->SetBinLabel(6,"Jets w/o S");
    h_tagging[jets_info]->GetXaxis()->SetBinLabel(7,"Events with both jets with S");
    h_tagging[jets_info]->GetXaxis()->SetBinLabel(8,"Events with both jets w/o S");

    h_general[n_K_ecal] = new TH1F("n_K_ecal","number of charged Kaons which reached calorimeter",5,-2,2);
    // h_general[nvtx_ctag] = new TH2F("nvtx_ctag",";ctag;nvtx",100,0,1,5,0,5);

    Hist2List();

}

void HistManager::Hist2List()
{

  for (auto ih : h1_K_reco) {
    hList_K_reco->Add(ih);
  }

  for (auto ih : h_general) {
    hList_general->Add(ih);
  }

  for (auto ih : h_PS) {
    hList_PS->Add(ih);
  }

  for (auto ih : h_tagging) {
    hList_tagging->Add(ih);
  }

}

void HistManager::WriteLists( TFile * output)
{
  output->cd(); 
  //Own histos
  hList_general->Write();
  
  TDirectory * d_reco = output->mkdir("Reconstructed");
    d_reco->cd();
    hList_K_reco->Write();
  TDirectory * d_PS = output->mkdir("PS");
    d_PS->cd();
    hList_PS->Write();

  TDirectory * d_tagging = output->mkdir("tagging");
    d_tagging->cd();
    hList_tagging->Write();
  
}