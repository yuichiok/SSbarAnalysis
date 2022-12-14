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

  Float_t bins_p[]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,64,68,72,80,90,100};
  Int_t   nbins_p=sizeof(bins_p)/sizeof(Float_t) - 1;

  Float_t bins_dEdx[200];
  bins_dEdx[0]=0.1;
  for(int i=1;i<200;i++) bins_dEdx[i]=bins_dEdx[i-1]+0.1/100.;
  Int_t nbins_dEdx=199;

  //////////////////
  //     TH1F     //
  //////////////////
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

    h1[gen_N_K_cos2]    = new TH1F("h_gen_N_K_cos2",";cos#theta;Generated N Kaons", 100,-1,1);
    h1[reco_N_K_cos2]   = new TH1F("h_reco_N_K_cos2",";cos#theta;Reconstructed N Kaons", 100,-1,1);
    h1[N_K_corr_cos2]   = new TH1F("h_N_K_corr_cos2",";cos#theta;Correctly Reconstructed N Kaons", 100,-1,1);

  // pq method
    h1_pq[acc_KK]      = new TH1F("h_acc_KK_cos",";Accepted K^{+}K^{-} cos#theta;N_{acc}",nbins_cos,bins_cos_fine);
    h1_pq[rej_KK]      = new TH1F("h_rej_KK_cos",";Rejected K^{+}K^{-} cos#theta;N_{rej}",nbins_cos,bins_cos_fine);

  // particle ratio
    h1_particle_ratio[K_rate_gen]  = new TH1F("h_K_rate_gen",";Ratio of Kaons / Event (gen);Entries",11,0,1.1);
    h1_particle_ratio[pi_rate_gen] = new TH1F("h_pi_rate_gen",";Ratio of Pions / Event (gen);Entries",11,0,1.1);
    h1_particle_ratio[p_rate_gen]  = new TH1F("h_p_rate_gen",";Ratio of Protons / Event (gen);Entries",11,0,1.1);

    h1_particle_ratio[K_rate_reco]  = new TH1F("h_K_rate_reco",";Ratio of Kaons / Event (reco);Entries",11,0,1.1);
    h1_particle_ratio[pi_rate_reco] = new TH1F("h_pi_rate_reco",";Ratio of Pions / Event (reco);Entries",11,0,1.1);
    h1_particle_ratio[p_rate_reco]  = new TH1F("h_p_rate_reco",";Ratio of Protons / Event (reco);Entries",11,0,1.1);


  //////////////////
  //     TH2F     //
  //////////////////
    h2[gen_K_p_cos]  = new TH2F("h2_gen_K_p_cos" ,";cos#theta;p (GeV)",cos_bin,-1,1,100,0,60);
    h2[reco_K_p_cos] = new TH2F("h2_reco_K_p_cos",";cos#theta;p (GeV)",cos_bin,-1,1,100,0,60);

    h2[nK_gen_reco]  = new TH2F("h2_nK_gen_reco" ,";N Kaons (Reco);N Kaons (Gen)",   50,0,50,50,0,50);
    h2[npi_gen_reco] = new TH2F("h2_npi_gen_reco",";N Pions (Reco);N Pions (Gen)",   50,0,50,50,0,50);
    h2[np_gen_reco]  = new TH2F("h2_np_gen_reco",";N Protons (Reco);N Protons (Gen)",50,0,50,50,0,50);
    
    h2[stable_cos]   = new TH2F("h2_stable_cos",";cos#theta;Stability",cos_bin,-1,1,50,0,1);
    h2[purity_cos]   = new TH2F("h2_purity_cos",";cos#theta;Purity",   cos_bin,-1,1,50,0,1);

  // particle ratio
    h2_particle_ratio_cos[K_rate_cos_gen]   = new TH2F("h2_K_rate_cos_gen",";qcos#theta;Ratio of Kaons / Event (gen)",40,-1,1,11,0,1.1);
    h2_particle_ratio_cos[pi_rate_cos_gen]  = new TH2F("h2_pi_rate_cos_gen",";qcos#theta;Ratio of Pions / Event (gen)",40,-1,1,11,0,1.1);
    h2_particle_ratio_cos[p_rate_cos_gen]   = new TH2F("h2_p_rate_cos_gen",";qcos#theta;Ratio of Protons / Event (gen)",40,-1,1,11,0,1.1);

    h2_particle_ratio_cos[K_rate_cos_reco]  = new TH2F("h2_K_rate_cos_reco",";qcos#theta;Ratio of Kaons / Event (reco)",40,-1,1,11,0,1.1);
    h2_particle_ratio_cos[pi_rate_cos_reco] = new TH2F("h2_pi_rate_cos_reco",";qcos#theta;Ratio of Pions / Event (reco)",40,-1,1,11,0,1.1);
    h2_particle_ratio_cos[p_rate_cos_reco]  = new TH2F("h2_p_rate_cos_reco",";qcos#theta;Ratio of Protons / Event (reco)",40,-1,1,11,0,1.1);

  // dEdx
    h2_dEdx[gen_K_dEdx_p]  = new TH2F("h2_gen_K_dEdx_p"  ,";p (GeV);#frac{dE}{dx}",nbins_p,bins_p,nbins_dEdx,bins_dEdx);
    h2_dEdx[gen_pi_dEdx_p] = new TH2F("h2_gen_pi_dEdx_p" ,";p (GeV);#frac{dE}{dx}",nbins_p,bins_p,nbins_dEdx,bins_dEdx);
    h2_dEdx[gen_p_dEdx_p]  = new TH2F("h2_gen_p_dEdx_p"  ,";p (GeV);#frac{dE}{dx}",nbins_p,bins_p,nbins_dEdx,bins_dEdx);

    h2_dEdx[gen_K_KdEdx_dist_cos]  = new TH2F("h2_gen_K_KdEdx_dist_cos"  ,";cos#theta;K #frac{dE}{dx} dist", 100, -1, 1, 100, -20, 20);
    h2_dEdx[gen_pi_KdEdx_dist_cos] = new TH2F("h2_gen_pi_KdEdx_dist_cos" ,";cos#theta;K #frac{dE}{dx} dist", 100, -1, 1, 100, -20, 20);
    h2_dEdx[gen_p_KdEdx_dist_cos]  = new TH2F("h2_gen_p_KdEdx_dist_cos"  ,";cos#theta;K #frac{dE}{dx} dist", 100, -1, 1, 100, -20, 20);

    h2_dEdx[reco_K_KdEdx_dist_cos] = new TH2F("h2_reco_K_KdEdx_dist_cos" ,";cos#theta;K #frac{dE}{dx} dist", 100, -1, 1, 100, -20, 20);

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

  for (auto ih : h1_particle_ratio) {
    hList1_particle_ratio->Add(ih);
  }

  for (auto ih : h2) {
    hList2->Add(ih);
  }
  
  for (auto ih : h2_particle_ratio_cos) {
    hList1_particle_ratio->Add(ih);
  }

  for (auto ih : h2_dEdx) {
    hList2_dEdx->Add(ih);
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

  TDirectory * d_particle_ratio = output->mkdir("particle_ratio");
    d_particle_ratio->cd();
    hList1_particle_ratio->Write();
    output->cd();

  TDirectory * d_dEdx = output->mkdir("dEdx");
    d_dEdx->cd();
    hList2_dEdx->Write();
    output->cd();

}