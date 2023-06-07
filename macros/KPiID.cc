#include <iostream>

void Normalize(TH1F *h)
{
  h->Scale( 1.0 / h->GetEntries() );
}

void StyleHist(TH1F *h, Color_t col)
{
  h->SetLineWidth(3);
  h->SetLineColor(col);
  h->SetFillStyle(3002);
  h->SetFillColor(col);
}

void SetKPiIDLabels(TH1F *h)
{
  h->GetXaxis()->SetBinLabel(1,"Pion");
  h->GetXaxis()->SetBinLabel(2,"Kaon");
  h->GetXaxis()->SetBinLabel(3,"Proton");
  h->GetXaxis()->SetBinLabel(4,"Electron");
  h->GetXaxis()->SetBinLabel(5,"Muon");
  h->GetXaxis()->SetBinLabel(6,"Others");
}

void KPiID()
{
  gStyle->SetOptStat(0);

  enum process {kUU, kDD, kSS};
  Color_t col[3] = {kBlue+2, kBlack, kRed+2};
  TString process_name[3] = {"uu","dd","ss"};
  TFile *files[3];

  TH1F *h_reco_K_pdgcheat[3];
  TH1F *h_reco_Pi_pdgcheat[3];

  TH1F * h_gen_reco_K_sep_cos[3];
  TH1F * h_gen_reco_Pi_sep_cos[3];

  for( int i=0; i<3; i++ )
  {
    TString str = "../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR." + process_name[i] + ".KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root";
    files[i] = new TFile(str.Data(),"READ");

    h_reco_K_pdgcheat[i]  = (TH1F*) files[i]->Get("h_reco_K_pdgcheat");
    h_reco_Pi_pdgcheat[i] = (TH1F*) files[i]->Get("h_reco_Pi_pdgcheat");

    h_gen_reco_K_sep_cos[i]  = (TH1F*) files[i]->Get("h_gen_reco_K_sep_cos");
    h_gen_reco_Pi_sep_cos[i] = (TH1F*) files[i]->Get("h_gen_reco_Pi_sep_cos");

    cout << process_name[i] << " K:  " << h_reco_K_pdgcheat[i]->GetBinContent(2) / h_reco_K_pdgcheat[i]->GetEntries() * 100 << "% ";
                       cout << "Pi: " << h_reco_Pi_pdgcheat[i]->GetBinContent(1) / h_reco_Pi_pdgcheat[i]->GetEntries() * 100 << "%" << endl;

    StyleHist(h_reco_K_pdgcheat[i],col[i]);
    StyleHist(h_reco_Pi_pdgcheat[i],col[i]);
    StyleHist(h_gen_reco_K_sep_cos[i],col[i]);
    StyleHist(h_gen_reco_Pi_sep_cos[i],col[i]);

    Normalize(h_reco_K_pdgcheat[i]);
    Normalize(h_reco_Pi_pdgcheat[i]);
    Normalize(h_gen_reco_K_sep_cos[i]);
    Normalize(h_gen_reco_Pi_sep_cos[i]);
  }

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  gPad->SetGrid(1,1);

  TLegend *leg = new TLegend(0.70,0.75,0.85,0.85);
  leg->SetMargin(0.7);
  for (int i=0; i<3; i++)
  {
    if(i==0){
      h_reco_K_pdgcheat[i]->SetTitle("Truth PID of Reconstructed Leading Kaon;;Ratio");
      h_reco_K_pdgcheat[i]->GetYaxis()->SetRangeUser(0,1);
      SetKPiIDLabels(h_reco_K_pdgcheat[i]);
      h_reco_K_pdgcheat[i]->Draw("h");
    }else{
      h_reco_K_pdgcheat[i]->Draw("hsame");
    }

    leg->SetLineColor(0);
    leg->AddEntry(h_reco_K_pdgcheat[i],process_name[i],"l");
  }
  
  c0->Draw();
  leg->Draw();


  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  gPad->SetGrid(1,1);

  TLegend *leg1 = new TLegend(0.70,0.75,0.85,0.85);
  leg->SetMargin(0.7);
  for (int i=0; i<3; i++)
  {
    if(i==0){
      h_reco_Pi_pdgcheat[i]->SetTitle("Truth PID of Reconstructed Leading Pion;;Ratio");
      h_reco_Pi_pdgcheat[i]->GetYaxis()->SetRangeUser(0,1);
      SetKPiIDLabels(h_reco_Pi_pdgcheat[i]);
      h_reco_Pi_pdgcheat[i]->Draw("h");
    }else{
      h_reco_Pi_pdgcheat[i]->Draw("hsame");
    }

    leg1->SetLineColor(0);
    leg1->AddEntry(h_reco_Pi_pdgcheat[i],process_name[i],"l");
  }
  
  c1->Draw();
  leg1->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  gPad->SetGrid(1,1);

  TLegend *leg2 = new TLegend(0.15,0.75,0.45,0.85);
  for (int i=0; i<3; i++)
  {
    if(i==0){
      h_gen_reco_K_sep_cos[i]->GetYaxis()->SetRangeUser(0,1);
      h_gen_reco_K_sep_cos[i]->Draw("h");
    }else{
      h_gen_reco_K_sep_cos[i]->Draw("hsame");
    }

    leg2->SetLineColor(0);
    leg2->AddEntry(h_gen_reco_K_sep_cos[i],process_name[i],"l");
  }
  
  c2->Draw();
  leg2->Draw();

  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  gPad->SetGrid(1,1);

  TLegend *leg3 = new TLegend(0.15,0.75,0.45,0.85);
  for (int i=0; i<3; i++)
  {
    if(i==0){
      h_gen_reco_Pi_sep_cos[i]->GetYaxis()->SetRangeUser(0,1);
      h_gen_reco_Pi_sep_cos[i]->Draw("h");
    }else{
      h_gen_reco_Pi_sep_cos[i]->Draw("hsame");
    }

    leg3->SetLineColor(0);
    leg3->AddEntry(h_gen_reco_Pi_sep_cos[i],process_name[i],"l");
  }
  
  c3->Draw();
  leg3->Draw();

}