#include <iostream>

template <class P1>
void StylePad(P1 *pad, Float_t t, Float_t b, Float_t r, Float_t l)
{
  pad->SetGrid(1,1);
  if(t) pad->SetTopMargin(t);
  if(b) pad->SetBottomMargin(b);
  if(r) pad->SetRightMargin(r);
  if(l) pad->SetLeftMargin(l);
  pad->Draw();
  pad->cd();

}

template <class Hist>
void Normalize(Hist *h)
{
  h->Scale( 1.0 / h->GetEntries() );
}

template <class HistType>
void StyleHist(HistType *h, Color_t col)
{
  h->SetLineWidth(3);
  h->SetLineColor(col);
  h->SetFillStyle(3002);
  h->SetFillColor(col);
}

void main_Particle_Ratio(TFile *file_uu, TFile *file_ss)
{
  gStyle->SetOptStat(0);

  Int_t n_particles = 3;

  TCanvas *c_particle_ratio = new TCanvas("c_particle_ratio","c_particle_ratio",2400,600);
  c_particle_ratio->Divide(n_particles,1);

  TH1F *uu_h1_particle_ratio[n_particles];
  uu_h1_particle_ratio[0] = (TH1F*) file_uu->Get("particle_ratio/h_K_rate_gen");
  uu_h1_particle_ratio[1] = (TH1F*) file_uu->Get("particle_ratio/h_pi_rate_gen");
  uu_h1_particle_ratio[2] = (TH1F*) file_uu->Get("particle_ratio/h_p_rate_gen");

  StyleHist(uu_h1_particle_ratio[0],kRed+1);
  StyleHist(uu_h1_particle_ratio[1],kBlue+1);
  StyleHist(uu_h1_particle_ratio[2],kGreen+1);

  TH1F *ss_h1_particle_ratio[n_particles];
  ss_h1_particle_ratio[0] = (TH1F*) file_ss->Get("particle_ratio/h_K_rate_gen");
  ss_h1_particle_ratio[1] = (TH1F*) file_ss->Get("particle_ratio/h_pi_rate_gen");
  ss_h1_particle_ratio[2] = (TH1F*) file_ss->Get("particle_ratio/h_p_rate_gen");

  StyleHist(ss_h1_particle_ratio[0],kRed+1);
  StyleHist(ss_h1_particle_ratio[1],kBlue+1);
  StyleHist(ss_h1_particle_ratio[2],kGreen+1);

  for ( int i=0; i<n_particles; i++ ){

    c_particle_ratio->cd(i+1);
    StylePad(gPad,0,0.15,0,0.17);

    Normalize(uu_h1_particle_ratio[i]);
    Normalize(ss_h1_particle_ratio[i]);
    ss_h1_particle_ratio[i]->SetLineStyle(2);
    uu_h1_particle_ratio[i]->GetYaxis()->SetRangeUser(0,1);

    uu_h1_particle_ratio[i]->Draw("h");
    ss_h1_particle_ratio[i]->Draw("hsame");

    if(i==0){
      TLegend *leg = new TLegend(0.6,0.76,0.8,0.85);
      leg->SetLineColor(0);
      leg->AddEntry(uu_h1_particle_ratio[i],"u#bar{u}","l");
      leg->AddEntry(ss_h1_particle_ratio[i],"s#bar{s}","l");
      leg->Draw();
    }

  }

}

void main_Particle_Ratio_cos_gen(TFile *file_uu, TFile *file_ss)
{
  gStyle->SetOptStat(0);

  const Int_t n_particles = 3;
  TString particle_types[n_particles] = {"K","pi","p"};
  Color_t colors[n_particles] = {kRed+1,kBlue+1,kGreen+1};

  TH2F *uu_h2_particle_ratio_cos_gen[n_particles];

  for ( int i=0; i<n_particles; i++ ) {
    uu_h2_particle_ratio_cos_gen[i] = (TH2F*) file_uu->Get(TString::Format("particle_ratio/h2_%s_rate_cos_gen",particle_types[i].Data()));
    StyleHist(uu_h2_particle_ratio_cos_gen[i],colors[i]);
  }

  TH2F *ss_h2_particle_ratio_cos_gen[n_particles];

  for ( int i=0; i<n_particles; i++ ) {
    ss_h2_particle_ratio_cos_gen[i] = (TH2F*) file_ss->Get(TString::Format("particle_ratio/h2_%s_rate_cos_gen",particle_types[i].Data()));
    StyleHist(ss_h2_particle_ratio_cos_gen[i],colors[i]);
  }


  // Draw 2D
  TCanvas *c_2d_particle_ratio_cos = new TCanvas("c_2d_particle_ratio_cos","c_2d_particle_ratio_cos",2400,600);
  c_2d_particle_ratio_cos->Divide(n_particles,2);

  for ( int i=0; i<n_particles; i++ ){

    c_2d_particle_ratio_cos->cd(i+1);
    gPad->SetLogz();

    StylePad(gPad,0,0.15,0.17,0.17);

    Normalize(uu_h2_particle_ratio_cos_gen[i]);
    uu_h2_particle_ratio_cos_gen[i]->GetYaxis()->SetRangeUser(0.1,1.1);
    uu_h2_particle_ratio_cos_gen[i]->Draw("colz");

  }

  for ( int i=0; i<n_particles; i++ ){

    c_2d_particle_ratio_cos->cd(i+4);
    gPad->SetLogz();

    StylePad(gPad,0,0.15,0.17,0.17);

    Normalize(ss_h2_particle_ratio_cos_gen[i]);
    ss_h2_particle_ratio_cos_gen[i]->GetYaxis()->SetRangeUser(0.1,1.1);
    ss_h2_particle_ratio_cos_gen[i]->Draw("colz");

  }

  // Draw 1D Projection


/*
  for ( int i=0; i<n_particles; i++ ){

    c_particle_ratio_lowcos->cd(i+1);
    StylePad(gPad,0,0.15,0,0.17);

    Normalize(uu_h1_particle_ratio_cos[i]);
    Normalize(ss_h1_particle_ratio_cos[i]);

    uu_h1_particle_ratio_lowcos[i]->Draw("h");
    ss_h1_particle_ratio_lowcos[i]->Draw("hsame");

    if(i==0){
      TLegend *leg = new TLegend(0.6,0.76,0.8,0.85);
      leg->SetLineColor(0);
      leg->AddEntry(uu_h1_particle_ratio_lowcos[i],"u#bar{u}","l");
      leg->AddEntry(ss_h1_particle_ratio_lowcos[i],"s#bar{s}","l");
      leg->Draw();
    }

  }
*/

}

void Particle_Ratio()
{
  try{
    TFile *file_uu = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.LPFOp15_pNaN.tpc0.hists.all.root","READ");
    TFile *file_ss = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.LPFOp15_pNaN.tpc0.hists.all.root","READ");
    main_Particle_Ratio(file_uu,file_ss);
    main_Particle_Ratio_cos_gen(file_uu,file_ss);
  }
  catch (int error_code) {
    switch ( error_code ){
      default:
        break;
    }
  }

}