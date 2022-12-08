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

void Normalize(TH1F *h)
{
  // h->Scale( 1.0 / h->GetEntries() );
  h->Scale( 1.0 / h->Integral(30,70) );
  // h->Scale( 1.0 / h->Integral(12,88) );
}

void StyleHist(TH1F *h, Color_t col)
{
  h->SetLineWidth(3);
  h->SetLineColor(col);
  h->SetFillStyle(3002);
  h->SetFillColor(col);
}

void main_N_Particles(TFile *file)
{
  gStyle->SetOptStat(0);

  Int_t n_particles = 3;

  TCanvas *c_n_particles = new TCanvas("c_n_particles","c_n_particles",2400,600);
  c_n_particles->Divide(n_particles,1);
  TPad *pad_n_particles[n_particles];

  TH2F *h2_n_particles[n_particles];
  h2_n_particles[0] = (TH2F*) file->Get("h2_nK_gen_reco");
  h2_n_particles[1] = (TH2F*) file->Get("h2_npi_gen_reco");
  h2_n_particles[2] = (TH2F*) file->Get("h2_np_gen_reco");

  for ( int i=0; i<n_particles; i++ ){

    c_n_particles->cd(i+1);
    StylePad(gPad,0,0.15,0.17,0.17);

    h2_n_particles[i]->Draw("colz");

  }



}

void N_Particles()
{
  try{
    TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.LPFOp15_pNaN.tpc0.hists.all.root","READ");
    main_N_Particles(file);
  }
  catch (int error_code) {
    switch ( error_code ){
      default:
        break;
    }
  }

}