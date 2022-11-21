#include <iostream>

void StylePad(TPad *pad, Float_t t, Float_t b, Float_t r, Float_t l)
{
  pad->SetGrid(1,1);
  if(t) pad->SetTopMargin(t);
  if(b) pad->SetBottomMargin(b);
  if(r) pad->SetRightMargin(r);
  if(l) pad->SetLeftMargin(l);
  pad->Draw();
  pad->cd();

}

double BetheBloch(const double *x, const double *pars){
  double mass=pars[5];
  double bg=x[0]/mass;
  double b=sqrt(bg*bg/(1.0+bg*bg));
  double tmax=pars[2]*TMath::Power(bg,2.0);   ///(1.0+pars[3]*g+pars[4]);

  return 1.350e-1*(0.5*pars[0]*TMath::Log(pars[1]*TMath::Power(bg,2.0)*tmax)-pars[3]*b*b-pars[4]*bg/2.0)/(b*b);
}

void tree2hist()
{

  Float_t bins_p[]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,36,40,44,48,52,56,60,64,68,72,80,90,100};
  Int_t   nbins_p=sizeof(bins_p)/sizeof(Float_t) - 1;

  Float_t bins_dEdx[200];
  bins_dEdx[0]=0.1;
  for(int i=1;i<200;i++) bins_dEdx[i]=bins_dEdx[i-1]+0.1/100.;
  Int_t nbins_dEdx=199;

  gStyle->SetOptStat(0);

  TFile *file    = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.LPFOp10_pNaN.tpc0.hists.all.root","READ");

  TTree* tree = (TTree*) file->Get("data");

  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  TPad *pad0 = new TPad("pad0", "pad0",0,0,1,1);
  StylePad(pad0,0,0.15,0,0.17);

  TH2F *h2_gen_K_dEdx_p  = new TH2F("h2_gen_K_dEdx_p"  ,";dE/dx;p (GeV)",nbins_p,bins_p,nbins_dEdx,bins_dEdx);
  TH2F *h2_gen_pi_dEdx_p = new TH2F("h2_gen_pi_dEdx_p" ,";dE/dx;p (GeV)",nbins_p,bins_p,nbins_dEdx,bins_dEdx);
  TH2F *h2_gen_p_dEdx_p  = new TH2F("h2_gen_p_dEdx_p"  ,";dE/dx;p (GeV)",nbins_p,bins_p,nbins_dEdx,bins_dEdx);

  Int_t n_gen_K_dEdx_p  = tree->Draw("vpfo_dedx:vpfo_cos >> h2_gen_K_dEdx_p", "abs(vpfo_pdgcheat)==321");
  Int_t n_gen_pi_dEdx_p = tree->Draw("vpfo_dedx:vpfo_cos >> h2_gen_pi_dEdx_p", "abs(vpfo_pdgcheat)==211");
  Int_t n_gen_p_dEdx_p = tree->Draw("vpfo_dedx:vpfo_cos >> h2_gen_p_dEdx_p", "abs(vpfo_pdgcheat)==2212");

  // TH2F *h2_gen_K_dEdx_p  = (TH2F*) file->Get("dEdx/h2_gen_K_dEdx_p");
  // TH2F *h2_gen_pi_dEdx_p = (TH2F*) file->Get("dEdx/h2_gen_pi_dEdx_p");
  // TH2F *h2_gen_p_dEdx_p  = (TH2F*) file->Get("dEdx/h2_gen_p_dEdx_p");

  h2_gen_K_dEdx_p->SetMarkerColor(kRed);
  h2_gen_pi_dEdx_p->SetMarkerColor(kBlue);
  h2_gen_p_dEdx_p->SetMarkerColor(kGreen);

  h2_gen_K_dEdx_p->SetLineColor(kRed);
  h2_gen_pi_dEdx_p->SetLineColor(kBlue);
  h2_gen_p_dEdx_p->SetLineColor(kGreen);

  h2_gen_K_dEdx_p->SetFillColor(kRed);
  h2_gen_pi_dEdx_p->SetFillColor(kBlue);
  h2_gen_p_dEdx_p->SetFillColor(kGreen);

  pad0->SetLogx();
  h2_gen_K_dEdx_p->SetTitle(";Track momentum [GeV];#frac{dE}{dx}[MeV]");
  h2_gen_K_dEdx_p->GetXaxis()->SetTitleOffset(1.5);
  h2_gen_K_dEdx_p->Draw("box");
  h2_gen_pi_dEdx_p->Draw("box same");
  h2_gen_p_dEdx_p->Draw("box same");

  // Kaon Bethe-Bloch formula
  std::vector< double > parskaon;
  parskaon.push_back(0.0792784);
  parskaon.push_back(3798.12);
  parskaon.push_back(4.06952e+07);
  parskaon.push_back(0.450671);
  parskaon.push_back(0.00050169);
  parskaon.push_back(0.493677); // mass

  TF1 *func = new TF1("func",BetheBloch,0.1,100,6) ;
  for (int i = 0; i < 6; ++i) func->SetParameter(i,parskaon[i]);
  func->SetLineColor(kBlack);
  func->Draw("same");

  TLegend *leg = new TLegend(0.5,0.76,0.75,0.85);
  leg->SetLineColor(0);
  leg->AddEntry(h2_gen_K_dEdx_p,"K^{#pm}","l");
  leg->AddEntry(h2_gen_pi_dEdx_p,"#pi^{#pm}","l");
  leg->AddEntry(h2_gen_p_dEdx_p,"p","l");
  leg->AddEntry(func,"K Bethe-Bloch formula","l");
  leg->Draw();

  TFile *outfile = new TFile("../rootfiles/plots/plot_dEdx.root");
  c0->Write();
  outfile->Close();

}