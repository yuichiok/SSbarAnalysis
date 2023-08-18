#include <iostream>
#include <vector>
#include "include/Styles.hh"
#include "include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector;

void BinNormal(TH1F *h)
{
  const Int_t nbins = h->GetNbinsX();
  for (int ibin=1; ibin<=nbins; ibin++){
    Float_t binc = h->GetBinContent(ibin);
    Float_t binw = h->GetBinWidth(ibin);
    Float_t bin_normal = binc / binw;
    h->SetBinContent(ibin,bin_normal);
  }

}

void Normalize(TH1F *h, Float_t norm_top)
{
  // h->Scale( 1.0 / h->GetEntries() );
  const Int_t nbins = h->GetNbinsX();
  Int_t nbins4 = nbins / 4;
  // Int_t int_high = (nbins / 2) + nbins4;
  // Int_t int_low  = (nbins / 2 + 1) - nbins4;
  // Int_t int_high = nbins-4;
  // Int_t int_low  = (nbins / 2);
  Int_t int_high = nbins-4;
  Int_t int_low  = (nbins / 2) + 4;
  h->Scale( norm_top / h->Integral(int_low,int_high) );
}

void NormalizeUU(TH1F *h, Float_t norm_top)
{
  const Int_t nbins = h->GetNbinsX();
  Int_t nbins4 = nbins / 4;
  Int_t int_high = nbins;
  Int_t int_low  = 1;

  cout << "uu integral = " << h->Integral(int_low,int_high) << endl;
  cout << "uu entries  = " << h->GetEntries() << endl;

  h->Scale( norm_top / h->Integral(int_low,int_high) );

}

void NormalizeGen(TH1F *h, Float_t norm_top)
{
  h->Scale( norm_top / h->GetEntries() );
}

void Normalize2Reco(TH1F *h_reco, TH1F *h_gen)
{
	double intCosReco = h_reco->Integral(10,90);
	double intCosGen  = h_gen->Integral(10,90);
  h_gen->Scale( intCosReco / intCosGen );
}

void StyleHist(TH1F *h, Color_t col)
{
  h->SetLineWidth(3);
  h->SetLineColor(col);
  h->SetFillStyle(3002);
  h->SetFillColor(col);
}

void StyleFunc(TF1 *f, Int_t styl, Color_t col)
{
  f->SetLineColor(col);
  f->SetLineWidth(5);
  f->SetLineStyle(styl);
}

void main_pq_BGFit( TFile *files[] )
{
  enum MixProcess {kUU,kDD,kUD};
  gStyle->SetOptStat(0);

  // gen uu/dd polar
  TH1F *h_gen_uu_qcos = (TH1F*) files[kUU]->Get("h_gen_q_qcos");
  TH1F *h_gen_dd_qcos = (TH1F*) files[kDD]->Get("h_gen_q_qcos");
  TH1F *h_gen_ud_qcos = (TH1F*) files[kUD]->Get("h_gen_q_qcos");

  StyleHist(h_gen_uu_qcos,kBlue+1);
  h_gen_uu_qcos->SetFillStyle(0);
  h_gen_uu_qcos->SetLineStyle(2);
  StyleHist(h_gen_dd_qcos,kGreen+2);
  h_gen_dd_qcos->SetFillStyle(0);
  h_gen_dd_qcos->SetLineStyle(2);
  StyleHist(h_gen_ud_qcos,kBlack);
  h_gen_ud_qcos->SetFillStyle(0);
  h_gen_ud_qcos->SetLineStyle(2);

  // reco us polar
  TH1F *h_reco_ud_Pi_scos  = (TH1F*) files[kUD]->Get("h_reco_Pi_scos");
  TH1F *h_reco_ud_Pi_qcos  = (TH1F*) files[kUD]->Get("h_reco_Pi_qcos");
  TH1F *h_cheat_ud_Pi_qcos = (TH1F*) files[kUD]->Get("h_cheat_Pi_qcos");

  // efficiency correction
  TH1F *h_reco_ud_Pi_scos_eff_corr  = resolutionCorrection(h_reco_ud_Pi_scos,"scos_corr",files[kUD]);
  TH1F *h_reco_ud_Pi_qcos_eff_corr  = resolutionCorrection(h_reco_ud_Pi_qcos,"qcos_corr",files[kUD]);
  TH1F *h_cheat_ud_Pi_qcos_eff_corr = resolutionCorrection(h_cheat_ud_Pi_qcos,"cheat_qcos_corr",files[kUD]);
  // TH1F *h_reco_ud_Pi_scos_eff_corr = (TH1F*) h_reco_ud_Pi_scos->Clone();
  // TH1F *h_reco_ud_Pi_qcos_eff_corr = (TH1F*) h_reco_ud_Pi_qcos->Clone();

  // used for pq correction
  TH1F *h_acc_PiPi_cos  = (TH1F*) files[kUD]->Get("pq/h_acc_PiPi_cos");
  TH1F *h_rej_PiPi_cos  = (TH1F*) files[kUD]->Get("pq/h_rej_PiPi_cos");

  TH1F *h_acc_PiPi_cos_eff_corr = resolutionCorrection(h_acc_PiPi_cos,"acc_corr",files[kUD]);
  TH1F *h_rej_PiPi_cos_eff_corr = resolutionCorrection(h_rej_PiPi_cos,"rej_corr",files[kUD]);
  // TH1F *h_acc_PiPi_cos_eff_corr = (TH1F*) h_acc_PiPi_cos->Clone();
  // TH1F *h_rej_PiPi_cos_eff_corr = (TH1F*) h_rej_PiPi_cos->Clone();

  StyleHist(h_reco_ud_Pi_scos_eff_corr,kBlue);
  h_reco_ud_Pi_scos_eff_corr->SetFillStyle(0);
  StyleHist(h_reco_ud_Pi_qcos_eff_corr,kRed+2);
  StyleHist(h_cheat_ud_Pi_qcos_eff_corr,kBlue);
  StyleHist(h_cheat_ud_Pi_qcos,kBlue);

  StyleHist(h_acc_PiPi_cos_eff_corr,kRed+2);
  StyleHist(h_rej_PiPi_cos_eff_corr,kBlue+2);

  const Int_t nbins = h_reco_ud_Pi_scos_eff_corr->GetNbinsX();

  TH1F *p_KK = new TH1F("p_KK", "p_KK", 50,0,1);
  p_KK->Sumw2();

  vector<Float_t> p_vec = GetP(h_acc_PiPi_cos_eff_corr, h_rej_PiPi_cos_eff_corr);

  for (unsigned i = 0; i < p_vec.size() / 2; i++)
  {
    p_KK->SetBinContent(nbins / 2 - i, p_vec.at(i));
    p_KK->SetBinError(nbins / 2 - i, p_vec.at(i + nbins / 2));
  }

  TH1F *h_reco_Pi_pq_cos = CorrectHist(h_reco_ud_Pi_qcos_eff_corr, p_vec);
  // TH1F *h_reco_Pi_pq_cos = (TH1F*) h_reco_ud_Pi_qcos_eff_corr->Clone();
  StyleHist(h_reco_Pi_pq_cos,kBlack);

  Normalize2Reco(h_reco_Pi_pq_cos,h_gen_ud_qcos);
  // Normalize2Reco(h_reco_Pi_pq_cos,h_cheat_ud_Pi_qcos);
  Normalize2Reco(h_gen_ud_qcos,h_cheat_ud_Pi_qcos);


  // Fit
  // Gen
  TF1 * f_gen_ud      = new TF1("f_gen_ud","[0]*(1+x*x)+[1]*x",-1.0,1.0);
  StyleFunc(f_gen_ud,2,kBlack);
  f_gen_ud->SetParNames("Su","Au");
  h_gen_ud_qcos->Fit("f_gen_ud","MNRS");

  // Reco
  TF1 * f_reco_ud      = new TF1("f_reco_ud","[0]*(1+x*x)+[1]*x",-0.8,0.8);
  StyleFunc(f_reco_ud,1,kRed);
  f_reco_ud->SetParNames("Su","Au");
  h_reco_Pi_pq_cos->Fit("f_reco_ud","MNRS");

  // Draw
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  TPad *pad1 = new TPad("pad1", "pad1",0,0,1,1);
  StylePad(pad1,0,0.12,0,0.12);

  h_reco_Pi_pq_cos->SetTitle(";cos#theta_{#pi};Entries");
  h_reco_Pi_pq_cos->GetYaxis()->SetRangeUser(0,30E3);

  h_reco_Pi_pq_cos->Draw("h");
  // h_reco_ud_Pi_scos_eff_corr->Draw("hsame");
  h_cheat_ud_Pi_qcos->Draw("hsame");
  // h_cheat_ud_Pi_qcos_eff_corr->Draw("hsame");
  h_gen_ud_qcos->Draw("hsame");
  f_gen_ud->Draw("same");
  f_reco_ud->Draw("same");

  cout << "Gen  Chi2 / ndf = " << f_gen_ud->GetChisquare() << " / " << f_gen_ud->GetNDF() << endl;
  cout << "Reco Chi2 / ndf = " << f_reco_ud->GetChisquare() << " / " << f_reco_ud->GetNDF() << endl;

  Float_t AFB_gen  = AFB_calculation(f_gen_ud);
  Float_t AFB_reco = AFB_calculation(f_reco_ud);

  cout << "Gen  AFB = " << AFB_gen << endl;
  cout << "Reco AFB = " << AFB_reco << endl;

  TLegend *leg = new TLegend(0.4,0.70,0.8,0.85);
  leg->SetLineColor(0);
  leg->AddEntry(h_reco_Pi_pq_cos,"Reco #pi angle","f");
  leg->AddEntry(h_gen_ud_qcos,"Gen #pi angle","f");
  leg->AddEntry(h_cheat_ud_Pi_qcos,"Cheated #pi angle","f");
  leg->AddEntry(f_reco_ud,"Reco #pi angle Fit","l");
  leg->AddEntry(f_gen_ud,"Gen #pi angle Fit","l");
  leg->Draw();

  // Plot p value
  TCanvas *c_pval = new TCanvas("c_pval","c_pval",800,800);
  TPad *pad_pval = new TPad("pad_pval", "pad_pval",0,0,1,1);
  StylePad(pad_pval,0,0.12,0,0.15);
  
  StyleHist(p_KK,kGreen+2);
  p_KK->SetTitle(";cos#theta_{#pi};p value");
  p_KK->GetYaxis()->SetRangeUser(0,1);
  p_KK->Draw("h");

}

void pq_method_PiBGFit_norm()
{
  TGaxis::SetMaxDigits(3);

  try
  {
    TString process[3] = {"uu","dd","ud"};
    TFile* files[3];
    for ( int i=0; i<3; i++ ){
      files[i] = new TFile(TString::Format("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.%s.KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root",process[i].Data()),"READ");
      if ( !files[i]->IsOpen() ) throw 0;
    }

    main_pq_BGFit( files );
  }
  catch ( int error_code ) {
    switch ( error_code ){
      default:
        cerr << "<< Error >>" << endl;
        break;
    }
  }

}