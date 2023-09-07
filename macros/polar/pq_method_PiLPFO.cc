#include <iostream>
#include <vector>

#include "../include/Styles.hh"
#include "../include/PolarTools.hh"

using std::cout; using std::endl;
using std::vector; using std::unordered_map;

TString prod_mode = "uu";
TString chiral     = "eR.pL";
TString LPFO_mode = "Pi";

TFile *file = new TFile("../../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h." + chiral + "." + prod_mode + ".KPiLPFO.PFOp15.LPFOp15_pNaN.tpc0.mix_uds.correctDist.all.root","READ");
enum correctionType {typeEfficiency, typeResolution, typeBoth, typeNone};

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

void Normalize(TH1F *h)
{
  // h->Scale( 1.0 / h->GetEntries() );
  const Int_t nbins = h->GetNbinsX();
  Int_t nbins4 = nbins / 4;
  Int_t int_high = (nbins / 2) + nbins4;
  Int_t int_low  = (nbins / 2 + 1) - nbins4;
  h->Scale( 1.0 / h->Integral(int_low,int_high) );
}

void Normalize2Gen(TH1F *h, TH1F *h_gen)
{
	double intCosReco = h->Integral(20,30);
	double intCosGen  = h_gen->Integral(20,30);
  h->Scale( intCosGen / intCosReco );
}


void main_pq()
{
  gStyle->SetOptStat(0);

  TString prodMode = getProductionMode(file);

  // reco and gen polar
  TH1F *h_gen_q_qcos    = (TH1F*) file->Get("gen/h_qcos");
  TH1F *h_reco_Pi_scos  = (TH1F*) file->Get("cos/h_Pi_scos");
  TH1F *h_reco_Pi_qcos  = (TH1F*) file->Get("cos/h_Pi_qcos");

  // used for pq correction
  TH1F *h_acc_PiPi_cos  = (TH1F*) file->Get("cos/h_Pi_acc_cos");
  TH1F *h_rej_PiPi_cos  = (TH1F*) file->Get("cos/h_Pi_rej_cos");

  // efficiency correction
  correctionType corrType = typeBoth;
  
  TH1F *h_reco_Pi_scos_eff_corr;
  TH1F *h_reco_Pi_qcos_eff_corr;

  TH1F *h_acc_PiPi_cos_eff_corr;
  TH1F *h_rej_PiPi_cos_eff_corr;

  switch (corrType)
  {
    default:
      break;

    case typeNone:
      h_reco_Pi_scos_eff_corr = (TH1F*) h_reco_Pi_scos->Clone();
      h_reco_Pi_qcos_eff_corr = (TH1F*) h_reco_Pi_qcos->Clone();
      h_acc_PiPi_cos_eff_corr = (TH1F*) h_acc_PiPi_cos->Clone();
      h_rej_PiPi_cos_eff_corr = (TH1F*) h_rej_PiPi_cos->Clone();
      break;    

    case typeEfficiency:
      h_reco_Pi_scos_eff_corr = efficiencyCorrection(h_reco_Pi_scos,"Pi",file);
      h_reco_Pi_qcos_eff_corr = efficiencyCorrection(h_reco_Pi_qcos,"Pi",file);
      h_acc_PiPi_cos_eff_corr = efficiencyCorrection(h_acc_PiPi_cos,"Pi",file);
      h_rej_PiPi_cos_eff_corr = efficiencyCorrection(h_rej_PiPi_cos,"Pi",file);
      break;
    
    case typeResolution:
      h_reco_Pi_scos_eff_corr = resolutionCorrection(h_reco_Pi_scos,"Pi",file);
      h_reco_Pi_qcos_eff_corr = resolutionCorrection(h_reco_Pi_qcos,"Pi",file);
      h_acc_PiPi_cos_eff_corr = resolutionCorrection(h_acc_PiPi_cos,"Pi",file);
      h_rej_PiPi_cos_eff_corr = resolutionCorrection(h_rej_PiPi_cos,"Pi",file);
      break;
    
    case typeBoth:
      TH1F* h_reco_Pi_scos_eff = efficiencyCorrection(h_reco_Pi_scos,"Pi",file);
      h_reco_Pi_scos_eff_corr = resolutionCorrection(h_reco_Pi_scos_eff,"Pi",file);

      TH1F* h_reco_Pi_qcos_eff = efficiencyCorrection(h_reco_Pi_qcos,"Pi",file);
      h_reco_Pi_qcos_eff_corr = resolutionCorrection(h_reco_Pi_qcos_eff,"Pi",file);

      TH1F* h_acc_PiPi_cos_eff = efficiencyCorrection(h_acc_PiPi_cos,"Pi",file);
      h_acc_PiPi_cos_eff_corr = resolutionCorrection(h_acc_PiPi_cos_eff,"Pi",file);

      TH1F* h_rej_PiPi_cos_eff = efficiencyCorrection(h_rej_PiPi_cos,"Pi",file);
      h_rej_PiPi_cos_eff_corr = resolutionCorrection(h_rej_PiPi_cos_eff,"Pi",file);
      
      break;

  }

  // if (isEffCorr)
  // {
  //   // scos & qcos
  //   h_reco_Pi_scos_eff_corr = efficiencyCorrection(h_reco_Pi_scos,"Pi",file);
  //   h_reco_Pi_qcos_eff_corr = efficiencyCorrection(h_reco_Pi_qcos,"Pi",file);
  //   // h_reco_Pi_scos_eff_corr = resolutionCorrection(h_reco_Pi_scos,"Pi",file);
  //   // h_reco_Pi_qcos_eff_corr = resolutionCorrection(h_reco_Pi_qcos,"Pi",file);

  //   // TH1F* h_reco_Pi_scos_eff = efficiencyCorrection(h_reco_Pi_scos,"Pi",file);
  //   // TH1F* h_reco_Pi_qcos_eff = efficiencyCorrection(h_reco_Pi_qcos,"Pi",file);
  //   // h_reco_Pi_scos_eff_corr = resolutionCorrection(h_reco_Pi_scos_eff,"Pi",file);
  //   // h_reco_Pi_qcos_eff_corr = resolutionCorrection(h_reco_Pi_qcos_eff,"Pi",file);

  //   // acc & rej
  //   h_acc_PiPi_cos_eff_corr = efficiencyCorrection(h_acc_PiPi_cos,"Pi",file);
  //   h_rej_PiPi_cos_eff_corr = efficiencyCorrection(h_rej_PiPi_cos,"Pi",file);
  //   // h_acc_PiPi_cos_eff_corr = resolutionCorrection(h_acc_PiPi_cos,"Pi",file);
  //   // h_rej_PiPi_cos_eff_corr = resolutionCorrection(h_rej_PiPi_cos,"Pi",file);

  //   // TH1F* h_acc_PiPi_cos_eff = efficiencyCorrection(h_acc_PiPi_cos,"Pi",file);
  //   // TH1F* h_rej_PiPi_cos_eff = efficiencyCorrection(h_rej_PiPi_cos,"Pi",file);
  //   // h_acc_PiPi_cos_eff_corr = resolutionCorrection(h_acc_PiPi_cos_eff,"Pi",file);
  //   // h_rej_PiPi_cos_eff_corr = resolutionCorrection(h_rej_PiPi_cos_eff,"Pi",file);

  // }else{
  //   h_reco_Pi_scos_eff_corr = (TH1F*) h_reco_Pi_scos->Clone();
  //   h_reco_Pi_qcos_eff_corr = (TH1F*) h_reco_Pi_qcos->Clone();

  //   h_acc_PiPi_cos_eff_corr = (TH1F*) h_acc_PiPi_cos->Clone();
  //   h_rej_PiPi_cos_eff_corr = (TH1F*) h_rej_PiPi_cos->Clone();
  // }

  StyleHist(h_gen_q_qcos,kGreen+1);
  h_gen_q_qcos->SetFillStyle(0);
  h_gen_q_qcos->SetLineStyle(2);

  StyleHist(h_reco_Pi_scos_eff_corr,kBlack);
  h_reco_Pi_scos_eff_corr->SetFillStyle(0);
  StyleHist(h_reco_Pi_qcos_eff_corr,kRed+2);
  StyleHist(h_acc_PiPi_cos_eff_corr,kRed+2);
  StyleHist(h_rej_PiPi_cos_eff_corr,kBlue+2);

  const Int_t nbins = h_reco_Pi_scos_eff_corr->GetNbinsX();

  TH1F *p_PiPi = new TH1F("p_PiPi", "p_PiPi", 50,0,1);
  p_PiPi->Sumw2();

  vector<Float_t> p_vec = GetP(h_acc_PiPi_cos_eff_corr, h_rej_PiPi_cos_eff_corr);

  for (unsigned i = 0; i < p_vec.size() / 2; i++)
  {
    p_PiPi->SetBinContent(nbins / 2 - i, p_vec.at(i));
    p_PiPi->SetBinError(nbins / 2 - i, p_vec.at(i + nbins / 2));
  }

  TH1F *h_reco_Pi_pq_cos = CorrectHist(prodMode, h_reco_Pi_qcos_eff_corr, p_vec);
  // TH1F *h_reco_Pi_pq_cos = (TH1F*) h_reco_Pi_qcos_eff_corr->Clone();
  StyleHist(h_reco_Pi_pq_cos,kBlue);

  Normalize2Gen(h_gen_q_qcos,h_reco_Pi_scos_eff_corr);
  // Normalize2Gen(h_gen_q_qcos,h_reco_Pi_pq_cos);

  // Fitting
  TF1 * f_gen = new TF1("f_gen","[0]*(1+x*x)+[1]*x",-0.8,0.8);
  f_gen->SetParNames("S","A");
  h_gen_q_qcos->Fit("f_gen","MNRS");
  cout << "Gen Chi2 / ndf = " << f_gen->GetChisquare() << " / " << f_gen->GetNDF() << endl;

  TF1 * f_reco = new TF1("f_reco","[0]*(1+x*x)+[1]*x",-0.8,0.8);
  f_reco->SetParNames("S","A");
  h_reco_Pi_pq_cos->Fit("f_reco","MNRS");
  cout << "Reco Chi2 / ndf = " << f_reco->GetChisquare() << " / " << f_reco->GetNDF() << endl;
  
  // output AFB
  Float_t AFB_gen  = AFB_calculation(f_gen);
  Float_t AFB_reco = AFB_calculation(f_reco);
  cout << "Gen  AFB = " << AFB_gen << endl;
  cout << "Reco AFB = " << AFB_reco << endl;


  // Draw polar angle
  TCanvas *c_polar = new TCanvas("c_polar","c_polar",800,800);
  TPad  *pad_polar = new TPad("pad_polar", "pad_polar",0,0,1,1);
  StylePad(pad_polar,0,0.12,0,0.15);

  h_reco_Pi_pq_cos->SetTitle(";cos#theta_{#pi^{-}};a.u.");
  h_reco_Pi_pq_cos->Draw("h");
  h_reco_Pi_qcos_eff_corr->Draw("hsame");
  h_reco_Pi_scos_eff_corr->Draw("hsame");
  h_gen_q_qcos->Draw("hsame");

  f_reco->Draw("same");

  TLegend *leg = new TLegend(0.3,0.7,0.75,0.85);
  leg->SetMargin(0.4);
  leg->SetLineColor(0);
  leg->AddEntry(h_gen_q_qcos,"Generated quark angle","l");
  leg->AddEntry(h_reco_Pi_scos_eff_corr,"Reconstructed #pi^{-} matched with quark angle","l");
  leg->AddEntry(h_reco_Pi_qcos_eff_corr,"Reconstructed #pi^{-}","l");
  leg->AddEntry(h_reco_Pi_pq_cos,"Reconstructed #pi^{-} (corrected)","l");
  leg->Draw();

  // Draw polar angle fit ratio
  TCanvas  *c_ratio = new TCanvas("c_ratio","c_ratio",800,800);
  TH1F *h_reco_Pi_pq_cos_subhist = new TH1F("h_reco_Pi_pq_cos_subhist",";LPFO Pion cos#theta; Entries",80,-0.8,0.8);
  for ( int ibin=1; ibin<=h_reco_Pi_pq_cos_subhist->GetNbinsX(); ibin++ )
  {
    Int_t recobin = ibin + 10;
    h_reco_Pi_pq_cos_subhist->SetBinContent(ibin,h_reco_Pi_pq_cos->GetBinContent(recobin));
    h_reco_Pi_pq_cos_subhist->SetBinError(ibin,h_reco_Pi_pq_cos->GetBinError(recobin));
  }
  h_reco_Pi_pq_cos_subhist->SetMarkerStyle(20);

  TF1 *f_reco_ratio = new TF1("f_reco_ratio","[0]*(1+x*x)+[1]*x",-0.8,0.8);
  f_reco_ratio->SetParNames("S","A");
  h_reco_Pi_pq_cos_subhist->Fit("f_reco_ratio");
  c_ratio->Clear();

  auto trp = new TRatioPlot(h_reco_Pi_pq_cos_subhist);
  trp->SetGraphDrawOpt("P");
  trp->SetSeparationMargin(0.0);
  trp->Draw();
  trp->GetLowerRefYaxis()->SetTitle("Ratio");
  trp->GetLowerRefGraph()->SetMinimum(-4);
  trp->GetLowerRefGraph()->SetMaximum(4);

  trp->GetUpperPad()->cd();
  TLegend *leg_trp = new TLegend(0.3,0.7,0.7,0.85);
  leg_trp->SetMargin(0.4);
  leg_trp->SetLineColor(0);
  leg_trp->AddEntry(h_reco_Pi_pq_cos_subhist,"Reconstructed #pi^{-} (corrected)","p");
  leg_trp->AddEntry(f_reco_ratio,"#frac{d#sigma}{dcos#theta} fit for LPFO","l");
  leg_trp->Draw();

  // Draw p value
  TCanvas *c_pval = new TCanvas("c_pval","c_pval",800,800);
  TPad  *pad_pval = new TPad("pad_pval", "pad_pval",0,0,1,1);
  StylePad(pad_pval,0,0.12,0,0.15);
  
  StyleHist(p_PiPi,kGreen+2);
  p_PiPi->SetTitle(";cos#theta_{#pi^{-}};p value");
  p_PiPi->GetYaxis()->SetRangeUser(0,1);
  p_PiPi->Draw("h");


  // Draw accepted and rejected
  TCanvas *c_acc_rej = new TCanvas("c_acc_rej","c_acc_rej",800,800);
  TGaxis::SetMaxDigits(3);
  TPad *pad_acc_rej = new TPad("pad_acc_rej", "pad_acc_rej",0,0,1,1);
  StylePad(pad_acc_rej,0,0.12,0,0.15);

  h_acc_PiPi_cos_eff_corr->SetTitle(";cos#theta_{#pi^{-}};Entries");

  h_acc_PiPi_cos_eff_corr->Draw("h");
  h_rej_PiPi_cos_eff_corr->Draw("hsame");

  TH1F * acc_full = (TH1F*) h_acc_PiPi_cos_eff_corr->Clone();
  TH1F * acc_add  = new TH1F("acc_add","acc_add",nbins,-1,1);

  for (int i = 1; i < nbins / 2 + 1; i++)
  {
    float accepted = acc_full->GetBinContent(nbins + 1 - i);
    accepted += acc_full->GetBinContent(i);
    acc_add->SetBinContent(i,accepted);
    acc_add->SetBinContent(nbins+1-i,accepted);
  }
  StyleHist(acc_add,kRed);
  acc_add->Draw("hsame");

  TLegend *leg2 = new TLegend(0.15,0.75,0.45,0.85);
  leg2->SetLineColor(0);
  leg2->AddEntry(h_acc_PiPi_cos_eff_corr,"N Accepted","l");
  leg2->AddEntry(h_rej_PiPi_cos_eff_corr,"N Rejected","l");
  leg2->AddEntry(acc_add,"N Accepted + opp. bin","l");
  leg2->Draw();

}

void pq_method_PiLPFO()
{
  try
  {
    if (!file->IsOpen()) return;

    main_pq();

  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }

}