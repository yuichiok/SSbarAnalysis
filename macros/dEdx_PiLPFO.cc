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

template <class T1>
void StyleHist1D(T1 *h, Color_t c)
{
  h->SetLineWidth(3);
  h->SetFillStyle(3002);
  h->SetMarkerColor(c);
  h->SetLineColor(c);
  h->SetFillColor(c);
}

template <class T2>
void StyleHist2D(T2 *h, Color_t c)
{
  h->SetMarkerColor(c);
  h->SetLineColor(c);
  h->SetFillColor(c);
}

template <class T3>
void StyleGraph(T3 *gr, Color_t c)
{
  gr->SetLineWidth(3);
  gr->SetMarkerSize(1);
  gr->SetMarkerStyle(2);
  gr->SetMarkerColor(c);
  gr->SetLineColor(c);
  gr->SetFillColor(c);
}

void Normalize(TH1F *h)
{
	double integral = h->Integral(1,h->GetNbinsX());
  h->Scale( 1 / integral );
}

double BetheBloch(const double *x, const double *pars){
  double mass=pars[5];
  double bg=x[0]/mass;
  double b=sqrt(bg*bg/(1.0+bg*bg));
  double tmax=pars[2]*TMath::Power(bg,2.0);   ///(1.0+pars[3]*g+pars[4]);

  return 1.350e-1*(0.5*pars[0]*TMath::Log(pars[1]*TMath::Power(bg,2.0)*tmax)-pars[3]*b*b-pars[4]*bg/2.0)/(b*b);
}

void dEdx_p(TFile *file)
{
  TCanvas *c0 = new TCanvas("c0","c0",800,800);
  TPad *pad0 = new TPad("pad0", "pad0",0,0,1,1);
  StylePad(pad0,0,0.15,0,0.17);

  TH2F *h2_gen_K_dEdx_p  = (TH2F*) file->Get("dEdx/h2_gen_K_reco_Pi_dEdx_p");
  TH2F *h2_gen_pi_dEdx_p = (TH2F*) file->Get("dEdx/h2_gen_pi_reco_Pi_dEdx_p");
  TH2F *h2_gen_p_dEdx_p  = (TH2F*) file->Get("dEdx/h2_gen_p_reco_Pi_dEdx_p");
  TH2F *h2_gen_e_dEdx_p  = (TH2F*) file->Get("dEdx/h2_gen_e_reco_Pi_dEdx_p");
  TH2F *h2_gen_mu_dEdx_p = (TH2F*) file->Get("dEdx/h2_gen_mu_reco_Pi_dEdx_p");

  StyleHist2D(h2_gen_K_dEdx_p,kRed);
  StyleHist2D(h2_gen_pi_dEdx_p,kBlue);
  StyleHist2D(h2_gen_p_dEdx_p,kGreen);
  StyleHist2D(h2_gen_e_dEdx_p,kViolet);
  StyleHist2D(h2_gen_mu_dEdx_p,kGray);

  pad0->SetLogx();
  h2_gen_pi_dEdx_p->SetTitle(";Track momentum [GeV];#frac{dE}{dx}[MeV]");
  h2_gen_pi_dEdx_p->GetXaxis()->SetTitleOffset(1.5);
  h2_gen_pi_dEdx_p->GetXaxis()->SetRangeUser(10,100);
  h2_gen_pi_dEdx_p->GetYaxis()->SetRangeUser(0.12,0.2);
  h2_gen_pi_dEdx_p->Draw("box");
  h2_gen_K_dEdx_p->Draw("box same");
  h2_gen_p_dEdx_p->Draw("box same");
  h2_gen_e_dEdx_p->Draw("box same");
  h2_gen_mu_dEdx_p->Draw("box same");

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
  leg->AddEntry(h2_gen_e_dEdx_p,"e^{#pm}","l");
  leg->AddEntry(h2_gen_mu_dEdx_p,"#mu^{#pm}","l");
  leg->AddEntry(func,"K Bethe-Bloch formula","l");
  leg->Draw();

}

void dEdx_dist_cos(TFile *file)
{
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  TPad *pad1 = new TPad("pad1", "pad1",0,0,1,1);
  StylePad(pad1,0,0.15,0,0.17);

  TH2F *h2_gen_K_KdEdx_dist_cos  = (TH2F*) file->Get("dEdx/h2_gen_K_reco_Pi_PidEdx_dist_cos");
  TH2F *h2_gen_pi_KdEdx_dist_cos = (TH2F*) file->Get("dEdx/h2_gen_pi_reco_Pi_PidEdx_dist_cos");
  TH2F *h2_gen_p_KdEdx_dist_cos  = (TH2F*) file->Get("dEdx/h2_gen_p_reco_Pi_PidEdx_dist_cos");
  TH2F *h2_gen_e_KdEdx_dist_cos  = (TH2F*) file->Get("dEdx/h2_gen_e_reco_Pi_PidEdx_dist_cos");
  TH2F *h2_gen_mu_KdEdx_dist_cos = (TH2F*) file->Get("dEdx/h2_gen_mu_reco_Pi_PidEdx_dist_cos");

  StyleHist2D(h2_gen_K_KdEdx_dist_cos,kRed);
  StyleHist2D(h2_gen_pi_KdEdx_dist_cos,kBlue);
  StyleHist2D(h2_gen_p_KdEdx_dist_cos,kGreen);
  StyleHist2D(h2_gen_e_KdEdx_dist_cos,kViolet);
  StyleHist2D(h2_gen_mu_KdEdx_dist_cos,kGray);

  h2_gen_pi_KdEdx_dist_cos->SetTitle(";cos#theta;dE/dx distance");
  h2_gen_pi_KdEdx_dist_cos->GetXaxis()->SetTitleOffset(1.5);
  h2_gen_pi_KdEdx_dist_cos->GetYaxis()->SetRangeUser(-5,5);
  h2_gen_pi_KdEdx_dist_cos->Draw("box");
  h2_gen_p_KdEdx_dist_cos->Draw("box same");
  h2_gen_K_KdEdx_dist_cos->Draw("box same");
  h2_gen_e_KdEdx_dist_cos->Draw("box same");
  h2_gen_mu_KdEdx_dist_cos->Draw("box same");

  TLegend *leg = new TLegend(0.5,0.76,0.75,0.85);
  leg->SetLineColor(0);
  leg->AddEntry(h2_gen_K_KdEdx_dist_cos,"K^{#pm}","l");
  leg->AddEntry(h2_gen_pi_KdEdx_dist_cos,"#pi^{#pm}","l");
  leg->AddEntry(h2_gen_p_KdEdx_dist_cos,"p","l");
  leg->AddEntry(h2_gen_e_KdEdx_dist_cos,"e^{#pm}","l");
  leg->AddEntry(h2_gen_mu_KdEdx_dist_cos,"#mu^{#pm}","l");
  leg->Draw();

}

void dEdx_dist_cos_proj(TFile *file)
{

  TH2F *hK  = (TH2F*) file->Get("dEdx/h2_gen_K_reco_Pi_PidEdx_dist_cos");
  TH2F *hpi = (TH2F*) file->Get("dEdx/h2_gen_pi_reco_Pi_PidEdx_dist_cos");
  TH2F *hp  = (TH2F*) file->Get("dEdx/h2_gen_p_reco_Pi_PidEdx_dist_cos");
  TH2F *he  = (TH2F*) file->Get("dEdx/h2_gen_e_reco_Pi_PidEdx_dist_cos");
  TH2F *hmu = (TH2F*) file->Get("dEdx/h2_gen_mu_reco_Pi_PidEdx_dist_cos");

  Int_t nslices = 3;

  TCanvas *c2 = new TCanvas("c2","c2",2400,600);
  c2->Divide(nslices,1);
  TPad *pad2s[nslices];

  TH1F * hK_proj[nslices]; 
  TH1F * hpi_proj[nslices];
  TH1F * hp_proj[nslices]; 
  TH1F * he_proj[nslices]; 
  TH1F * hmu_proj[nslices]; 
  
  for ( int islice=0; islice < nslices; islice++ ){

    c2->cd(islice+1);
    StylePad(gPad,0,0.15,0,0.17);

    Int_t binL = 93 - islice * 25 + 1;
    // Int_t binL = 3 + islice * 25 + 1;
    Int_t binH = binL + 1;

    hK_proj[islice]  = (TH1F*) hK->ProjectionY(TString::Format("hK_proj_%d",islice).Data(),binL,binH);
    hpi_proj[islice] = (TH1F*) hpi->ProjectionY(TString::Format("hpi_proj_%d",islice).Data(),binL,binH);
    hp_proj[islice]  = (TH1F*) hp->ProjectionY(TString::Format("hp_proj_%d",islice).Data(),binL,binH);
    he_proj[islice]  = (TH1F*) hp->ProjectionY(TString::Format("he_proj_%d",islice).Data(),binL,binH);
    hmu_proj[islice] = (TH1F*) hp->ProjectionY(TString::Format("hmu_proj_%d",islice).Data(),binL,binH);

    StyleHist1D(hK_proj[islice],kRed);
    StyleHist1D(hpi_proj[islice],kBlue);
    StyleHist1D(hp_proj[islice],kGreen+1);
    StyleHist1D(he_proj[islice],kViolet);
    StyleHist1D(hmu_proj[islice],kGray);

    Float_t binL_low  = hK->GetXaxis()->GetBinLowEdge(binL);
    Float_t binH_high = hK->GetXaxis()->GetBinLowEdge(binH+1);

    hpi_proj[islice]->SetTitle(TString::Format("Slice %.2f < cos#theta < %.2f;dE/dx distance [MeV];a.u.",binL_low,binH_high).Data());
    if(islice==0) {
      hpi_proj[islice]->GetYaxis()->SetRangeUser(0,0.4E4);
    }else{
      hpi_proj[islice]->GetYaxis()->SetRangeUser(0,0.3E4);
    }
    hpi_proj[islice]->GetXaxis()->SetRangeUser(-5,5);
    hpi_proj[islice]->Draw("h");
    hK_proj[islice]->Draw("hsame");
    hp_proj[islice]->Draw("hsame");
    he_proj[islice]->Draw("hsame");
    hmu_proj[islice]->Draw("hsame");

    if ( islice == 0 ){
      TLegend *leg = new TLegend(0.25,0.7,0.45,0.85);
      leg->SetLineColor(0);
      leg->AddEntry(hK_proj[islice],"K^{#pm}","l");
      leg->AddEntry(hpi_proj[islice],"#pi^{#pm}","l");
      leg->AddEntry(hp_proj[islice],"p","l");
      leg->AddEntry(he_proj[islice],"e^{#pm}","l");
      leg->AddEntry(hmu_proj[islice],"#mu^{#pm}","l");
      leg->Draw();
    }

  }

}

void dEdx_dist_cos_eff_pur(TFile *file)
{
  vector<TGraph*> graphs;

  TH2F *hK_gen  = (TH2F*) file->Get("dEdx/h2_gen_K_PidEdx_dist_cos");
  TH2F *hpi_gen = (TH2F*) file->Get("dEdx/h2_gen_pi_PidEdx_dist_cos");
  TH2F *hp_gen  = (TH2F*) file->Get("dEdx/h2_gen_p_PidEdx_dist_cos");
  TH2F *he_gen  = (TH2F*) file->Get("dEdx/h2_gen_e_PidEdx_dist_cos");
  TH2F *hmu_gen = (TH2F*) file->Get("dEdx/h2_gen_mu_PidEdx_dist_cos");

  TH2F *hK_sel  = (TH2F*) file->Get("dEdx/h2_gen_K_reco_Pi_PidEdx_dist_cos");
  TH2F *hpi_sel = (TH2F*) file->Get("dEdx/h2_gen_pi_reco_Pi_PidEdx_dist_cos");
  TH2F *hp_sel  = (TH2F*) file->Get("dEdx/h2_gen_p_reco_Pi_PidEdx_dist_cos");
  TH2F *he_sel  = (TH2F*) file->Get("dEdx/h2_gen_e_reco_Pi_PidEdx_dist_cos");
  TH2F *hmu_sel = (TH2F*) file->Get("dEdx/h2_gen_mu_reco_Pi_PidEdx_dist_cos");

  Int_t NbinsX = hpi_sel->GetNbinsX();
  Int_t NbinsY = hpi_sel->GetNbinsY();

  Float_t x[NbinsX], eff[NbinsX], pur[NbinsX];

  for ( int ibin=1; ibin <= NbinsX; ibin++ ){

    Int_t iarr = ibin - 1;

    x[iarr] = hK_sel->GetXaxis()->GetBinCenter(ibin);

    Float_t npions_all = hpi_gen->Integral(ibin,ibin,1,NbinsY);

    Float_t nkaons     = hK_sel->Integral(ibin,ibin,1,NbinsY);
    Float_t npions     = hpi_sel->Integral(ibin,ibin,1,NbinsY);
    Float_t nprotons   = hp_sel->Integral(ibin,ibin,1,NbinsY);
    Float_t nelectrons = he_sel->Integral(ibin,ibin,1,NbinsY);
    Float_t nmuons     = hmu_sel->Integral(ibin,ibin,1,NbinsY);

    eff[iarr] = npions / npions_all;
    pur[iarr] = npions / (nkaons + npions + nprotons + nelectrons + nmuons) ;

  }

  TGraph *gr_eff = new TGraph(NbinsX, x, eff);
  TGraph *gr_pur = new TGraph(NbinsX, x, pur);
  TGraph *gr_eff_pur = new TGraph(NbinsX, eff, pur);

  TCanvas * c_gr_eff = new TCanvas("c_gr_eff","c_gr_eff",800,800);
  gr_eff->SetTitle(";cos#theta;Efficiency");
  gr_eff->Draw("AC*");

  TCanvas * c_gr_pur = new TCanvas("c_gr_pur","c_gr_pur",800,800);
  gr_pur->SetTitle(";cos#theta;Purity");
  gr_pur->GetYaxis()->SetRangeUser(0,1);
  gr_pur->Draw("AC*");

}

void dEdx_PiLPFO()
{
  gStyle->SetOptStat(0);

  // TFile *uu_file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.PiLPFO.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");
  // TFile *uu_file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.ss.PiLPFO.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");

  TFile *file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.uu.KPiLPFO.distPi0.PFOp15.LPFOp15_pNaN.tpc0.hists.all.root","READ");

  dEdx_p(file);
  dEdx_dist_cos(file);
  dEdx_dist_cos_proj(file);
  dEdx_dist_cos_eff_pur(file);

}