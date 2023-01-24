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

Double_t func(Double_t *val, Double_t *par)
{
  Float_t accepted = val[1];
  Float_t rejected = val[0];
  Double_t a = 1;
  Double_t b = -1;
  Double_t c = rejected / (2 * (accepted + rejected));
  Double_t p = (0.5 / a) * (-b + sqrt(b * b - 4 * a * c));

  return p;
}

Double_t line(Double_t *val, Double_t *par)
{
  Float_t  x = val[0];
  Double_t y = 109.0*x / 91.0;

  return y;
}

void p_visual()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
  TCanvas *c0 = new TCanvas("c1","c1",800,800);
  gPad->SetGrid(1,1);
  StylePad(gPad,0,0,0.17,0.1);

  auto f   = new TF2("f",func,4000,21000,4000,21000);
  auto lin = new TF2("f",line,4000,21000,4000,21000);
  f->SetTitle(";Rejected;Accepted;p value");
  f->GetZaxis()->SetTitleOffset(1.5);
  f->Draw("colz");
  // f->Draw("surf2");
  // lin->Draw("surf2");


}