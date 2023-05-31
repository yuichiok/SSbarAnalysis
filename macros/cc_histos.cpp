void cc_histos(char dz, string mode) {
    TCanvas* canvas = new TCanvas("c","");
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(10);
    gStyle->SetPadColor(10);
    canvas->SetFrameLineWidth(2);
    canvas->SetFrameBorderMode(0);
    canvas->SetFrameBorderSize(0);
    canvas->SetLeftMargin(0.1);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.15);

    TVirtualPad *pad1 = canvas->cd(1);
    pad1->SetGrid(1,1);

    TFile* file = new TFile("../rootfiles/merged/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.cc.PFOp5.yevhenii.all.root");
    TH1F* h1 = (TH1F*)file->Get(Form("/Own/PS/%c0_P_%s",dz,mode.c_str()));
    TH1F* h2 = (TH1F*)file->Get(Form("/Own/PS/%c0_S_%s",dz,mode.c_str()));
        
    h1->Scale(1./h1->GetEntries());
    h1->SetTitle(Form("%c_{0} %s prong",dz,mode.c_str()));
    h1->GetXaxis()->SetTitle(Form("%c_{0}, [mm]",dz));
    h1->SetFillColor(kBlue);
    h1->SetFillStyle(3001);
    h1->SetLineColor(kBlue);
    h1->Draw("HIST");
    h2->Scale(1./h2->GetEntries());
    h2->SetLineColor(kRed);
    h2->SetFillStyle(3001);
    h2->SetFillColor(kRed);
    h2->Draw("HIST SAME");

    auto legend = new TLegend(0.7,0.6,0.90,0.80);
    legend->SetMargin(0.2);
    legend->AddEntry(h1, Form("P: mean=%.3f, #sigma=%.3f",h1->GetMean(),h1->GetStdDev()), "l");
    legend->AddEntry(h2, Form("S: mean=%.3f, #sigma=%.3f",h2->GetMean(),h2->GetStdDev()), "l");
    legend->Draw();
}