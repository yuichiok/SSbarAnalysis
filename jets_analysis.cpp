void jets_analysis() {
    TFile* file = new TFile("/home/padniuk/Desktop/SSbarAnalysis/rootfiles/tmp_root/output.root");
    TH1F* h1 = (TH1F*)file->Get("/Own/tagging/jets_info");
    cout<<"Events: = "<<h1->GetBinContent(1)<<endl;
    cout<<"Jet 1: = "<<h1->GetBinContent(2)<<endl;
    cout<<"Jet 2: = "<<h1->GetBinContent(3)<<endl;
    cout<<"not Jet 1/2: = "<<h1->GetBinContent(4)<<endl;
    cout<<"Jets with S: = "<<h1->GetBinContent(5)<<endl;
    cout<<"Jets w/out S: = "<<h1->GetBinContent(6)<<endl;
    cout<<"Events with both jets: = "<<h1->GetBinContent(7)<<endl;
    
    cout<<h1->GetBinContent(1)<<"\t"<<h1->GetBinContent(2)+h1->GetBinContent(3)+h1->GetBinContent(4)<<"\t"<<h1->GetBinContent(5)+h1->GetBinContent(6)<<"\t"<<endl;
} 
