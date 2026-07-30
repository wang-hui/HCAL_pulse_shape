int DLPHIN_illustrate () {

    auto hs = new THStack("hs","");

    auto h1 = new TH1F("h1","", 8, -0.5, 7.5);
    auto h2 = new TH1F("h2","", 8, -0.5, 7.5);
    auto h3 = new TH1F("h3","", 8, -0.5, 7.5);
    auto h4 = new TH1F("h4","", 8, -0.5, 7.5);
    auto h5 = new TH1F("h5","", 8, -0.5, 7.5);
    auto h6 = new TH1F("h6","", 8, -0.5, 7.5);

    auto tot = new TH1F("tot","", 8, -0.5, 7.5);

    std::vector<float> shape = {0.000371486, 0.620554, 0.202694, 0.0549495, 0.0238179};
    //std::vector<float> amplatude = {0.2, 0.1, 1.0, 0.2, 0.1, 0.15};
    std::vector<float> amplatude = {0.02, 0.01, 1.0, 0.02, 0.01, 0.015};

    for(int nbin = 0; nbin <= 4; nbin++) {
        h1->SetBinContent(nbin+1, shape.at(nbin)*amplatude.at(0));
        h2->SetBinContent(nbin+2, shape.at(nbin)*amplatude.at(1));
        h3->SetBinContent(nbin+3, shape.at(nbin)*amplatude.at(2));
        h4->SetBinContent(nbin+4, shape.at(nbin)*amplatude.at(3));
        h5->SetBinContent(nbin+5, shape.at(nbin)*amplatude.at(4));
        h6->SetBinContent(nbin+6, shape.at(nbin)*amplatude.at(5));
    }

    h1->SetLineStyle(4);
    h2->SetLineStyle(4);
    h3->SetLineStyle(4);
    h4->SetLineStyle(4);
    h5->SetLineStyle(4);
    h6->SetLineStyle(4);

    tot->SetMarkerStyle(20);
    tot->SetMarkerColor(kRed);

    hs->Add(h1);
    hs->Add(h2);
    hs->Add(h3);
    //hs->Add(h4);
    //hs->Add(h5);
    hs->Add(h6);

    auto cs = new TCanvas("cs","cs",600,300);
    cs->SetBottomMargin(0.2);
    cs->SetLeftMargin(0.15);

    auto leg = new TLegend(0.18,0.75,0.35,0.85);

    hs->Draw("hist");
    hs->GetXaxis()->SetTitle("Time sample");
    hs->GetXaxis()->SetTitleSize(0.1);
    hs->GetXaxis()->SetTitleOffset(1.0);
    hs->GetXaxis()->SetLabelSize(0.1);
    hs->GetXaxis()->SetNdivisions(8);

    hs->GetYaxis()->SetTitle("A.U.");
    hs->GetYaxis()->SetTitleSize(0.1);
    hs->GetYaxis()->SetTitleOffset(0.6);
    hs->GetYaxis()->SetLabelSize(0.1);
    hs->SetMaximum(0.8);

    for(int n = 1; n <=8; n++) {
        auto y = ((TH1*)hs->GetStack()->Last())->GetBinContent(n);
        cout << y << ", ";
        tot->SetBinContent(n, y);
    }
    tot->Draw("psame");

    leg->AddEntry(tot, "QIE data", "p");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.1);
    leg->Draw("same");

    cs->SaveAs("DLPHIN_cartoon.png");

    return 0;
}
