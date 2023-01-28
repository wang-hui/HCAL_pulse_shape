int M3_illustrate () {

    auto hs = new THStack("hs","");
    auto leg = new TLegend(0.7,0.6,0.9,0.9);

    auto h1 = new TH1F("h1","", 8, -0.5, 7.5);
    h1->SetFillColor(kAzure+1);
    auto h2 = new TH1F("h2","", 8, -0.5, 7.5);
    h2->SetFillColor(kYellow-10);
    auto h3 = new TH1F("h3","", 8, -0.5, 7.5);
    h3->SetFillColor(kPink+1);
    auto base = new TH1F("base","", 8, -0.5, 7.5);
    base->SetFillColor(kGray);

    std::vector<float> shape = {0.660921, 0.156969, 0.0469796};
    std::vector<float> amplatude = {2, 5, 1};

    for(int nbin = 1; nbin <= 8; nbin++) {
        if(nbin >= 3 && nbin <= 5){
            h1->SetBinContent(nbin, shape.at(nbin-3)*amplatude.at(0));
            h2->SetBinContent(nbin+1, shape.at(nbin-3)*amplatude.at(1));
            h3->SetBinContent(nbin+2, shape.at(nbin-3)*amplatude.at(2));
        }
        base->SetBinContent(nbin, 0.1);
    }

    hs->Add(base);
    hs->Add(h1);
    hs->Add(h2);
    hs->Add(h3);

    leg->AddEntry(base, "B", "f");
    leg->AddEntry(h1, "#mu_{SOI-1}", "f");
    leg->AddEntry(h2, "#mu_{SOI}", "f");
    leg->AddEntry(h3, "#mu_{SOI+1}", "f");

    auto cs = new TCanvas("cs","cs",600,400);
    cs->SetFrameLineWidth(0);   // remove top and right frame

    hs->Draw();
    hs->GetXaxis()->SetTitle("Time sample");
    hs->GetXaxis()->SetTitleSize(0.05);
    hs->GetXaxis()->SetLabelSize(0.05);
    hs->GetXaxis()->SetNdivisions(8);

    hs->GetYaxis()->SetLabelColor(0);   // hide lable
    hs->GetYaxis()->SetAxisColor(0);    // hide axis and tick
    hs->GetYaxis()->SetTickLength(0);   // such that the first white tick doesn't affect existing paint

    leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
    leg->Draw("same");

    for(int x = 0; x <= 3; x++) {
        auto l1 = new TLine(x+1.5, -0.2, x+1.5, 4);
        l1->SetLineStyle(2);
        l1->SetLineWidth(2);
        l1->Draw("same");
    }

    std::vector<TString> TSvec = {"#bf{A_{SOI-1}}", "#bf{A_{SOI}}", "#bf{A_{SOI+1}}"};
    for(int n = 0; n <=2; n++) {
        auto y = ((TH1*)hs->GetStack()->Last())->GetBinContent(n+3);
        auto t1 = new TLatex(n+1.65, y+0.2, TSvec.at(n));
        t1->Draw("same");
    }

    cs->SaveAs("method3_cartoon.png");
    cs->SaveAs("method3_cartoon.pdf");

    return 0;
}
