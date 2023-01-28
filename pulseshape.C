//void computeHPDShape(float ts1, float ts2, float ts3, float thpd, float tpre, float wd1, float wd2, float wd3, Shape& tmphpdShape_) 
void computeHPDShape(float ts1, float ts2, float ts3, float thpd, float tpre, float wd1, float wd2, float wd3, TH1D* &h1)
{
  // pulse shape time constants in ns
  /*
  const float ts1  = 8.;          // scintillation time constants : 1,2,3
  const float ts2  = 10.;           
  const float ts3  = 29.3;         
  const float thpd = 4.;          // HPD current collection drift time
  const float tpre = 9.;          // preamp time constant (refit on TB04 data)
  
  const float wd1 = 2.;           // relative weights of decay exponents 
  const float wd2 = 0.7;
  const float wd3 = 1.;
*/
  // pulse shape components over a range of time 0 ns to 255 ns in 1 ns steps
  unsigned int nbin = 256;
  //tmphpdShape_.setNBin(nbin);
  std::vector<float> ntmp(nbin, 0.0);  // zeroing output pulse shape
  std::vector<float> nth(nbin, 0.0);   // zeroing HPD drift shape
  std::vector<float> ntp(nbin, 0.0);   // zeroing Binkley preamp shape
  std::vector<float> ntd(nbin, 0.0);   // zeroing Scintillator decay shape

  unsigned int i, j, k;
  float norm;

  // HPD starts at I and rises to 2I in thpd of time
  norm = 0.0;
  for (j = 0; j < thpd && j < nbin; j++) {
    nth[j] = 1.0 + ((float)j) / thpd;
    norm += nth[j];
  }
  // normalize integrated current to 1.0
  for (j = 0; j < thpd && j < nbin; j++) {
    nth[j] /= norm;
  }

  // Binkley shape over 6 time constants
  norm = 0.0;
  for (j = 0; j < 6 * tpre && j < nbin; j++) {
    ntp[j] = ((float)j) * exp(-((float)(j * j)) / (tpre * tpre));
    norm += ntp[j];
  }
  // normalize pulse area to 1.0
  for (j = 0; j < 6 * tpre && j < nbin; j++) {
    ntp[j] /= norm;
  }

  // ignore stochastic variation of photoelectron emission
  // <...>

  // effective tile plus wave-length shifter decay time over 4 time constants
  unsigned int tmax = 6 * (int)ts3;

  norm = 0.0;
  for (j = 0; j < tmax && j < nbin; j++) {
    ntd[j] = wd1 * exp(-((float)j) / ts1) + wd2 * exp(-((float)j) / ts2) + wd3 * exp(-((float)j) / ts3);
    norm += ntd[j];
  }
  // normalize pulse area to 1.0
  for (j = 0; j < tmax && j < nbin; j++) {
    ntd[j] /= norm;
  }

  unsigned int t1, t2, t3, t4;
  for (i = 0; i < tmax && i < nbin; i++) {
    t1 = i;
    //    t2 = t1 + top*rand;
    // ignoring jitter from optical path length
    t2 = t1;
    for (j = 0; j < thpd && j < nbin; j++) {
      t3 = t2 + j;
      for (k = 0; k < 4 * tpre && k < nbin; k++) {  // here "4" is set deliberately,
        t4 = t3 + k;                                // as in test fortran toy MC ...
        if (t4 < nbin) {
          unsigned int ntb = t4;
          ntmp[ntb] += ntd[i] * nth[j] * ntp[k];
        }
      }
    }
  }

  // normalize for 1 GeV pulse height
  norm = 0.;
  for (i = 0; i < nbin; i++) {
    norm += ntmp[i];
  }

  for (i = 0; i < nbin; i++) {
    ntmp[i] /= norm;
  }

  //for (i = 0; i < nbin; i++) {
  //for (i = 0; i < 250; i++) {
  for (i = 0; i < 125; i++) {
    //h1->SetBinContent(i+1, ntmp[i]);
    h1->SetBinContent(i+18, ntmp[i]);
  }

}



void pulseshape()
{ 
 float SiPMShapeData2018[250] = 
      {5.22174e-12, 7.04852e-10, 3.49584e-08, 7.78029e-07, 9.11847e-06, 6.39666e-05, 0.000297587, 0.000996661,
       0.00256618,  0.00535396,  0.00944073,  0.0145521,   0.020145,    0.0255936,   0.0303632,   0.0341078,
       0.0366849,   0.0381183,   0.0385392,   0.0381327,   0.0370956,   0.0356113,   0.0338366,   0.0318978,
       0.029891,    0.0278866,   0.0259336,   0.0240643,   0.0222981,   0.0206453,   0.0191097,   0.0176902,
       0.0163832,   0.0151829,   0.0140826,   0.0130752,   0.0121533,   0.01131,     0.0105382,   0.00983178,
       0.00918467,  0.00859143,  0.00804709,  0.0075471,   0.00708733,  0.00666406,  0.00627393,  0.00591389,
       0.00558122,  0.00527344,  0.00498834,  0.00472392,  0.00447837,  0.00425007,  0.00403754,  0.00383947,
       0.00365465,  0.00348199,  0.00332052,  0.00316934,  0.00302764,  0.0028947,   0.00276983,  0.00265242,
       0.00254193,  0.00243785,  0.00233971,  0.00224709,  0.0021596,   0.00207687,  0.0019986,   0.00192447,
       0.00185421,  0.00178756,  0.0017243,   0.00166419,  0.00160705,  0.00155268,  0.00150093,  0.00145162,
       0.00140461,  0.00135976,  0.00131696,  0.00127607,  0.00123699,  0.00119962,  0.00116386,  0.00112963,
       0.00109683,  0.0010654,   0.00103526,  0.00100634,  0.000978578, 0.000951917, 0.000926299, 0.000901672,
       0.000877987, 0.000855198, 0.00083326,  0.000812133, 0.000791778, 0.000772159, 0.000753242, 0.000734994,
       0.000717384, 0.000700385, 0.000683967, 0.000668107, 0.000652779, 0.00063796,  0.000623629, 0.000609764,
       0.000596346, 0.000583356, 0.000570777, 0.000558592, 0.000546785, 0.00053534,  0.000524243, 0.000513481,
       0.00050304,  0.000492907, 0.000483072, 0.000473523, 0.000464248, 0.000455238, 0.000446483, 0.000437974,
       0.0004297,   0.000421655, 0.00041383,  0.000406216, 0.000398807, 0.000391595, 0.000384574, 0.000377736,
       0.000371076, 0.000364588, 0.000358266, 0.000352104, 0.000346097, 0.00034024,  0.000334528, 0.000328956,
       0.00032352,  0.000318216, 0.000313039, 0.000307986, 0.000303052, 0.000298234, 0.000293528, 0.000288931,
       0.000284439, 0.00028005,  0.000275761, 0.000271567, 0.000267468, 0.000263459, 0.000259538, 0.000255703,
       0.000251951, 0.00024828,  0.000244688, 0.000241172, 0.00023773,  0.000234361, 0.000231061, 0.00022783,
       0.000224666, 0.000221566, 0.000218528, 0.000215553, 0.000212636, 0.000209778, 0.000206977, 0.00020423,
       0.000201537, 0.000198896, 0.000196307, 0.000193767, 0.000191275, 0.000188831, 0.000186432, 0.000184079,
       0.000181769, 0.000179502, 0.000177277, 0.000175092, 0.000172947, 0.000170841, 0.000168772, 0.000166741,
       0.000164745, 0.000162785, 0.000160859, 0.000158967, 0.000157108, 0.00015528,  0.000153484, 0.000151719,
       0.000149984, 0.000148278, 0.000146601, 0.000144951, 0.000143329, 0.000141734, 0.000140165, 0.000138622,
       0.000137104, 0.00013561,  0.000134141, 0.000132695, 0.000131272, 0.000129871, 0.000128493, 0.000127136,
       0.000125801, 0.000124486, 0.000123191, 0.000121917, 0.000120662, 0.000119426, 0.000118209, 0.00011701,
       0.000115829, 0.000114665, 0.000113519, 0.00011239,  0.000111278, 0.000110182, 0.000109102, 0.000108037,
       0.000106988, 0.000105954, 0.000104935, 0.00010393,  0.000102939, 0.000101963, 0.000101,    0.000100051,
       9.91146e-05, 9.81915e-05, 9.7281e-05,  9.63831e-05, 9.54975e-05, 9.46239e-05, 9.37621e-05, 9.2912e-05,
       9.20733e-05, 9.12458e-05};

  //TH1D *h1_2018he = new TH1D("h1_2018he", "h1_2018he", 250, 50, 300);
  TH1D *h1_2018he = new TH1D("h1_2018he", "h1_2018he", 125, 50, 175);
  for(unsigned int i=1; i<=250; i++)
  {
    //h1_2018he->SetBinContent(i, SiPMShapeData2018[i-1]);
    h1_2018he->SetBinContent(i+14, SiPMShapeData2018[i-1]);
  }

  // HB
  //TH1D *h1_2018hb = new TH1D("h1_2018hb", "h1_2018hb", 250, 50, 300);
  TH1D *h1_2018hb = new TH1D("h1_2018hb", "h1_2018hb", 125, 50, 175);
  //  HPD Shape  Version 3 for CMSSW 5. Nov 2011  (RECO and MC separately)
  float ts1 = 8.;
  float ts2 = 19.;
  float ts3 = 29.3;
  float thpd = 4.0;
  float tpre = 9.0;
  float wd1 = 2.0;
  float wd2 = 0.7;
  float wd3 = 0.32; 
  computeHPDShape(ts1, ts2, ts3, thpd, tpre, wd1, wd2, wd3, h1_2018hb);
 
  TH1D* h1_2018he_rebin = (TH1D*)h1_2018he->Clone("h1_2018he_rebin");
  TH1D* h1_2018hb_rebin = (TH1D*)h1_2018hb->Clone("h1_2018hb_rebin");
  h1_2018he_rebin->Rebin(25);
  h1_2018hb_rebin->Rebin(25);
  
  cout << "he: " << h1_2018he_rebin->GetBinContent(2) << endl;
  cout << "hb: " << h1_2018hb_rebin->GetBinContent(2) << endl;

  int n_bins = h1_2018he_rebin->GetNbinsX();
  for (int i = 1; i <= n_bins; i++) {
    cout << h1_2018he_rebin->GetBinContent(i) << ", ";
  }

  
  TCanvas *c = new TCanvas("c","c",600,400);
  c->cd(1);

  h1_2018he_rebin->SetTitle("");
  h1_2018he_rebin->SetStats(0);
  h1_2018he_rebin->SetMaximum(1);
  h1_2018he_rebin->SetLineColor(kBlack);
  h1_2018he_rebin->SetFillColor(kYellow-10);
  h1_2018he_rebin->SetLineWidth(3);
  h1_2018he_rebin->Draw("hist");
  h1_2018he_rebin->GetXaxis()->SetLabelSize(0.04);
  h1_2018he_rebin->GetXaxis()->SetLabelOffset(0.01);
  h1_2018he_rebin->GetXaxis()->SetTitle("Time [ns]");
  h1_2018he_rebin->GetXaxis()->SetTitleSize(0.05); 
  //h1_2018he_rebin->GetXaxis()->SetTitleOffset(1.1); 
  h1_2018he_rebin->GetXaxis()->SetNdivisions(-505); 
  h1_2018he_rebin->GetYaxis()->SetLabelSize(0.04);
  h1_2018he_rebin->GetYaxis()->SetLabelOffset(0.01);
  h1_2018he_rebin->GetYaxis()->SetTitle("A.U.");
  h1_2018he_rebin->GetYaxis()->SetTitleSize(0.05); 
  //h1_2018he_rebin->GetYaxis()->SetTitleOffset(1.1); 
  h1_2018he_rebin->GetYaxis()->SetRangeUser(0, 1.4); 

  h1_2018he->Scale(25);
  h1_2018he->SetLineColor(kRed);
  h1_2018he->SetLineWidth(3);
  h1_2018he->Draw("hist same");

  // legend
  TLegend *leg = new TLegend(0.4, 0.7, 0.9, 0.9);
  leg->SetNColumns(1);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextAlign(12);
  leg->SetTextSize(0.04);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->AddEntry(h1_2018he, " Pulse shape with 1 ns granularity", "l");
  leg->AddEntry(h1_2018he_rebin, " Pulse integrated over 25 ns", "f");
  leg->Draw();

  // labels
  float textSize = 0.05;
  TLatex *TexEnergyLumi = new TLatex(0.9,0.92,"#font[42]{(13 TeV)}");
  TexEnergyLumi->SetNDC();
  TexEnergyLumi->SetTextSize(textSize);
  TexEnergyLumi->SetTextAlign (31);
  TexEnergyLumi->SetLineWidth(2);

  //TLatex *TexCMS = new TLatex(0.1,0.92,"CMS #font[52]{Preliminary} 2018");
  TLatex *TexCMS = new TLatex(0.13,0.81,"CMS");
  TexCMS->SetNDC();
  TexCMS->SetTextSize(textSize+0.03);
  TexCMS->SetLineWidth(2);
  
  TLatex *TexPrel = new TLatex(0.13,0.80,"#font[52]{Preliminary}");
  TexPrel->SetNDC();
  TexPrel->SetTextSize(textSize);
  TexPrel->SetLineWidth(2);

  TexEnergyLumi->Draw("same");
  TexCMS->Draw("same");
  //TexPrel->Draw("same");

  c->SaveAs("pulseshape.png");
  c->SaveAs("pulseshape.pdf");
}
