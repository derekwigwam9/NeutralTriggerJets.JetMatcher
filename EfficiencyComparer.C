// 'EfficiencyComparer.C'
// Derek Anderson
//
// Reads in a efficiency from a text file and plots it on top of a
// specified histogram.  Format is as so:
//
//   'x y dY'
//
// Where 'x' is the bin center, 'y' the bin content, and 'dY' is
// the error on 'y'.


#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

using namespace std;



void EfficiencyComparer() {

  // open files
  TFile *fOut = new TFile("effRecoR03.hJetXgJet.pTbinOneQt0151500.d10m7y2018.root", "recreate");
  TFile *fEff = new TFile("pp200r9embed.forCollabMeetingJul2018_qTtest15.pTbinOne.et9vz55.r03a02rm1chrg.dr03q15.root", "read");

  // grab comparison hist
  TH1D *hEffG = (TH1D*) fEff -> Get("hEfficiencyAll");
  if (!hEffG) {
    cerr << "PANIC: couldn't grab comparison histogram!" << endl;
    return;
  }
  hEffG -> SetName("hEffG");
  hEffG -> SetTitle("Matching efficiency, R = 0.3");


  vector<Double_t> b;
  vector<Double_t> c;
  vector<Double_t> w;
  vector<Double_t> e;

  // stream in data
  ifstream eff("input/hJetGroomed.effReco6080r03.d9m7y2018.txt");
  if (!eff) {
    cerr << "PANIC: input stream couldn't be opened!" << std::endl;
    return;
  }
  else {

    b.clear();
    c.clear();
    w.clear();
    e.clear();

    Double_t x  = 0.;
    Double_t y  = 0.;
    Double_t dY = 0.;
    while (eff) {
      eff >> x;
      eff >> y;
      eff >> dY;
      b.push_back(x);
      c.push_back(y);
      w.push_back(0.);
      e.push_back(dY);
    }

  }  // end stream

  // store data in arrays
  const Int_t nBin = (Int_t) b.size();
  Double_t bin[nBin];
  Double_t cnt[nBin];
  Double_t siz[nBin];
  Double_t err[nBin];
  for (Int_t i = 0; i < nBin; ++i) {
    bin[i] = b[i];
    cnt[i] = c[i];
    siz[i] = w[i];
    err[i] = e[i];
  }


  // create TGraph
  TGraphErrors *gEffH = new TGraphErrors(nBin, bin, cnt, siz, err);
  gEffH -> SetNameTitle("gEffH", "Matching efficiency, R = 0.3");

  // create histogram
  const Int_t    bins = 62;
  const Double_t bin1 = -0.75;
  const Double_t bin2 = 30.25;
  TH1D *hEffH = new TH1D("hEffH", "Matching efficiency, R = 0.3", bins, bin1, bin2);
  for (Int_t i = 0; i < nBin; ++i) {
    const Int_t iBin = hEffH -> FindBin(bin[i]);
    hEffH -> SetBinContent(iBin, cnt[i]);
    hEffH -> SetBinError(iBin, err[i]);
  }

  // create plot
  const Int_t   fTxt    = 42;
  const Int_t   fAln    = 12;
  const Int_t   fCnt    = 1;
  const Int_t   fColL   = 0;
  const Int_t   fCol[2] = {810, 1};
  const Int_t   fLin[2] = {1, 1};
  const Int_t   fFil[2] = {3001, 0};
  const Int_t   fMar[2] = {1, 8};
  const Float_t fLab    = 0.03;
  const Float_t fOffX   = 1.;

  hEffH -> SetLineStyle(fLin[0]);
  hEffH -> SetLineColor(fCol[0]);
  hEffH -> SetFillStyle(fFil[0]);
  hEffH -> SetFillColor(fCol[0]);
  hEffH -> SetMarkerStyle(fMar[0]);
  hEffH -> SetMarkerColor(fCol[0]);
  hEffH -> SetTitleFont(fTxt);
  hEffH -> GetXaxis() -> SetTitle("p_{T}^{jet} [GeV/c]");
  hEffH -> GetXaxis() -> SetTitleFont(fTxt);
  hEffH -> GetXaxis() -> SetTitleOffset(fOffX);
  hEffH -> GetXaxis() -> SetLabelFont(fTxt);
  hEffH -> GetXaxis() -> SetLabelSize(fLab);
  hEffH -> GetXaxis() -> CenterTitle(fCnt);
  hEffH -> GetYaxis() -> SetTitle("#epsilon_{match}(p_{T}^{jet})");
  hEffH -> GetYaxis() -> SetTitleFont(fTxt);
  hEffH -> GetYaxis() -> SetLabelFont(fTxt);
  hEffH -> GetYaxis() -> SetLabelSize(fLab);
  hEffH -> GetYaxis() -> CenterTitle(fCnt);

  hEffG -> SetLineStyle(fLin[1]);
  hEffG -> SetLineColor(fCol[1]);
  hEffG -> SetFillStyle(fFil[1]);
  hEffG -> SetFillColor(fCol[1]);
  hEffG -> SetMarkerStyle(fMar[1]);
  hEffG -> SetMarkerColor(fCol[1]);
  hEffG -> SetTitleFont(fTxt);
  hEffG -> GetXaxis() -> SetTitle("p_{T}^{jet} [GeV/c]");
  hEffG -> GetXaxis() -> SetTitleFont(fTxt);
  hEffG -> GetXaxis() -> SetTitleOffset(fOffX);
  hEffG -> GetXaxis() -> SetLabelFont(fTxt);
  hEffG -> GetXaxis() -> SetLabelSize(fLab);
  hEffG -> GetXaxis() -> CenterTitle(fCnt);
  hEffG -> GetYaxis() -> SetTitle("#epsilon_{match}(p_{T}^{jet})");
  hEffG -> GetYaxis() -> SetTitleFont(fTxt);
  hEffG -> GetYaxis() -> SetLabelFont(fTxt);
  hEffG -> GetYaxis() -> SetLabelSize(fLab);
  hEffG -> GetYaxis() -> CenterTitle(fCnt);

  TLegend *leg = new TLegend(0.1, 0.1, 0.3, 0.3);
  leg -> SetFillColor(fColL);
  leg -> SetLineColor(fColL);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(hEffH, "h^{#pm} jet, AuAu (60-80%)");
  leg -> AddEntry(hEffG, "#pi^{0}-jet, pp");

  TCanvas *cEff = new TCanvas("cEff", "", 750, 750);
  hEffH -> Draw("C E6");
  hEffG -> Draw("same");
  leg   -> Draw();
  fOut  -> cd();
  cEff  -> Write();
  cEff  -> Close();


  // save and close files
  fOut  -> cd();
  gEffH -> Write();
  hEffH -> Write();
  hEffG -> Write();
  fOut  -> Close();
  fEff  -> cd();
  fEff  -> Close();

}

// End ------------------------------------------------------------------------
