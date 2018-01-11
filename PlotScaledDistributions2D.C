// 'PlotScaledDistributions2D.C'
// Derek Anderson
// 10.04.2017
//
// Use this to scale some distributions, add em' up,
// and compare them to another distribution.
//
// NOTE: if 'embedNorm' is set to -1., then it will
//       scale each element of 'nEmbedTrg' by the
//       corresponding scale and add them up to give
//       the embedding normalization.
//
// NOTE: if 'fScale' is set to 5, it will calculate
//       the weights corresponding to each pTparton
//       bin and weight the distributions accordingly.


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;

static const UInt_t fScale(5);
static const UInt_t nDist(6);
static const UInt_t nHist(nDist + 1);
static const UInt_t nTotal(11);



void PlotScaledDistributions2D() {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning plot script..." << endl;


  // i/o parameters
  const TString sOutput("responseUnnormalized.pTcorr.r03a02rm1full.d11m12y2017.root");
  const TString sInput[nHist] = {"pp200r12pt9.match.r03a02rm1full.d11m12y2017.root", "pp200r12pt11.match.r03a02rm1full.d11m12y2017.root", "pp200r12pt15.match.r03a02rm1full.d11m12y2017.root", "pp200r12pt20.match.r03a02rm1full.d11m12y2017.root", "pp200r12pt25.match.r03a02rm1full.d11m12y2017.root", "pp200r12pt35.match.r03a02rm1full.d11m12y2017.root", "output/Tests/check.moreEvtsMatched.r03a02rm1chrg.d4m11y2017.root"};

  // histogram parameters
  const TString sHist[nHist] = {"hResponsePtc", "hResponsePtc", "hResponsePtc", "hResponsePtc", "hResponsePtc", "hResponsePtc", "hResponsePtc"};
  const TString sTitleD("Response matrix, p_{T}^{corr} (Pythia)");
  const TString sTitleS("Response matrix, p_{T}^{corr} (embedding)");
  const TString sTitleX("p_{T}^{corr}(detector)");
  const TString sTitleY("p_{T}^{corr}(particle)");

  // scales and misc. parameters
  const Double_t fudge[nTotal]    = {1.228, 1.051, 1.014, 1., 1., 1., 1., 1., 1., 1., 1.};
  const Double_t xSctn[nTotal]    = {9.005805, 1.461907, 0.3544354, 0.1513760, 2.488644e-02, 5.845846e-03, 2.304880e-03, 3.426618e-04, 4.562988e-05, 9.738044e-06, 5.019977e-07};
  const Double_t nEvts[nTotal]    = {2100295., 600300., 600300., 300289., 300289., 300289., 160295., 100302., 80293., 76303., 23307.};
  const Double_t nTrgEmbed[nDist] = {5., 21., 99., 287., 715., 652.};
  const Double_t dataNorm(100000.);
  const Double_t embedNorm(-1.);
  const Bool_t   doIntNorm(false);
  const Bool_t   doNormer(false);
  const Bool_t   doRebin(false);
  const UInt_t   nRebin(5);

  // constants
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fColLeg(0);
  const UInt_t  fStyLeg(0);
  const UInt_t  fLog(1);
  const UInt_t  fGrid(0);
  const UInt_t  nDec(2);
  const Float_t weight(1.);
  const Float_t fLblSize(0.02);
  const Float_t xPavN[2]  = {0.5, 0.7};
  const Float_t yPavN[2]  = {0.1, 0.3};
  const Float_t xPadD[2]  = {0., 0.5};
  const Float_t xPadS[2]  = {0.5, 1.};
  const Float_t yPad[2]   = {0., 1.};
  const Float_t width(1400);
  const Float_t height(700);
  const TString sSum("hSum");
  const TString sSumN("hSumNorm");
  const TString sCan("cHist");
  const TString sCanN("cNorm");
  const TString sPadD("pData");
  const TString sPadDN("pDataNorm");
  const TString sPadS("pSum");
  const TString sPadSN("pSumNorm");


  // open files and grab histograms
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput[nHist];
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    fInput[iHist]  = new TFile(sInput[iHist].Data(), "read");
    if (!fInput[iHist]) {
      cerr << "PANIC: couldn't open input file " << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Files opened." << endl;

  TH2D *hHist[nHist];
  TH2D *hNorm[nHist];
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    hHist[iHist] = (TH2D*) fInput[iHist] -> Get(sHist[iHist].Data()) -> Clone();
    hNorm[iHist] = (TH2D*) fInput[iHist] -> Get(sHist[iHist].Data()) -> Clone();
    if (!hHist[iHist]) {
      cerr << "PANIC: couldn't grab histogram " << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Histograms grabbed." << endl;

  // calculate scale set 5
  Double_t scales5[nDist];
  UInt_t   iStart = nTotal - nDist;
  for (UInt_t iDist = iStart; iDist < nTotal; iDist++) {
    const Double_t binLumi   = nEvts[iDist] / xSctn[iDist];
    const Double_t binWeight = 1. / (fudge[iDist] * binLumi);
    scales5[iDist - iStart]  = binWeight;
  }

  // select scales
  Double_t *scales;
  switch (fScale) {
    case 1:
      scales = scales1;
      break;
    case 2:
      scales = scales2;
      break;
    case 3:
      scales = scales3;
      break;
    case 4:
      scales = scales4;
      break;
    case 5:
      scales = scales5;
      break;
  }

  // norm check
  Double_t normer = 1.;
  if (doNormer) {
    TFile *fNormer = new TFile("output/scaleFactors.eTtrg.hotTowerCheck3.d16m10y2017.root", "read");
    TH1D  *hEmbedN = (TH1D*) fNormer -> Get("hSumUnormalized");
    TH1D  *hDataN  = (TH1D*) fNormer -> Get("hData");

    const Double_t embInt = hEmbedN -> Integral();
    const Double_t datInt = hDataN  -> Integral();
    normer = datInt / embInt;
    fNormer -> Close();
    fOutput -> cd();
  }

  // scale histograms
  Double_t nTrgTotal(0.);
  Double_t nTrgScale[nDist];
  for (UInt_t iHist = 0; iHist < nDist; iHist++) {
    if (iHist < nDist) {
      hHist[iHist] -> Scale(scales[iHist]);
      hNorm[iHist] -> Scale(scales[iHist]);
      hHist[iHist] -> Scale(normer);
      hNorm[iHist] -> Scale(normer);
      nTrgScale[iHist]  = nTrgEmbed[iHist] * scales[iHist] * normer;
      nTrgTotal        += nTrgScale[iHist];
    }
  }
  if (embedNorm != -1.) nTrgTotal = embedNorm;

  // normalize histograms
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    if (iHist == nDist)
      hNorm[iHist] -> Scale(1. / dataNorm);
    else
      hNorm[iHist] -> Scale(1. / nTrgTotal);
  }
  cout << "    Normalization:\n"
       << "      data  -- " << dataNorm <<"\n"
       << "      embed -- " << nTrgTotal
       << endl;
  
  fOutput -> cd();
  cout << "    Histograms scaled." << endl;


  // sum histograms
  const UInt_t  nBinsX = hHist[nDist] -> GetNbinsX();
  const UInt_t  nBinsY = hHist[nDist] -> GetNbinsY();
  const Float_t xBin1  = hHist[nDist] -> GetXaxis() -> GetBinLowEdge(1);
  const Float_t xBin2  = hHist[nDist] -> GetXaxis() -> GetBinLowEdge(nBinsX + 1);
  const Float_t yBin1  = hHist[nDist] -> GetYaxis() -> GetBinLowEdge(1);
  const Float_t yBin2  = hHist[nDist] -> GetYaxis() -> GetBinLowEdge(nBinsY + 1);

  TH2D *hSum  = new TH2D(sSum.Data(), "", nBinsX, xBin1, xBin2, nBinsY, yBin1, yBin2);
  TH2D *hSumN = new TH2D(sSumN.Data(), "", nBinsX, xBin1, xBin2, nBinsY, yBin1, yBin2);
  hSum  -> Sumw2();
  hSumN -> Sumw2();
  for (UInt_t iDist = 0; iDist < nDist; iDist++) {
    hSum  -> Add(hHist[iDist]);
    hSumN -> Add(hHist[iDist]);
  }
  hSumN -> Scale(1. / nTrgTotal);
  cout << "    Histograms summed." << endl;

  // quick fix
  const Double_t dIntN = hNorm[nDist] -> Integral();
  const Double_t eIntN = hSumN        -> Integral();
  if (doIntNorm) {
    hNorm[nDist] -> Scale(1. / dIntN);
    hSumN        -> Scale(1. / eIntN);
  }


  // set styles
  hHist[nDist] -> SetTitle(sTitleD.Data());
  hHist[nDist] -> SetTitleFont(fTxt);
  hHist[nDist] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hHist[nDist] -> GetXaxis() -> SetTitleFont(fTxt);
  hHist[nDist] -> GetXaxis() -> CenterTitle(fCnt);
  hHist[nDist] -> GetXaxis() -> SetLabelSize(fLblSize);
  hHist[nDist] -> GetYaxis() -> SetTitle(sTitleY.Data());
  hHist[nDist] -> GetYaxis() -> SetTitleFont(fTxt);
  hHist[nDist] -> GetYaxis() -> CenterTitle(fCnt);
  hHist[nDist] -> GetYaxis() -> SetLabelSize(fLblSize);
  hHist[nDist] -> GetZaxis() -> SetLabelSize(fLblSize);
  // normalized data
  hNorm[nDist] -> SetTitle(sTitleD.Data());
  hNorm[nDist] -> SetTitleFont(fTxt);
  hNorm[nDist] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hNorm[nDist] -> GetXaxis() -> SetTitleFont(fTxt);
  hNorm[nDist] -> GetXaxis() -> CenterTitle(fCnt);
  hNorm[nDist] -> GetXaxis() -> SetLabelSize(fLblSize);
  hNorm[nDist] -> GetYaxis() -> SetTitle(sTitleY.Data());
  hNorm[nDist] -> GetYaxis() -> SetTitleFont(fTxt);
  hNorm[nDist] -> GetYaxis() -> CenterTitle(fCnt);
  hNorm[nDist] -> GetYaxis() -> SetLabelSize(fLblSize);
  hNorm[nDist] -> GetZaxis() -> SetLabelSize(fLblSize);
  // summed distribution
  hSum -> SetTitle(sTitleS.Data());
  hSum -> SetTitleFont(fTxt);
  hSum -> GetXaxis() -> SetTitle(sTitleX.Data());
  hSum -> GetXaxis() -> SetTitleFont(fTxt);
  hSum -> GetXaxis() -> CenterTitle(fCnt);
  hSum -> GetXaxis() -> SetLabelSize(fLblSize);
  hSum -> GetYaxis() -> SetTitle(sTitleY.Data());
  hSum -> GetYaxis() -> SetTitleFont(fTxt);
  hSum -> GetYaxis() -> CenterTitle(fCnt);
  hSum -> GetYaxis() -> SetLabelSize(fLblSize);
  hSum -> GetZaxis() -> SetLabelSize(fLblSize);
  // summed (normalized) distribution
  hSumN -> SetTitle(sTitleS.Data());
  hSumN -> SetTitleFont(fTxt);
  hSumN -> GetXaxis() -> SetTitle(sTitleX.Data());
  hSumN -> GetXaxis() -> SetTitleFont(fTxt);
  hSumN -> GetXaxis() -> CenterTitle(fCnt);
  hSumN -> GetXaxis() -> SetLabelSize(fLblSize);
  hSumN -> GetYaxis() -> SetTitle(sTitleY.Data());
  hSumN -> GetYaxis() -> SetTitleFont(fTxt);
  hSumN -> GetYaxis() -> CenterTitle(fCnt);
  hSumN -> GetYaxis() -> SetLabelSize(fLblSize);
  hSumN -> GetZaxis() -> SetLabelSize(fLblSize);
  cout << "    Styles set." << endl;


  // make plots
  TCanvas *cHist = new TCanvas(sCan.Data(), "", width, height);
  TPad    *pData = new TPad(sPadD.Data(), "", xPadD[0], yPad[0], xPadD[1], yPad[1]);
  TPad    *pSum  = new TPad(sPadS.Data(), "", xPadS[0], yPad[0], xPadS[1], yPad[1]);
  pData        -> SetGrid(fGrid, fGrid);
  pData        -> SetLogz(fLog);
  pSum         -> SetGrid(fGrid, fGrid);
  pSum         -> SetLogz(fLog);
  pData        -> Draw();
  pSum         -> Draw();
  pData        -> cd();
  hHist[nDist] -> Draw("colz");
  pSum         -> cd();
  hSum         -> Draw("colz");
  cHist        -> Write();
  cHist        -> Close();

  TCanvas *cNorm  = new TCanvas(sCanN.Data(), "", width, height);
  TPad    *pDataN = new TPad(sPadDN.Data(), "", xPadD[0], yPad[0], xPadD[1], yPad[1]);
  TPad    *pSumN  = new TPad(sPadSN.Data(), "", xPadS[0], yPad[0], xPadS[1], yPad[1]);
  pDataN       -> SetGrid(fGrid, fGrid);
  pDataN       -> SetLogz(fLog);
  pSumN        -> SetGrid(fGrid, fGrid);
  pSumN        -> SetLogz(fLog);
  pDataN       -> Draw();
  pSumN        -> Draw();
  pDataN       -> cd();
  hNorm[nDist] -> Draw("colz");
  pSumN        -> cd();
  hSumN        -> Draw("colz");
  cNorm        -> Write();
  cNorm        -> Close();


  // set names of data and sum
  hSum         -> SetName("hSum");
  hSumN        -> SetName("hSumPerTrigger");
  hHist[nDist] -> SetName("hData");
  hNorm[nDist] -> SetName("hDataPerTrigger");


  // close files
  fOutput      -> cd();
  hSum         -> Write();
  hSumN        -> Write();
  hHist[nDist] -> Write();
  hNorm[nDist] -> Write();
  fOutput      -> Close();
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    fInput[iHist] -> cd();
    fInput[iHist] -> Close();
  }
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
