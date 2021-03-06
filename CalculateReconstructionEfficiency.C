// 'CalculateReconstructionEfficiency.C'
// Derek Anderson
// 05.15.2018
//
// Calculate the jet reconstruction efficiency
// from the output of 'MatchJets.C'
//
// NOTE: Configuration = 0, RFF configuration
//       Configuration = 1, FF configuration
//       Configuration = 2, both (1st RFF, then FF)


#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;


// constants
static const UInt_t   NHistFF(7);
static const UInt_t   NHistRFF(8);
static const UInt_t   NHistAll(NHistRFF + NHistFF);
static const UInt_t   NTotal(10);
static const Double_t NormsFF[NHistFF]   = {2., 25., 111., 357., 2039., 7688., 4063.};
static const Double_t NormsRFF[NHistRFF] = {1., 3., 26., 152., 515., 2350., 6681., 3549.};
static const Double_t NormsAll[NHistAll] = {1., 3., 26., 152., 515., 2350., 6681., 3549., 2., 25., 111., 357., 2039., 7688., 4063.};
static const Double_t WeightsFF[NTotal]  = {1.0, 3.361596e-01, 1.401161e-01, 1.337302e-01, 2.895246e-02, 1.042577e-02, 8.294575e-03, 2.064352e-03, 8.088693e-05, 1.417116e-05};
static const Double_t WeightsRFF[NTotal] = {1.0, 3.501425e-01, 1.395103e-01, 1.326444e-01, 2.801546e-02, 1.031377e-02, 8.210314e-03, 1.985107e-03, 8.054588e-05, 1.449037e-05};

// options
static const UInt_t NHist(NHistAll);
static const UInt_t Configuration(2);
static const Bool_t UseRootTrgs(true);
static const Bool_t DoIntNorm(false);
static const Bool_t VariableBins(true);



void CalculateReconstructionEfficiency() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning reconstruction efficiency calcualation..." << endl;

  // io parameters
  const TString sOut("pp200r9eff.weightingCheckNoWeight.pTbinVar.et920vz55.r03a02rm1chrg.root");
  const TString sIn[NHist] = {"output/CollabMeetingJul2018/Et920/pp200r9pt4rff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt5rff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt7rff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt9rff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt11rff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt15rff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt25rff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt35rff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt5ff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt7ff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt9ff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt11ff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt15ff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt25ff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root", "output/CollabMeetingJul2018/Et920/pp200r9pt35ff.forCollabMeetingJul2018.pTbinVar.et920vz55.r03a02rm1chrg.dr03q15.root"};
  const TString sPar[NHist] = {"EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr", "EventInfo/hParPtCorr"};
  const TString sDet[NHist] = {"EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr", "EventInfo/hDetPtCorr"};
  const TString sTrg[NHist] = {"EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD", "EventInfo/hRefmultD"};


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fIn[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fIn[iHist] = new TFile(sIn[iHist].Data(), "read");
    if (!fIn[iHist]) {
      cerr << "PANIC: couldn't open input file no. " << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // get histograms
  TH1D *hPar[NHist];
  TH1D *hDet[NHist];
  TH1D *hTrg[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hPar[iHist] = (TH1D*) fIn[iHist] -> Get(sPar[iHist].Data());
    hDet[iHist] = (TH1D*) fIn[iHist] -> Get(sDet[iHist].Data());
    hTrg[iHist] = (TH1D*) fIn[iHist] -> Get(sTrg[iHist].Data());
    if (!hPar[iHist] || !hDet[iHist] || !hTrg[iHist]) {
      cerr << "PANIC: couldn't grab a histogram!\n"
           << "       hPar[" << iHist << "] = " << hPar[iHist] << ", hDet[" << iHist << "] = " << hDet[iHist] << ", hTrg[" << iHist << "] = " << hTrg[iHist]
           << endl;
      return;
    }
  }
  fOut -> cd();
  cout << "    Grabbed histograms." << endl;


  // scale histograms
  Double_t normerFF(0.);
  Double_t normerRFF(0.);
  Double_t normerAll(0.);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    switch (Configuration) {
      case 0:
        hPar[iHist] -> Scale(WeightsRFF[(NTotal - NHist) + iHist]);
        hDet[iHist] -> Scale(WeightsRFF[(NTotal - NHist) + iHist]);
        if (UseRootTrgs)
          normerRFF += ((Double_t) hTrg[iHist] -> GetEntries()) * WeightsRFF[(NTotal - NHist) + iHist];
        else
          normerRFF += NormsRFF[iHist] * WeightsRFF[(NTotal - NHist) + iHist];
        break;
      case 1:
        hPar[iHist] -> Scale(WeightsFF[(NTotal - NHist) + iHist]);
        hDet[iHist] -> Scale(WeightsFF[(NTotal - NHist) + iHist]);
        if (UseRootTrgs)
          normerFF += ((Double_t) hTrg[iHist] -> GetEntries()) * WeightsFF[(NTotal - NHist) + iHist];
        else
          normerFF += NormsFF[iHist] * WeightsFF[(NTotal - NHist) + iHist];
        break;
      case 2:
        const Bool_t isRFF = (iHist < NHistRFF);
        if (isRFF) {
          hPar[iHist] -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHist]);
          hDet[iHist] -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHist]);
          if (UseRootTrgs) {
            normerRFF += ((Double_t) hTrg[iHist] -> GetEntries()) * WeightsRFF[(NTotal - NHistRFF) + iHist];
            normerAll += ((Double_t) hTrg[iHist] -> GetEntries()) * WeightsRFF[(NTotal - NHistRFF) + iHist];
          }
          else {
            normerRFF += NormsAll[iHist] * WeightsRFF[(NTotal - NHistRFF) + iHist];
            normerAll += NormsAll[iHist] * WeightsRFF[(NTotal - NHistRFF) + iHist];
          }
        }
        else {
          hPar[iHist] -> Scale(WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)]);
          hDet[iHist] -> Scale(WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)]);
          if (UseRootTrgs) {
            normerFF  += ((Double_t) hTrg[iHist] -> GetEntries()) * WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)];
            normerAll += ((Double_t) hTrg[iHist] -> GetEntries()) * WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)];
          }
          else {
            normerFF  += NormsAll[iHist] * WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)];
            normerAll += NormsAll[iHist] * WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)];
          }
        }
        break;
    }
  }
  cout << "    Scaled histograms.\n"
       << "      norm(RFF, FF, All) = (" << normerRFF << ", " << normerFF << ", " << normerAll << ")"
       << endl;


  // sum histograms
  TH1D *hSumP = (TH1D*) hPar[0] -> Clone();
  TH1D *hSumD = (TH1D*) hDet[0] -> Clone();
  hSumP -> SetName("hSumPar");
  hSumD -> SetName("hSumDet");
  hSumP -> Reset("ICE");
  hSumD -> Reset("ICE");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    switch (Configuration) {
      const Double_t weightFF  = 1.;
      const Double_t weightRFF = 1.;
      case 0:
        hSumP -> Add(hPar[iHist], weightRFF);
        hSumD -> Add(hDet[iHist], weightRFF);
        break;
      case 1:
        hSumP -> Add(hPar[iHist], weightFF);
        hSumD -> Add(hDet[iHist], weightFF);
        break;
      case 2:
        const Bool_t isRFF = (iHist < NHistRFF);
        if (isRFF) {
          hSumP -> Add(hPar[iHist], weightRFF);
          hSumD -> Add(hDet[iHist], weightRFF);
        }
        else {
          hSumP -> Add(hPar[iHist], weightFF);
          hSumD -> Add(hDet[iHist], weightFF);
        }
        break;
    }
  }
  cout << "    Summed histograms." << endl;


  // normalize histograms
  if (VariableBins) {
    const UInt_t nBinP = hSumP -> GetNbinsX();
    const UInt_t nBinD = hSumD -> GetNbinsX();
    for (UInt_t iBinP = 1; iBinP < (nBinP + 1); iBinP++) {
      const Double_t pWidth = hSumP -> GetBinWidth(iBinP);
      const Double_t pVal   = hSumP -> GetBinContent(iBinP);
      const Double_t pErr   = hSumP -> GetBinError(iBinP);
      hSumP -> SetBinContent(iBinP, pVal / pWidth);
      hSumP -> SetBinError(iBinP, pErr / pWidth);
    }
    for (UInt_t iBinD = 1; iBinD < (nBinD + 1); iBinD++) {
      const Double_t dWidth = hSumD -> GetBinWidth(iBinD);
      const Double_t dVal   = hSumD -> GetBinContent(iBinD);
      const Double_t dErr   = hSumD -> GetBinError(iBinD);
      hSumD -> SetBinContent(iBinD, dVal / dWidth);
      hSumD -> SetBinError(iBinD, dErr / dWidth);
    }
  }

  if (DoIntNorm) {
    const Double_t pIntegral = hSumP -> Integral();
    const Double_t dIntegral = hSumD -> Integral();
    hSumP -> Scale(1. / pIntegral);
    hSumD -> Scale(1. / dIntegral);
  }
  else {
    switch (Configuration) {
      case 0:
        hSumP -> Scale(1. / normerRFF);
        hSumD -> Scale(1. / normerRFF);
        break;
      case 1:
        hSumP -> Scale(1. / normerFF);
        hSumD -> Scale(1. / normerFF);
        break;
      case 2:
        hSumP -> Scale(1. / normerAll);
        hSumD -> Scale(1. / normerAll);
        break;
    }
  }
  cout << "    Normalized histograms." << endl;


  // calculate efficiency
  TH1D *hEff = (TH1D*) hSumP -> Clone();
  hEff -> SetName("hEff");
  hEff -> Reset();
  hEff -> Divide(hSumD, hSumP, 1., 1.);
  cout << "    Calculated efficiency." << endl;


  // set styles
  const UInt_t  cSum[2] = {858, 898};
  const UInt_t  mSum[2] = {7, 4};
  const UInt_t  cEff(1);
  const UInt_t  mEff(8);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.1);
  const Float_t fLbl(0.02);
  const TString sSum("Recoil jet p_{T}^{corr}");
  const TString sSumX("p_{T}^{corr}(MC)");
  const TString sSumY("(1/N^{trg}_{eff}) dN^{jet}/dp_{T}^{corr}");
  const TString sEff("Reconstruction efficiency");
  const TString sEffX("p_{T}^{corr}(MC)");
  const TString sEffY("#epsilon(p_{T}^{corr})");
  hSumP -> SetLineColor(cSum[0]);
  hSumP -> SetMarkerColor(cSum[0]);
  hSumP -> SetMarkerStyle(mSum[0]);
  hSumP -> SetTitle(sSum.Data());
  hSumP -> SetTitleFont(fTxt);
  hSumP -> GetXaxis() -> SetTitle(sSumX.Data());
  hSumP -> GetXaxis() -> SetTitleFont(fTxt);
  hSumP -> GetXaxis() -> SetTitleOffset(fOffX);
  hSumP -> GetXaxis() -> SetLabelSize(fLbl);
  hSumP -> GetXaxis() -> CenterTitle(fCnt);
  hSumP -> GetYaxis() -> SetTitle(sSumY.Data());
  hSumP -> GetYaxis() -> SetTitleFont(fTxt);
  hSumP -> GetYaxis() -> SetTitleOffset(fOffY);
  hSumP -> GetYaxis() -> SetLabelSize(fLbl);
  hSumP -> GetYaxis() -> CenterTitle(fCnt);
  hSumD -> SetLineColor(cSum[1]);
  hSumD -> SetMarkerColor(cSum[1]);
  hSumD -> SetMarkerStyle(mSum[1]);
  hSumD -> SetTitle(sSum.Data());
  hSumD -> SetTitleFont(fTxt);
  hSumD -> GetXaxis() -> SetTitle(sSumX.Data());
  hSumD -> GetXaxis() -> SetTitleFont(fTxt);
  hSumD -> GetXaxis() -> SetTitleOffset(fOffX);
  hSumD -> GetXaxis() -> SetLabelSize(fLbl);
  hSumD -> GetXaxis() -> CenterTitle(fCnt);
  hSumD -> GetYaxis() -> SetTitle(sSumY.Data());
  hSumD -> GetYaxis() -> SetTitleFont(fTxt);
  hSumD -> GetYaxis() -> SetTitleOffset(fOffY);
  hSumD -> GetYaxis() -> SetLabelSize(fLbl);
  hSumD -> GetYaxis() -> CenterTitle(fCnt);
  hEff  -> SetLineColor(cEff);
  hEff  -> SetMarkerColor(cEff);
  hEff  -> SetMarkerStyle(mEff);
  hEff  -> SetTitle(sEff.Data());
  hEff  -> SetTitleFont(fTxt);
  hEff  -> GetXaxis() -> SetTitle(sEffX.Data());
  hEff  -> GetXaxis() -> SetTitleFont(fTxt);
  hEff  -> GetXaxis() -> SetTitleOffset(fOffX);
  hEff  -> GetXaxis() -> SetLabelSize(fLbl);
  hEff  -> GetXaxis() -> CenterTitle(fCnt);
  hEff  -> GetYaxis() -> SetTitle(sEffY.Data());
  hEff  -> GetYaxis() -> SetTitleFont(fTxt);
  hEff  -> GetYaxis() -> SetTitleOffset(fOffY);
  hEff  -> GetYaxis() -> SetLabelSize(fLbl);
  hEff  -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;


  // make legend
  const UInt_t  cFill(0);
  const UInt_t  cLine(0);
  const UInt_t  fAlign(12);
  const Float_t xyLeg[4]   = {0.1, 0.1, 0.3, 0.3};
  const TString sLabels[2] = {"particle", "detector"};

  TLegend *lSum = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  lSum -> SetFillColor(cFill);
  lSum -> SetLineColor(cLine);
  lSum -> SetTextFont(fTxt);
  lSum -> SetTextAlign(fAlign);
  lSum -> AddEntry(hSumP, sLabels[0].Data());
  lSum -> AddEntry(hSumD, sLabels[1].Data());
  cout << "    Made legends." << endl;


  // draw plot
  const UInt_t  width(1500);
  const UInt_t  height(750);
  const UInt_t  grid(0);
  const UInt_t  log(1);
  const Float_t xyPadEff[4] = {0., 0., 0.5, 1.};
  const Float_t xyPadSum[4] = {0.5, 0., 1., 1.};

  TCanvas *cPlot = new TCanvas("cEff", "", width, height);
  TPad    *pEff  = new TPad("pEff", "", xyPadEff[0], xyPadEff[1], xyPadEff[2], xyPadEff[3]);
  TPad    *pSum  = new TPad("pSum", "", xyPadSum[0], xyPadSum[1], xyPadSum[2], xyPadSum[3]);
  pEff  -> SetGrid(grid, grid);
  pSum  -> SetGrid(grid, grid);
  pSum  -> SetLogy(log);
  cPlot -> cd();
  pEff  -> Draw();
  pSum  -> Draw();
  pEff  -> cd();
  hEff  -> Draw();
  pSum  -> cd();
  hSumP -> Draw();
  hSumD -> Draw("same");
  lSum  -> Draw();
  cPlot -> Write();
  cPlot -> Close();
  cout << "    Drew plot." << endl;


  // close files
  fOut  -> cd();
  hSumP -> Write();
  hSumD -> Write();
  hEff  -> Write();
  fOut  -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fIn[iHist] -> cd();
    fIn[iHist] -> Close();
  }
  cout << "  Finished calculation!\n" << endl;

}

// End ------------------------------------------------------------------------
