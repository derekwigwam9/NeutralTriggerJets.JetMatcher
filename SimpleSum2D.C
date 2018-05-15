// 'SimpleSum2D.C'
// Derek Anderson
// 02.22.2018
//
// Use this to sum a set of histograms with
// a given set of weights from embedding.
//
// NOTE: Configuration = 0, RFF configuration
//       Configuration = 1, FF configuration


#include <iostream>
#include "TH2.h"
#include "TFile.h"
#include "TString.h"

using namespace std;


// constants
static const UInt_t   NHistRFF(8);
static const UInt_t   NHistFF(7);
static const UInt_t   NTotal(10);
static const Double_t NormsRFF[NHistRFF] = {1., 5., 36., 204., 628., 2807., 7547., 3855.};
static const Double_t NormsFF[NHistFF]   = {2., 37., 163., 463., 2419., 8762., 4455.};
static const Double_t WeightsRFF[NTotal] = {1.0, 3.501425e-01, 1.395103e-01, 1.326444e-01, 2.801546e-02, 1.031377e-02, 8.210314e-03, 1.985107e-03, 8.054588e-05, 1.449037e-05};
static const Double_t WeightsFF[NTotal]  = {1.0, 3.361596e-01, 1.401161e-01, 1.337302e-01, 2.895246e-02, 1.042577e-02, 8.294575e-03, 2.064352e-03, 8.088693e-05, 1.417116e-05};

// options
static const UInt_t NHist(7);
static const UInt_t Configuration(1);
static const Bool_t DoIntNorm(false);



void SimpleSum2D() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning sum script..." << endl;

  // io parameters
  const TString  sOut("pp200r9ff.pTcorrResponse.et9vz55had.r03a02rm1chrg.dr03q15.d15m5y2018.root");
  const TString  sIn[NHist]   = {/*"output/pp200r9pt4ff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", */"output/pp200r9pt5ff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt7ff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt9ff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt11ff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt15ff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt25ff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt35ff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root"};
  const TString  sHist[NHist] = {/*"hResponsePtc", */"hResponsePtc", "hResponsePtc", "hResponsePtc", "hResponsePtc", "hResponsePtc", "hResponsePtc", "hResponsePtc"};


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
  TH2D *hHist[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hHist[iHist] = (TH2D*) fIn[iHist] -> Get(sHist[iHist].Data());
    if (!hHist[iHist]) {
      cerr << "PANIC: couldn't grab histogram no. " << iHist << "!" << endl;
      return;
    }
  }
  fOut -> cd();
  cout << "    Grabbed histograms." << endl;


  // scale histograms
  Double_t normer(0.);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    switch (Configuration) {
      case 0:
        hHist[iHist] -> Scale(WeightsRFF[(NTotal - NHist) + iHist]);
        normer += NormsRFF[iHist] * WeightsRFF[(NTotal - NHist) + iHist];
        break;
      case 1:
        hHist[iHist] -> Scale(WeightsFF[(NTotal - NHist) + iHist]);
        normer += NormsFF[iHist] * WeightsFF[(NTotal - NHist) + iHist];
        break;
    }
  }
  cout << "    Scaled histograms." << endl;

  // sum histograms
  const UInt_t  nBinsX = hHist[0] -> GetNbinsX();
  const UInt_t  nBinsY = hHist[0] -> GetNbinsY();
  const Float_t xBin1  = hHist[0] -> GetXaxis() -> GetBinLowEdge(1);
  const Float_t xBin2  = hHist[0] -> GetXaxis() -> GetBinLowEdge(nBinsX + 1);
  const Float_t yBin1  = hHist[0] -> GetYaxis() -> GetBinLowEdge(1);
  const Float_t yBin2  = hHist[0] -> GetYaxis() -> GetBinLowEdge(nBinsY + 1);

  TH2D *hSum  = new TH2D("hSum", "", nBinsX, xBin1, xBin2, nBinsY, yBin1, yBin2);
  TH2D *hNorm = new TH2D("hNorm", "", nBinsX, xBin1, xBin2, nBinsY, yBin1, yBin2);
  hSum  -> Sumw2();
  hNorm -> Sumw2();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hSum  -> Add(hHist[iHist]);
    hNorm -> Add(hHist[iHist]);
  }

  // normalize sum
  if (DoIntNorm) {
    const Double_t integral = hNorm -> Integral();
    hNorm -> Scale(1. / integral);
    cout << "    Normalized histogram.\n"
         << "      normalization = " << integral
         << endl;
  }
  else {
    hNorm -> Scale(1. / normer);
    cout << "    Normalized histogram.\n"
         << "      normalization = " << normer
         << endl;
  }


  // save and close
  fOut  -> cd();
  hSum  -> Write();
  hNorm -> Write();
  fOut  -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fIn[iHist] -> cd();
    fIn[iHist] -> Close();
  }
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
