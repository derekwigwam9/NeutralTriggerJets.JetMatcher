// 'SimpleSum.C'
// Derek Anderson
// 01.17.2018
//
// Use this to sum a set of histograms with
// a given set of weights from embedding.
//
// NOTE: Configuration = 0, RFF configuration
//       Configuration = 1, FF configuration


#include <iostream>
#include "TH1.h"
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
static const UInt_t NHist(8);
static const UInt_t Configuration(0);
static const Bool_t DoIntNorm(false);



void SimpleSum() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning sum script..." << endl;

  // io parameters
  const TString  sOut("pp200r9rff.pTcorrMC.et9vz55had.r03a02rm1chrg.dr03q15.d15m5y2018.root");
  const TString  sIn[NHist]   = {"output/pp200r9pt4rff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt5rff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt7rff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt9rff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt11rff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt15rff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt25rff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root", "output/pp200r9pt35rff.matchNoNorm.et920vz55had.r03a02rm1chrg.dr03q15.root"};
  const TString  sHist[NHist] = {"ParticleJets/hJetPtCorrP", "ParticleJets/hJetPtCorrP", "ParticleJets/hJetPtCorrP", "ParticleJets/hJetPtCorrP", "ParticleJets/hJetPtCorrP", "ParticleJets/hJetPtCorrP", "ParticleJets/hJetPtCorrP", "ParticleJets/hJetPtCorrP"};


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
  TH1D *hHist[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hHist[iHist] = (TH1D*) fIn[iHist] -> Get(sHist[iHist].Data());
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
  const UInt_t  nBins = hHist[0] -> GetNbinsX();
  const Float_t xBin1 = hHist[0] -> GetBinLowEdge(1);
  const Float_t xBin2 = hHist[0] -> GetBinLowEdge(nBins + 1);

  TH1D *hSum  = new TH1D("hSum", "", nBins, xBin1, xBin2);
  TH1D *hNorm = new TH1D("hNorm", "", nBins, xBin1, xBin2);
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
