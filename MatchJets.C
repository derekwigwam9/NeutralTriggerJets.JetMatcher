// 'MatchJets.C'
// Nihar Sahoo, Derek Anderson
// 12.01.2017

#include <vector>
#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TDirectory.h"

using namespace std;


// filepaths
static const TString SOutDefault("test.root");
static const TString SGntDefault("../JetMaker/mc/pp200r9pt35rff.particle.r03rm1chrg.root");
static const TString SMuDefault("../JetMaker/mudst/pp200r9pt35rff.et920vz55had.r03rm1chrg.root");

// jet parameters
static const Double_t MinJetPt = 0.2;
static const Double_t MinArea  = 0.2;  // R03: 0.2, R04: 0.5, R05: 0.65, R07: 1.2
static const Double_t Rcut     = 0.3;  // Rcut = Rjet



void MatchJets(const TString gPath=SGntDefault, const TString uPath=SMuDefault, const TString oPath=SOutDefault, Bool_t inBatchMode=false) {

  cout << "\n  Beginning match script!" << endl;


  // event constants
  const Double_t MaxVz     = 55.;
  const Double_t MaxTrgEta = 0.9;
  const Double_t MinTrgEt  = 9.;
  const Double_t MaxTrgEt  = 20.;
  const Double_t MinTrgTsp = 3.;
  const Double_t MaxTrgTsp = 100.;

  // matching constants
  const Double_t Qmin    = 0.15;  // fraction of jet pT must be above Qmin
  const Double_t Qmax    = 1.5;   // fraction of jet pT must be below Qmax
  const Double_t HardCut = 10.;   // jets w/ pT>HardCut are considered 'hard'
  const Double_t Gcut    = 0.;    // pTgnt must be above this
  const Double_t Ucut    = 0.;    // pTmu must be above this

  // misc. constants
  const Double_t pi  = TMath::Pi();
  const Double_t rDf = pi / 4.;
  cout << "    Opening files and grabbing trees..." << endl;


  // open files
  TFile *fGnt = new TFile(gPath, "read");
  TFile *fMu  = new TFile(uPath, "read");
  TFile *fOut = new TFile(oPath, "recreate");
  if (!fGnt) {
    cerr << "PANIC: Geant input file could not be opened!" << endl;
    assert(fGnt);
  }
  if (!fMu) {
    cerr << "PANIC: MuDst input file could not be opened!" << endl;
    assert(fMu);
  }

  // get trees
  TTree *tGnt;
  TTree *tMu;
  if (fGnt && fMu) {
    fGnt -> GetObject("JetTree", tGnt);
    fMu  -> GetObject("JetTree", tMu);
  }


  // declare Geant event leaves
  cout << "    Setting branch addresses..." << endl;
  Int_t    gEventIndex    = 0;
  Int_t    gRunId         = 0;
  Int_t    gNJets         = 0;
  Double_t gPartonicPt    = 0.;
  Double_t gRefmult       = 0.;
  Double_t gTSP           = 0.;
  Double_t gTrgEta        = 0.;
  Double_t gTrgPhi        = 0.;
  Double_t gTrgEt         = 0.;
  Double_t gRho           = 0.;
  Double_t gSigma         = 0.;
  Double_t gVz            = 0.;
  // declare Geant jet leaves
  vector<Double_t> *gJetPt     = 0;
  vector<Double_t> *gJetNCons  = 0;
  vector<Double_t> *gJetIndex  = 0;
  vector<Double_t> *gJetPtCorr = 0;
  vector<Double_t> *gJetPhi    = 0;
  vector<Double_t> *gJetEta    = 0;
  vector<Double_t> *gJetE      = 0;
  vector<Double_t> *gJetArea   = 0;
  // declare Geant constituent leaves
  vector<vector<Double_t> > *gJetConsPt  = 0;
  vector<vector<Double_t> > *gJetConsEta = 0;
  vector<vector<Double_t> > *gJetConsPhi = 0;
  vector<vector<Double_t> > *gJetConsE   = 0;
  
  // declare MuDst event leaves
  Int_t    uEventIndex    = 0;
  Int_t    uRunId         = 0;
  Int_t    uNJets         = 0;
  Double_t uPartonicPt    = 0.;
  Double_t uRefmult       = 0.;
  Double_t uTSP           = 0.;
  Double_t uTrgEta        = 0.;
  Double_t uTrgPhi        = 0.;
  Double_t uTrgEt         = 0.;
  Double_t uRho           = 0.;
  Double_t uSigma         = 0.;
  Double_t uVz            = 0.;
  // declare MuDst jet leaves
  vector<Double_t> *uJetPt     = 0;
  vector<Double_t> *uJetNCons  = 0;
  vector<Double_t> *uJetIndex  = 0;
  vector<Double_t> *uJetPtCorr = 0;
  vector<Double_t> *uJetEta    = 0;
  vector<Double_t> *uJetPhi    = 0;
  vector<Double_t> *uJetE      = 0;
  vector<Double_t> *uJetArea   = 0;
  // declare MuDst constituent leaves
  vector<vector<Double_t> > *uJetConsPt  = 0;
  vector<vector<Double_t> > *uJetConsEta = 0;
  vector<vector<Double_t> > *uJetConsPhi = 0;
  vector<vector<Double_t> > *uJetConsE   = 0;


  // declare Geant branches
  TBranch *bEventIndexG    = 0;
  TBranch *bRunIdG         = 0;
  TBranch *bNJetsG         = 0;
  TBranch *bRefmultG       = 0;
  TBranch *bPartonicPtG    = 0;
  TBranch *bTspG           = 0;
  TBranch *bTrgEtaG        = 0;
  TBranch *bTrgPhiG        = 0;
  TBranch *bTrgEtG         = 0;
  TBranch *bRhoG           = 0;
  TBranch *bSigmaG         = 0;
  TBranch *bVzG            = 0;
  TBranch *bJetPtG         = 0;
  TBranch *bJetNConsG      = 0;
  TBranch *bJetIndexG      = 0;
  TBranch *bJetPtCorrG     = 0;
  TBranch *bJetEtaG        = 0;
  TBranch *bJetPhiG        = 0;
  TBranch *bJetEG          = 0;
  TBranch *bJetAreaG       = 0;
  TBranch *bJetConsPtG     = 0;
  TBranch *bJetConsEtaG    = 0;
  TBranch *bJetConsPhiG    = 0;
  TBranch *bJetConsEG      = 0;

  // declare MuDst branches
  TBranch *bEventIndexU    = 0;
  TBranch *bRunIdU         = 0;
  TBranch *bNJetsU         = 0;
  TBranch *bPartonicPtU    = 0;
  TBranch *bRefmultU       = 0;
  TBranch *bTspU           = 0;
  TBranch *bTrgEtaU        = 0;
  TBranch *bTrgPhiU        = 0;
  TBranch *bTrgEtU         = 0;
  TBranch *bRhoU           = 0;
  TBranch *bSigmaU         = 0;
  TBranch *bVzU            = 0;
  TBranch *bJetPtU         = 0;
  TBranch *bJetNConsU      = 0;
  TBranch *bJetIndexU      = 0;
  TBranch *bJetPtCorrU     = 0;
  TBranch *bJetEtaU        = 0;
  TBranch *bJetPhiU        = 0;
  TBranch *bJetEU          = 0;
  TBranch *bJetAreaU       = 0;
  TBranch *bJetConsPtU     = 0;
  TBranch *bJetConsEtaU    = 0;
  TBranch *bJetConsPhiU    = 0;
  TBranch *bJetConsEU      = 0;


  // set Geant branches
  tGnt -> SetBranchAddress("eventIndex", &gEventIndex, &bEventIndexG);
  tGnt -> SetBranchAddress("RunId", &gRunId, &bRunIdG);
  tGnt -> SetBranchAddress("Refmult", &gRefmult, &bRefmultG);
  tGnt -> SetBranchAddress("NJets", &gNJets, &bNJetsG);
  tGnt -> SetBranchAddress("PartonicPt", &gPartonicPt, &bPartonicPtG);
  tGnt -> SetBranchAddress("TSP", &gTSP, &bTspG);
  tGnt -> SetBranchAddress("TrgEta", &gTrgEta, &bTrgEtaG);
  tGnt -> SetBranchAddress("TrgPhi", &gTrgPhi, &bTrgPhiG);
  tGnt -> SetBranchAddress("TrgEt", &gTrgEt, &bTrgEtG);
  tGnt -> SetBranchAddress("Rho", &gRho, &bRhoG);
  tGnt -> SetBranchAddress("Sigma", &gSigma, &bSigmaG);
  tGnt -> SetBranchAddress("Vz", &gVz,&bVzG);
  tGnt -> SetBranchAddress("JetIndex", &gJetIndex, &bJetIndexG);
  tGnt -> SetBranchAddress("JetPt", &gJetPt, &bJetPtG);
  tGnt -> SetBranchAddress("JetNCons", &gJetNCons, &bJetNConsG);
  tGnt -> SetBranchAddress("JetPtCorr", &gJetPtCorr, &bJetPtCorrG);
  tGnt -> SetBranchAddress("JetEta", &gJetEta, &bJetEtaG);
  tGnt -> SetBranchAddress("JetPhi",&gJetPhi, &bJetPhiG); 
  tGnt -> SetBranchAddress("JetE", &gJetE, &bJetEG); 
  tGnt -> SetBranchAddress("JetArea",&gJetArea, &bJetAreaG);
  tGnt -> SetBranchAddress("JetConsPt", &gJetConsPt, &bJetConsPtG);
  tGnt -> SetBranchAddress("JetConsEta", &gJetConsEta, &bJetConsEtaG);
  tGnt -> SetBranchAddress("JetConsPhi", &gJetConsPhi, &bJetConsPhiG);
  tGnt -> SetBranchAddress("JetConsE", &gJetConsE, &bJetConsEG);

  // set MuDst branches
  tMu -> SetBranchAddress("eventIndex", &uEventIndex, &bEventIndexU);
  tMu -> SetBranchAddress("RunId", &uRunId, &bRunIdU);
  tMu -> SetBranchAddress("Refmult", &uRefmult, &bRefmultU);
  tMu -> SetBranchAddress("NJets", &uNJets, &bNJetsU);
  tMu -> SetBranchAddress("PartonicPt", &uPartonicPt, &bPartonicPtU);
  tMu -> SetBranchAddress("TSP", &uTSP, &bTspU);
  tMu -> SetBranchAddress("TrgEta", &uTrgEta, &bTrgEtaU);
  tMu -> SetBranchAddress("TrgPhi", &uTrgPhi, &bTrgPhiU);
  tMu -> SetBranchAddress("TrgEt", &uTrgEt, &bTrgEtU);
  tMu -> SetBranchAddress("Rho", &uRho, &bRhoU);
  tMu -> SetBranchAddress("Sigma", &uSigma, &bSigmaU);
  tMu -> SetBranchAddress("Vz", &uVz, &bVzU);
  tMu -> SetBranchAddress("JetIndex", &uJetIndex, &bJetIndexU);
  tMu -> SetBranchAddress("JetPt", &uJetPt, &bJetPtU);
  tMu -> SetBranchAddress("JetNCons", &uJetNCons, &bJetNConsU);
  tMu -> SetBranchAddress("JetPtCorr", &uJetPtCorr, &bJetPtCorrU);
  tMu -> SetBranchAddress("JetEta", &uJetEta, &bJetEtaU);
  tMu -> SetBranchAddress("JetPhi",&uJetPhi, &bJetPhiU); 
  tMu -> SetBranchAddress("JetE", &uJetE, &bJetEU); 
  tMu -> SetBranchAddress("JetArea",&uJetArea, &bJetAreaU);
  tMu -> SetBranchAddress("JetConsPt", &uJetConsPt, &bJetConsPtU);
  tMu -> SetBranchAddress("JetConsEta", &uJetConsEta, &bJetConsEtaU);
  tMu -> SetBranchAddress("JetConsPhi", &uJetConsPhi, &bJetConsPhiU);
  tMu -> SetBranchAddress("JetConsE", &uJetConsE, &bJetConsEU);



  cout << "    Creating histograms..." << endl;

  const Int_t nJetTypes   = 7;
  const Int_t nMatchTypes = 5;
  TH1D *hEfficiencyA;
  TH1D *hEfficiencyPt;
  TH2D *hResponseA;
  TH2D *hResponseAn;
  TH2D *hResponsePt;
  TH2D *hResponsePtN;
  TH2D *hResponsePtc;
  TH2D *hResponsePtcN;
  // event histograms
  TH1D *hRefmultG;
  TH1D *hRefmultU;
  TH1D *hNumJetsG;
  TH1D *hNumJetsU;
  TH1D *hNumToMatch;
  TH1D *hNumMatched;
  TH1D *hNumHardG;
  TH1D *hNumHardU;
  TH1D *hGeantArea;
  TH1D *hMatchArea;
  TH1D *hGeantPtCorr;
  TH1D *hMatchPtCorr;
  // jet histograms  [0='G',1='U',2='C',3='M',4='J',5='Y',6='N']
  TH1D *hJetArea[nJetTypes];
  TH1D *hJetEta[nJetTypes];
  TH1D *hJetPhi[nJetTypes];
  TH1D *hJetPt[nJetTypes];
  TH1D *hJetPtCorr[nJetTypes];
  TH2D *hJetPhiVsEta[nJetTypes];
  // matching histograms [0='U',1='C',2='M',3='J',4='Y']
  TH1D *hJetQt[nMatchTypes];
  TH1D *hJetDr[nMatchTypes];
  TH1D *hJetS[nMatchTypes];
  TH1D *hJetDp[nMatchTypes];
  TH2D *hJetQtVsDr[nMatchTypes];
  TH2D *hJetSvsDr[nMatchTypes];

  // histogram parameters
  const Int_t    nM = 200;
  const Int_t    nN = 100;
  const Int_t    nA = 500;
  const Int_t    nH = 100;
  const Int_t    nF = 360;
  const Int_t    nP = 110;
  const Int_t    nQ = 600;
  const Int_t    nR = 600;
  const Int_t    nS = 600;
  const Int_t    nD = 100;
  const Double_t m1 = 0.;
  const Double_t m2 = 200.;
  const Double_t n1 = 0.;
  const Double_t n2 = 100;
  const Double_t a1 = 0.;
  const Double_t a2 = 5.;
  const Double_t h1 = -5.;
  const Double_t h2 = 5.;
  const Double_t f1 = -2.*pi;
  const Double_t f2 = 2.*pi;
  const Double_t p1 = -10.;
  const Double_t p2 = 100.;
  const Double_t q1 = 0.;
  const Double_t q2 = 3.;
  const Double_t r1 = 0.;
  const Double_t r2 = 3.;
  const Double_t s1 = 0.;
  const Double_t s2 = 3.;
  const Double_t d1 = -50.;
  const Double_t d2 = 50.;
  hEfficiencyA    = new TH1D("hEfficiencyA", "Efficiency, #epsilon(A_{jet}) = N_{match}(A_{jet})/N_{geant}(A_{jet})", nA, a1, a2);
  hEfficiencyPt   = new TH1D("hEfficiencyPt", "Efficiency, #epsilon(p_{T}^{jet}) = N_{match}(p_{T}^{jet})/N_{geant}(p_{T}^{jet})", nP, p1, p2);
  hResponseA      = new TH2D("hResponseA", "Response matrix, jet area; match; geant", nA, a1, a2, nA, a1, a2);
  hResponseAn     = new TH2D("hResponseAn", "Response matrix, jet area (normalized); match; geant", nA, a1, a2, nA, a1, a2);
  hResponsePt     = new TH2D("hResponsePt", "Response matrix, jet p_{T}; match; geant", nP, p1, p2, nP, p1, p2);
  hResponsePtN    = new TH2D("hResponsePtN", "Response matrix, jet p_{T} (normalized); match; geant", nP, p1, p2, nP, p1, p2);
  hResponsePtc    = new TH2D("hResponsePtc", "Response matrix, jet p_{T}^{corr}; match; geant", nP, p1, p2, nP, p1, p2);
  hResponsePtcN   = new TH2D("hResponsePtcN", "Response matrix, jet p_{T}^{corr} (normalized); match; geant", nP, p1, p2, nP, p1, p2);
  // event histograms
  hRefmultG       = new TH1D("hRefmultG", "Geant Refmult", nM, m1, m2);
  hRefmultU       = new TH1D("hRefmultU", "MuDst Refmult", nM, m1, m2);
  hNumJetsG       = new TH1D("hNumJetsG", "no. of jets, Geant", nN, n1, n2);
  hNumJetsU       = new TH1D("hNumJetsU", "No. of jets, MuDst", nN, n1, n2);
  hNumToMatch     = new TH1D("hNumToMatch", "No. of jets to match (ie. no. of jets w/ pT > pTcut)", nN, n1, n2);
  hNumMatched     = new TH1D("hNumMatched", "No. of jets matched", nN, n1, n2);
  hNumHardG       = new TH1D("hNumHardG", "No. of jets w/ p_{T} above a threshold, Geant", nN, n1, n2);
  hNumHardU       = new TH1D("hNumHardU", "No. of jets w/ p_{T} above a threshold, MuDst", nN, n1, n2);
  hGeantArea      = new TH1D("hGeantArea", "Total no. of jets to match per A_{jet} bin (for efficiency)", nA, a1, a2);
  hMatchArea      = new TH1D("hMatchArea", "Total no. of jets matched per A_{jet} bin (for efficiency)", nA, a1, a2);
  hGeantPtCorr    = new TH1D("hGeantPtCorr", "Total no. of jets to match per p_{T}^{corr} bin (for efficiency)", nP, p1, p2);
  hMatchPtCorr    = new TH1D("hMatchPtCorr", "Total no. of jets matched per p_{T}^{corr} bin (for efficiency)", nP, p1, p2);
  // geant jets
  hJetArea[0]     = new TH1D("hJetAreaG", "Jet area, Geant", nA, a1, a2);
  hJetEta[0]      = new TH1D("hJetEtaG", "Jet eta, Geant", nH, h1, h2);
  hJetPhi[0]      = new TH1D("hJetPhiG", "Jet phi, Geant", nF, f1, f2);
  hJetPt[0]       = new TH1D("hJetPtG", "Jet p_{T}, Geant", nP, p1, p2);
  hJetPtCorr[0]   = new TH1D("hJetPtCorrG", "Jet p_{T}^{corr}, Geant", nP, p1, p2);
  hJetPhiVsEta[0] = new TH2D("hJetPhiVsEtaG", "Jet #varphi vs. #eta, Geant", nH, h1, h2, nF, f1, f2);
  // mudst jets
  hJetArea[1]     = new TH1D("hJetAreaU", "Jet area, MuDst", nA, a1, a2);
  hJetEta[1]      = new TH1D("hJetEtaU", "Jet eta, MuDst", nH, h1, h2);
  hJetPhi[1]      = new TH1D("hJetPhiU", "Jet phi, MuDst", nF, f1, f2);
  hJetPt[1]       = new TH1D("hJetPtU", "Jet p_{T}, MuDst", nP, p1, p2);
  hJetPtCorr[1]   = new TH1D("hJetPtCorrU", "Jet p_{T}^{corr}, MuDst", nP, p1, p2);
  hJetPhiVsEta[1] = new TH2D("hJetPhiVsEtaU", "Jet #varphi vs. #eta, MuDst", nH, h1, h2, nF, f1, f2);
  hJetQt[0]       = new TH1D("hJetQtU", "Jet q_{T}, MuDst (normalization different!)", nQ, q1, q2);
  hJetDr[0]       = new TH1D("hJetDrU", "Jet #Deltar, MuDst (normalization different!)", nR, r1, r2);
  hJetS[0]        = new TH1D("hJetSu", "Jet s=A_{#mu}/A_{geant}, MuDst (normalization different!)", nS, s1, s2);
  hJetDp[0]       = new TH1D("hJetDpU", "Jet #Deltap_{T}=p_{T}^{#mu}-p_{T}^{geant}, MuDst (normalization different!)", nD, d1, d2);
  hJetQtVsDr[0]   = new TH2D("hJetQtVsDrU", "Jet q_{T} vs. #Deltar, MuDst (normalization different!); #Deltar; q_{T}", nR, r1, r2, nQ, q1, q2);
  hJetSvsDr[0]    = new TH2D("hJetSvsDrU", "Jet s vs. #Deltar, MuDst (normalization different!); #Deltar; s", nR, r1, r2, nS, s1, s2);
  // candidate matches
  hJetArea[2]     = new TH1D("hJetAreaC", "Jet area, candidates", nA, a1, a2);
  hJetEta[2]      = new TH1D("hJetEtaC", "Jet eta, candidates", nH, h1, h2);
  hJetPhi[2]      = new TH1D("hJetPhiC", "Jet phi, candidates", nF, f1, f2);
  hJetPt[2]       = new TH1D("hJetPtC", "Jet p_{T}, candidates", nP, p1, p2);
  hJetPtCorr[2]   = new TH1D("hJetPtCorrC", "Jet p_{T}^{corr}, candidates", nP, p1, p2);
  hJetPhiVsEta[2] = new TH2D("hJetPhiVsEtaC", "Jet #varphi vs. #eta, candidates", nH, h1, h2, nF, f1, f2);
  hJetQt[1]       = new TH1D("hJetQtC", "Jet q_{T}, candidates", nQ, q1, q2);
  hJetDr[1]       = new TH1D("hJetDrC", "Jet #Deltar, candidates", nR, r1, r2);
  hJetS[1]        = new TH1D("hJetSc", "Jet s=A_{cand.}/A_{geant}, candidates", nS, s1, s2);
  hJetDp[1]       = new TH1D("hJetDpC", "Jet #Deltap_{T}=p_{T}^{cand.}-p_{T}^{geant}, candidates (normalization different!)", nD, d1, d2);
  hJetQtVsDr[1]   = new TH2D("hJetQtVsDrC", "Jet q_{T} vs. #Deltar, candidates; #Deltar; q_{T}", nR, r1, r2, nQ, q1, q2);
  hJetSvsDr[1]    = new TH2D("hJetSvsDrC", "Jet s vs. #Deltar, candidates; #Deltar; s", nR, r1, r2, nS, s1, s2);
  // matches
  hJetArea[3]     = new TH1D("hJetAreaM", "Jet area, matches", nA, a1, a2);
  hJetEta[3]      = new TH1D("hJetEtaM", "Jet eta, matches", nH, h1, h2);
  hJetPhi[3]      = new TH1D("hJetPhiM", "Jet phi, matches", nF, f1, f2);
  hJetPt[3]       = new TH1D("hJetPtM", "Jet p_{T}, matches", nP, p1, p2);
  hJetPtCorr[3]   = new TH1D("hJetPtCorrM", "Jet p_{T}^{corr}, matches", nP, p1, p2);
  hJetPhiVsEta[3] = new TH2D("hJetPhiVsEtaM", "Jet #varphi vs. #eta, matches", nH, h1, h2, nF, f1, f2);
  hJetQt[2]       = new TH1D("hJetQtM", "Jet q_{T}, matches", nQ, q1, q2);
  hJetDr[2]       = new TH1D("hJetDrM", "Jet #Deltar, matches", nR, r1, r2);
  hJetS[2]        = new TH1D("hJetSm", "Jet s=A_{match}/A_{geant}, matches", nS, s1, s2);
  hJetDp[2]       = new TH1D("hJetDpM", "Jet #Deltap_{T}=p_{T}^{match}-p_{T}^{geant}, matches (normalization different!)", nD, d1, d2);
  hJetQtVsDr[2]   = new TH2D("hJetQtVsDrM", "Jet q_{T} vs. #Deltar, matches", nR, r1, r2, nQ, q1, q2);
  hJetSvsDr[2]    = new TH2D("hJetSvsDrM", "Jet s vs. #Deltar, matches", nR, r1, r2, nS, s1, s2);
  // junk (mudst jets that weren't matched)
  hJetArea[4]     = new TH1D("hJetAreaJ", "Jet area, junk", nA, a1, a2);
  hJetEta[4]      = new TH1D("hJetEtaJ", "Jet eta, junk", nH, h1, h2);
  hJetPhi[4]      = new TH1D("hJetPhiJ", "Jet phi, junk", nF, f1, f2);
  hJetPt[4]       = new TH1D("hJetPtJ", "Jet p_{T}, junk", nP, p1, p2);
  hJetPtCorr[4]   = new TH1D("hJetPtCorrJ", "Jet p_{T}^{corr}, junk", nP, p1, p2);
  hJetPhiVsEta[4] = new TH2D("hJetPhiVsEtaJ", "Jet #varphi vs. #eta, junk", nH, h1, h2, nF, f1, f2);
  hJetQt[3]       = new TH1D("hJetQtJ", "Jet q_{T}, junk (normalization different!)", nQ, q1, q2);
  hJetDr[3]       = new TH1D("hJetDrJ", "Jet #Deltar, junk (normalization different!)", nR, r1, r2);
  hJetS[3]        = new TH1D("hJetSj", "Jet s=A_{junk}/A_{geant}, junk (normalization different!)", nS, s1, s2);
  hJetDp[3]       = new TH1D("hJetDpJ", "Jet #Deltap_{T}=p_{T}^{junk}-p_{T}^{geant}, junk (normalization different!)", nD, d1, d2);
  hJetQtVsDr[3]   = new TH2D("hJetQtVsDrJ", "Jet q_{T} vs. #Deltar, junk (normalization different!)", nR, r1, r2, nQ, q1, q2);
  hJetSvsDr[3]    = new TH2D("hJetSvsDrJ", "Jet s vs. #Deltar, junk (normalization different!)", nR, r1, r2, nS, s1, s2);
  // mystery (jets w/ dR > Rjet and |qT-1|<.1)
  hJetArea[5]     = new TH1D("hJetAreaY", "Jet area, mystery", nA, a1, a2);
  hJetEta[5]      = new TH1D("hJetEtaY", "Jet eta, mystery", nH, h1, h2);
  hJetPhi[5]      = new TH1D("hJetPhiY", "Jet phi, mystery", nF, f1, f2);
  hJetPt[5]       = new TH1D("hJetPtY", "Jet p_{T}, mystery", nP, p1, p2);
  hJetPtCorr[5]   = new TH1D("hJetPtCorrY", "Jet p_{T}^{corr}, mystery", nP, p1, p2);
  hJetPhiVsEta[5] = new TH2D("hJetPhiVsEtaY", "Jet #varphi vs. #eta, mystery", nH, h1, h2, nF, f1, f2);
  hJetQt[4]       = new TH1D("hJetQtY", "Jet q_{T}, mystery", nQ, q1, q2);
  hJetDr[4]       = new TH1D("hJetDrY", "Jet #Deltar, mystery", nR, r1, r2);
  hJetS[4]        = new TH1D("hJetSy", "Jet s=A_{?}/A_{geant}, mystery", nS, s1, s2);
  hJetDp[4]       = new TH1D("hJetDpY", "Jet #Deltap_{T}=p_{T}^{?}-p_{T}^{geant}, mystery (normalization different!)", nD, d1, d2);
  hJetQtVsDr[4]   = new TH2D("hJetQtVsDrY", "Jet q_{T} vs. #Deltar, mystery", nR, r1, r2, nQ, q1, q2);
  hJetSvsDr[4]    = new TH2D("hJetSvsDrY", "Jet s vs. #Deltar, mystery", nR, r1, r2, nS, s1, s2);
  // not matches (geant jets that weren't matched)
  hJetArea[6]     = new TH1D("hJetAreaN", "Jet area, (geant) not matches", nA, a1, a2);
  hJetEta[6]      = new TH1D("hJetEtaN", "Jet eta, (geant) not matches", nH, h1, h2);
  hJetPhi[6]      = new TH1D("hJetPhiN", "Jet phi, (geant) not matches", nF, f1, f2);
  hJetPt[6]       = new TH1D("hJetPtN", "Jet p_{T}, (geant) not matches", nP, p1, p2);
  hJetPtCorr[6]   = new TH1D("hJetPtCorrN", "Jet p_{T}^{corr}, (geant) not matches", nP, p1, p2);
  hJetPhiVsEta[6] = new TH2D("hJetPhiVsEtaN", "Jet #varphi vs. #eta, (geant) not matches", nH, h1, h2, nF, f1, f2);

  // errors
  hEfficiencyA   -> Sumw2();
  hEfficiencyPt  -> Sumw2();
  hResponseA     -> Sumw2();
  hResponseAn    -> Sumw2();
  hResponsePt    -> Sumw2();
  hResponsePtN   -> Sumw2();
  hResponsePtc   -> Sumw2();
  hResponsePtcN  -> Sumw2();
  hRefmultG      -> Sumw2();
  hRefmultU      -> Sumw2();
  hNumJetsG      -> Sumw2();
  hNumJetsU      -> Sumw2();
  hNumToMatch    -> Sumw2();
  hNumMatched    -> Sumw2();
  hNumHardG      -> Sumw2();
  hNumHardU      -> Sumw2();
  hGeantArea     -> Sumw2();
  hMatchArea     -> Sumw2();
  hGeantPtCorr   -> Sumw2();
  hMatchPtCorr   -> Sumw2();
  for (Int_t i = 0; i < nJetTypes; i++) {
    hJetArea[i]   -> Sumw2();
    hJetEta[i]    -> Sumw2();
    hJetPhi[i]    -> Sumw2();
    hJetPt[i]     -> Sumw2();
    hJetPtCorr[i] -> Sumw2();
  }
  for (Int_t i = 0; i < nMatchTypes; i++) {
    hJetQt[i] -> Sumw2();
    hJetDr[i] -> Sumw2();
    hJetS[i]  -> Sumw2();
    hJetDp[i] -> Sumw2();
  }


  // check to make sure there are a reasonable no. of events
  Int_t fEvt  = 0;
  Int_t nEvts = 0;
  Int_t gEvts = (Int_t) tGnt -> GetEntries();
  Int_t uEvts = (Int_t) tMu  -> GetEntries();
  if (gEvts < uEvts) {
    cerr << "WARNING: There are less particle-level events than detector-level!\n"
         << "         Please double-check that everything is in order...\n"
         << "         nParticle = " << gEvts << ", nDetector = " << uEvts
         << endl;
    if (gEvts > uEvts) fEvt = 1;
    if (gEvts < uEvts) fEvt = 2;
    nEvts = TMath::Min(gEvts, uEvts);
  }
  else {
    fEvt  = 1;
    nEvts = uEvts;
  }

  // create map of tree-index to event / run no.
  Int_t gMap[gEvts][3];
  Int_t uMap[uEvts][3];
  for (Int_t iGnt = 0; iGnt < gEvts; iGnt++) {
    tGnt -> GetEntry(iGnt);
    gMap[iGnt][0] = iGnt;
    gMap[iGnt][1] = gEventIndex;
    gMap[iGnt][2] = gRunId;
  }
  for (Int_t iDst = 0; iDst < uEvts; iDst++) {
    tMu -> GetEntry(iDst);
    uMap[iDst][0] = iDst;
    uMap[iDst][1] = uEventIndex;
    uMap[iDst][2] = uRunId;
  }


  // vector for matching
  vector<Int_t> matchIndices;
  cout << "    Beginning event loop..." << endl;

  // event loop
  Int_t nShift   = 0;
  Int_t nFound   = 0;
  Int_t nTrig    = 0;
  Int_t nBytesG  = 0;
  Int_t nBytesU  = 0;
  Int_t breakVal = 0;
  for (Int_t i = 0; i < nEvts; i++) {

    // locate event in geant or mudst tree
    Int_t iGntTree = -1;
    Int_t iDstTree = -1;
    switch (fEvt) {
      case 1:
        iDstTree = i;
        tMu -> GetEntry(iDstTree);
        for (Int_t iGnt = 0; iGnt < gEvts; iGnt++) {
          const Bool_t evtMatch = (gMap[iGnt][1] == uEventIndex);
          const Bool_t runMatch = (gMap[iGnt][2] == uRunId);
          if (evtMatch && runMatch) {
            iGntTree = gMap[iGnt][0];
            break;
          }
        }
        break;
      case 2:
        iGntTree = i;
        tGnt -> GetEntry(iGntTree);
        for (Int_t iDst = 0; iDst < uEvts; iDst++) {
          const Bool_t evtMatch = (uMap[iDst][1] == gEventIndex);
          const Bool_t runMatch = (uMap[iDst][2] == gRunId);
          if (evtMatch && runMatch) {
            iDstTree = uMap[iDst][0];
            break;
          }
        }
        break;
    }  // end swich case
    const Bool_t didNotFindGnt = (iGntTree == -1);
    const Bool_t didNotFindDst = (iDstTree == -1);
    if (didNotFindGnt || didNotFindDst) continue;

    // load entries
    Int_t gBytes = tGnt -> GetEntry(iGntTree);
    Int_t uBytes = tMu  -> GetEntry(iDstTree);
    if (gBytes < 0) {
      cerr << "ERROR: problem with Geant event " << i << "...\n" << endl;
      breakVal = 1;
      break;
    }
    if (uBytes < 0) {
      cerr << "ERROR: problem with MuDst event " << i + nShift << "..." << endl;
      breakVal = 1;
      break;
    }

    // should be same run and event
    Bool_t isSameEvent = (gEventIndex == uEventIndex);
    Bool_t isSameRun   = (gRunId == uRunId);
    if (!isSameEvent || !isSameRun) {
      cerr << "PANIC: event index and run ID are NOT the same! Stopped at i = " << i << "\n"
           << "       GeantEvt = " << gEventIndex << ", MuDstEvt = " << uEventIndex << "\n"
           << "       GeantRun = " << gRunId << ", MuDstRun = " << uRunId
           << endl;
      breakVal = 1;
      break;
    }
    nFound++;

    nBytesG += gBytes;
    nBytesU += uBytes;
    if (!inBatchMode) {
      cout << "      Processing event " << i+1 << "/" << nEvts << "...\r" << flush;
      if (i+1 == nEvts) cout << endl;
    }
    else 
      cout << "      Processing event " << i+1 << "/" << nEvts << "..." << endl;


    // trigger info
    const Double_t vZtrg  = uVz;
    const Double_t eTtrg  = uTrgEt;
    const Double_t hTrg   = uTrgEta;
    const Double_t tspTrg = TMath::Abs(uTSP);

    // trigger cuts
    const Bool_t isInVzCut  = (TMath::Abs(vZtrg) < MaxVz);
    const Bool_t isInEtaCut = (TMath::Abs(hTrg) < MaxTrgEta);
    const Bool_t isInEtCut  = ((eTtrg > MinTrgEt) && (eTtrg < MaxTrgEt));
    const Bool_t isInTspCut = ((tspTrg >= MinTrgTsp) && (tspTrg <= MaxTrgTsp));
    if (!isInVzCut || !isInEtaCut || !isInEtCut || !isInTspCut) continue;
    nTrig++;


    // Geant jet loop
    Int_t nHardG   = 0;
    Int_t nToMatch = 0;
    Int_t nMatched = 0;
    Int_t nGjets   = (Int_t) gJetEta -> size();
    Int_t nUjets   = (Int_t) uJetEta -> size();
    for (Int_t j = 0; j < nGjets; j++) {

      const Double_t gA   = gJetArea   -> at(j);
      const Double_t gH   = gJetEta    -> at(j);
      const Double_t gF   = gJetPhi    -> at(j);
      const Double_t gPt  = gJetPt     -> at(j);
      //const Double_t gPtc = gJetPtCorr -> at(j);
      const Double_t gPtc = gPt - (gRho * gA);  // quick fix [11.27.2017] 
      if (gPt > HardCut)
        ++nHardG;

      // calculate delta phi
      Double_t gDf = gF - gTrgPhi;
      if (gDf < ((-1. * pi) / 2.)) gDf += (2. * pi);
      if (gDf > ((3. * pi) / 2.))  gDf -= (2. * pi);
      const Double_t gDfCut    = TMath::Abs(gDf - pi);
      const Bool_t   isRecoilG = (gDfCut < rDf);

      if (gPt < MinJetPt)
        continue;
      if (gA < MinArea)
        continue;
      if (!isRecoilG)
        continue;

      // fill Geant histograms
      hJetArea[0]     -> Fill(gA);
      hJetEta[0]      -> Fill(gH);
      hJetPhi[0]      -> Fill(gF);
      hJetPt[0]       -> Fill(gPt);
      hJetPtCorr[0]   -> Fill(gPtc);
      hJetPhiVsEta[0] -> Fill(gH, gF);


      if (gPt < Gcut)
        continue;
      else {
        hGeantArea   -> Fill(gA);
        hGeantPtCorr -> Fill(gPtc);
        ++nToMatch;
      }

      // match jets ['b' for best]
      Int_t    bIndex = 0.;
      Double_t bH     = 0.;
      Double_t bPt    = 0.;
      Double_t bPtc   = 0.;
      Double_t bF     = 0.;
      Double_t bA     = 0.;
      Double_t bQt    = 0.;
      Double_t bS     = 0.;
      Double_t bDp    = 0.;
      Double_t bDr    = 999.;

      // MuDst jet loop
      Bool_t isMatched = false;
      for (Int_t k = 0; k < nUjets; k++) {

        const Double_t uA   = uJetArea   -> at(k);
        const Double_t uH   = uJetEta    -> at(k);
        const Double_t uF   = uJetPhi    -> at(k);
        const Double_t uPt  = uJetPt     -> at(k);
        //const Double_t uPtc = uJetPtCorr -> at(k);
        const Double_t uPtc = uPt - (uRho * uA);  // quick fix [11.27.2017]

        // calculate delta phi
        Double_t uDf = uF - uTrgPhi;
        if (uDf < ((-1. * pi) / 2.)) uDf += (2. * pi);
        if (uDf > ((3. * pi) / 2.))  uDf -= (2. * pi);
        const Double_t uDfCut    = TMath::Abs(uDf - pi);
        const Bool_t   isRecoilU = (uDfCut < rDf);

        Bool_t isInAcceptance = true;
        if ((uPt < MinJetPt) || (uA < MinArea) || (!isRecoilU))
          isInAcceptance = false;


        // match jets
        Double_t qT = uPt / gPt;
        Double_t s  = uA / gA;
        Double_t dP = uPt - gPt;
        Double_t dH = uH - gH;
        Double_t dF = uF - gF;
        Double_t dR = sqrt(dH*dH + dF*dF);

        if (uPt < Ucut)
          continue;

        // fill MuDst histograms
        hJetQt[0]     -> Fill(qT);
        hJetDr[0]     -> Fill(dR);
        hJetS[0]      -> Fill(s);
        hJetDp[0]     -> Fill(dP);
        hJetQtVsDr[0] -> Fill(dR, qT);
        hJetSvsDr[0]  -> Fill(dR, s);


        Bool_t isBetter = false;
        Bool_t isInRcut = (dR < Rcut);
        Bool_t isInQcut = ((qT > Qmin) && (qT < Qmax));
        if (isInRcut && isInQcut && isInAcceptance) {
          isMatched = true;
          isBetter  = ((dR < bDr) && (qT > bQt));

          // fill candidate histograms
          hJetArea[2]     -> Fill(uA);
          hJetEta[2]      -> Fill(uH);
          hJetPhi[2]      -> Fill(uF);
          hJetPt[2]       -> Fill(uPt);
          hJetPtCorr[2]   -> Fill(uPtc);
          hJetPhiVsEta[2] -> Fill(uH, uF);
          hJetQt[1]       -> Fill(qT);
          hJetDr[1]       -> Fill(dR);
          hJetS[1]        -> Fill(s);
          hJetDp[1]       -> Fill(dP);
          hJetQtVsDr[1]   -> Fill(dR, qT);
          hJetSvsDr[1]    -> Fill(dR, s);
        }
        else {
          // fill junk histograms
          hJetQt[3]     -> Fill(qT);
          hJetDr[3]     -> Fill(dR);
          hJetS[3]      -> Fill(s);
          hJetDp[3]     -> Fill(dP);
          hJetQtVsDr[3] -> Fill(dR, qT);
          hJetSvsDr[3]  -> Fill(dR, s);
        }

        // fill mystery histograms
        Double_t qCut      = TMath::Abs(qT - 1);
        Bool_t   isNearOne = (qCut < 0.1);
        if (!isInRcut && isNearOne) {
          hJetArea[5]     -> Fill(uA);
          hJetEta[5]      -> Fill(uH);
          hJetPhi[5]      -> Fill(uF);
          hJetPt[5]       -> Fill(uPt);
          hJetPtCorr[5]   -> Fill(uPtc);
          hJetPhiVsEta[5] -> Fill(uH, uF);
          hJetQt[4]       -> Fill(qT);
          hJetDr[4]       -> Fill(dR);
          hJetS[4]        -> Fill(s);
          hJetDp[4]       -> Fill(dP);
          hJetQtVsDr[4]   -> Fill(dR, qT);
          hJetSvsDr[4]    -> Fill(dR, s);
        }

        // check if candidate is best match
        if (isMatched && isBetter) {
          bIndex = k;
          bA     = uA;
          bH     = uH;
          bF     = uF;
          bPt    = uPt;
          bPtc   = uPtc;
          bQt    = qT;
          bS     = s;
          bDp    = dP;
          bDr    = dR;
        }

      }  // end MuDst jet loop


      // fill match histograms
      if (isMatched) {
        hResponseA      -> Fill(bA, gA);
        hResponseAn     -> Fill(bA, gA);
        hResponsePt     -> Fill(bPt, gPt);
        hResponsePtN    -> Fill(bPt, gPt);
        hResponsePtc    -> Fill(bPtc, gPtc);
        hResponsePtcN   -> Fill(bPtc, gPtc);
        hMatchArea      -> Fill(bA);
        hMatchPtCorr    -> Fill(bPtc);
        hJetArea[3]     -> Fill(bA);
        hJetEta[3]      -> Fill(bH);
        hJetPhi[3]      -> Fill(bF);
        hJetPt[3]       -> Fill(bPt);
        hJetPtCorr[3]   -> Fill(bPtc);
        hJetPhiVsEta[3] -> Fill(bH, bF);
        hJetQt[2]       -> Fill(bQt);
        hJetDr[2]       -> Fill(bDr);
        hJetS[2]        -> Fill(bS);
        hJetDp[2]       -> Fill(bDp);
        hJetQtVsDr[2]   -> Fill(bDr, bQt);
        hJetSvsDr[2]    -> Fill(bDr, bS);
        matchIndices.push_back(bIndex);
        ++nMatched;
      }
      else {
        hJetArea[6]     -> Fill(gA);
        hJetEta[6]      -> Fill(gH);
        hJetPhi[6]      -> Fill(gF);
        hJetPt[6]       -> Fill(gPt);
        hJetPtCorr[6]   -> Fill(gPtc);
        hJetPhiVsEta[6] -> Fill(gH, gF);
      }

    }  // end Geant jet loop

    Int_t matchSize = (Int_t) matchIndices.size();
    if (nMatched != matchSize) {
      cerr << "ERROR: matchIndices did something weird in event " << i << "..." << endl;
      breakVal = 1;
      break;
    }


    // MuDst (detector) jet loop
    Int_t nHardU = 0;
    for (Int_t j = 0; j < nUjets; j++) {

      Double_t uA   = uJetArea   -> at(j);
      Double_t uH   = uJetEta    -> at(j);
      Double_t uF   = uJetPhi    -> at(j);
      Double_t uPt  = uJetPt     -> at(j);
      //Double_t uPtc = uJetPtCorr -> at(j);
      Double_t uPtc = uPt - (uRho * uA);  // quick fix [11.27.2017]
      if (uPt > HardCut)
        ++nHardU;

      // calculate delta phi
      Double_t uDf = uF - uTrgPhi;
      if (uDf < ((-1. * pi) / 2.)) uDf += (2. * pi);
      if (uDf > ((3. * pi) / 2.))  uDf -= (2. * pi);
      const Double_t uDfCut    = TMath::Abs(uDf - pi);
      const Bool_t   isRecoilU = (uDfCut < rDf);

      if (uPt < MinJetPt)
        continue;
      if (uA < MinArea)
        continue;
      if (!isRecoilU)
        continue;

      // fill MuDst histograms
      hJetArea[1]     -> Fill(uA);
      hJetEta[1]      -> Fill(uH);
      hJetPhi[1]      -> Fill(uF);
      hJetPt[1]       -> Fill(uPt);
      hJetPtCorr[1]   -> Fill(uPtc);
      hJetPhiVsEta[1] -> Fill(uH, uF);


      // check if MuDst jet matches Geant jet
      Bool_t isMatch = false;
      for (Int_t k = 0; k < nMatched; k++) {
        Int_t m = matchIndices.at(k);
        if (j == m) {
          isMatch = true;
          break;
        }
      }

      // fill junk histograms
      if (!isMatch) {
        hJetArea[4]     -> Fill(uA);
        hJetEta[4]      -> Fill(uH);
        hJetPhi[4]      -> Fill(uF);
        hJetPt[4]       -> Fill(uPt);
        hJetPtCorr[4]   -> Fill(uPtc);
        hJetPhiVsEta[4] -> Fill(uH, uF);
      }

    }  // end MuDst jet loop


    // fill event histograms
    hRefmultG   -> Fill(gRefmult);
    hRefmultU   -> Fill(uRefmult);
    hNumJetsG   -> Fill(gNJets);
    hNumJetsU   -> Fill(uNJets);
    hNumToMatch -> Fill(nToMatch);
    hNumMatched -> Fill(nMatched);
    hNumHardG   -> Fill(nHardG);
    hNumHardU   -> Fill(nHardU);

    matchIndices.clear();

  }  // end event loop


  if (breakVal == 1) {
    cerr << "ERROR: Occured during event loop!\n"
         << "       Aborting program!"
         << endl;
    assert(0);
  }
  else
    cout << "    Event loop finished!\n"
         << "      nFound = " << nFound << "\n"
         << "      nTrig  = " << nTrig
         << endl;


  // calculate efficiency
  cout << "    Calculating efficiency..." << endl;
  TH1D *hGeantA  = (TH1D*) hJetArea[0] -> Clone();
  TH1D *hMatchA  = (TH1D*) hJetArea[3] -> Clone();
  TH1D *hGeantPt = (TH1D*) hJetPt[0]   -> Clone();
  TH1D *hMatchPt = (TH1D*) hJetPt[3]   -> Clone();
  hEfficiencyA  -> Divide(hMatchA, hGeantA, 1., 1.);
  hEfficiencyPt -> Divide(hMatchPt, hGeantPt, 1., 1.);


  cout << "    Normalizing histograms..." << endl;

  // bin widths
  const Double_t mBin   = (m2 - m1) / nM;
  const Double_t nBin   = (n2 - n1) / nN;
  const Double_t aBin   = (a2 - a1) / nA;
  const Double_t hBin   = (h2 - h1) / nH;
  const Double_t fBin   = (f2 - f1) / nF;
  const Double_t pBin   = (p2 - p1) / nP;
  // overall normalizations
  const Double_t mNorm  = nEvts * mBin;
  const Double_t nNorm  = nEvts * nBin;
  const Double_t aNorm  = nEvts * aBin;
  const Double_t hNorm  = nEvts * hBin;
  const Double_t fNorm  = nEvts * fBin;
  const Double_t pNorm  = nEvts * pBin;
  const Double_t hfNorm = hNorm * fBin;
  // normalize histograms
  hRefmultG   -> Scale(1. / mNorm);
  hRefmultU   -> Scale(1. / mNorm);
  hNumJetsG   -> Scale(1. / nNorm);
  hNumJetsU   -> Scale(1. / nNorm);
  hNumToMatch -> Scale(1. / nNorm);
  hNumMatched -> Scale(1. / nNorm);
  hNumHardG   -> Scale(1. / nNorm);
  hNumHardU   -> Scale(1. / nNorm);
  for (Int_t i = 0; i < nJetTypes; i++) {
    hJetArea[i]     -> Scale(1. / aNorm);
    hJetEta[i]      -> Scale(1. / hNorm);
    hJetPhi[i]      -> Scale(1. / fNorm);
    hJetPt[i]       -> Scale(1. / pNorm);
    hJetPtCorr[i]   -> Scale(1. / pNorm);
    hJetPhiVsEta[i] -> Scale(1. / hfNorm);
  }

  // normalize response matrices
  const Int_t nAbinsX = hResponseAn   -> GetNbinsX();
  const Int_t nAbinsY = hResponseAn   -> GetNbinsY();
  const Int_t nPbinsX = hResponsePtN  -> GetNbinsX();
  const Int_t nPbinsY = hResponsePtN  -> GetNbinsY();
  const Int_t nCbinsX = hResponsePtcN -> GetNbinsX();
  const Int_t nCbinsY = hResponsePtcN -> GetNbinsY();
  for (Int_t i = 1; i < nAbinsY+1; i++) {
    const Double_t aNorm = hResponseAn -> Integral(1, nAbinsX, i, i);
    if (aNorm == 0.) continue;

    for (Int_t j = 1; j < nAbinsX+1; j++) {
      const Double_t oldCnt = hResponseAn -> GetBinContent(j, i);
      const Double_t oldErr = hResponseAn -> GetBinError(j, i);
      const Double_t xWidth = hResponseAn -> GetXaxis() -> GetBinWidth(j);
      const Double_t yWidth = hResponseAn -> GetYaxis() -> GetBinWidth(i);
      const Double_t dArea  = xWidth * yWidth;
      const Double_t newCnt = (oldCnt / aNorm) * dArea;
      const Double_t newErr = (oldErr / aNorm) * dArea;
      hResponseAn -> SetBinContent(j, i, newCnt);
      hResponseAn -> SetBinError(j, i, newErr);
    }
  }
  for (Int_t i = 1; i < nPbinsY+1; i++) {
    const Double_t pNorm = hResponsePtN -> Integral(1, nPbinsX, i, i);
    if (pNorm == 0.) continue;

    for (Int_t j = 1; j < nPbinsX+1; j++) {
      const Double_t oldCnt = hResponsePtN -> GetBinContent(j, i);
      const Double_t oldErr = hResponsePtN -> GetBinError(j, i);
      const Double_t xWidth = hResponsePtN -> GetXaxis() -> GetBinWidth(j);
      const Double_t yWidth = hResponsePtN -> GetYaxis() -> GetBinWidth(i);
      const Double_t dArea  = xWidth * yWidth;
      const Double_t newCnt = (oldCnt / pNorm) * dArea;
      const Double_t newErr = (oldErr / pNorm) * dArea;
      hResponsePtN -> SetBinContent(j, i, newCnt);
      hResponsePtN -> SetBinError(j, i, newErr);
    }
  }
  for (Int_t i = 1; i < nCbinsY+1; i++) {
    const Double_t cNorm = hResponsePtcN -> Integral(1, nCbinsX, i, i);
    if (cNorm == 0.) continue;

    for (Int_t j = 1; j < nCbinsX+1; j++) {
      const Double_t oldCnt = hResponsePtcN -> GetBinContent(j, i);
      const Double_t oldErr = hResponsePtcN -> GetBinError(j, i);
      const Double_t xWidth = hResponsePtcN -> GetXaxis() -> GetBinWidth(j);
      const Double_t yWidth = hResponsePtcN -> GetYaxis() -> GetBinWidth(i);
      const Double_t dArea  = xWidth * yWidth;
      const Double_t newCnt = (oldCnt / cNorm) * dArea;
      const Double_t newErr = (oldErr / cNorm) * dArea;
      hResponsePtcN -> SetBinContent(j, i, newCnt);
      hResponsePtcN -> SetBinError(j, i, newErr);
    }
  }


  // create directory structure
  const Int_t nDir = nJetTypes + 1;
  TDirectory *dir[nDir];
  dir[7] = (TDirectory*) fOut -> mkdir("EventInfo");
  dir[0] = (TDirectory*) fOut -> mkdir("GeantJets");
  dir[1] = (TDirectory*) fOut -> mkdir("MuDstJets");
  dir[2] = (TDirectory*) fOut -> mkdir("Candidates");
  dir[3] = (TDirectory*) fOut -> mkdir("MatchJets");
  dir[4] = (TDirectory*) fOut -> mkdir("JunkJets");
  dir[5] = (TDirectory*) fOut -> mkdir("Mystery");
  dir[6] = (TDirectory*) fOut -> mkdir("NotMatches");

  // write and close output
  fOut          -> cd();
  hEfficiencyA  -> Write();
  hEfficiencyPt -> Write();
  hResponseA    -> Write();
  hResponseAn   -> Write();
  hResponsePt   -> Write();
  hResponsePtN  -> Write();
  hResponsePtc  -> Write();
  hResponsePtcN -> Write();
  dir[7]        -> cd();
  hRefmultG     -> Write();
  hRefmultU     -> Write();
  hNumJetsG     -> Write();
  hNumJetsU     -> Write();
  hNumToMatch   -> Write();
  hNumMatched   -> Write();
  hNumHardG     -> Write();
  hNumHardU     -> Write();
  hGeantArea    -> Write();
  hMatchArea    -> Write();
  hGeantPtCorr  -> Write();
  hMatchPtCorr  -> Write();
  for (Int_t i = 0; i < nJetTypes; i++) {
    dir[i]          -> cd();
    hJetArea[i]     -> Write();
    hJetEta[i]      -> Write();
    hJetPhi[i]      -> Write();
    hJetPt[i]       -> Write();
    hJetPtCorr[i]   -> Write();
    hJetPhiVsEta[i] -> Write();
  }
  for (Int_t i = 0; i < nMatchTypes; i++) {
    Int_t iDir = i + 1;
    dir[iDir]     -> cd();
    hJetQt[i]     -> Write();
    hJetDr[i]     -> Write();
    hJetS[i]      -> Write();
    hJetDp[i]     -> Write();
    hJetQtVsDr[i] -> Write();
    hJetSvsDr[i]  -> Write();
  }
  fOut -> cd();
  fOut -> Close();

  // close input
  fGnt -> cd();
  fGnt -> Close();
  fMu  -> cd();
  fMu  -> Close();


  cout << "  Matching script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
