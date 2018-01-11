// 'MakePlots.C'
// Derek Anderson
// 12.06.2016
//
// Use this to make plots from matching output.

#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;

// jet parameters
static const Int_t    Nrm    = 1;
static const Int_t    Type   = 1;
static const Bool_t   Ghosts = false;
static const Double_t Rjet   = 0.3;
// plot parameters
static const Double_t A1 = 0.;
static const Double_t A2 = 1.;
static const Double_t P1 = -5.;
static const Double_t P2 = 60.;
static const Double_t Q1 = 0.;
static const Double_t Q2 = 5.;
static const Double_t R1 = 0.;
static const Double_t R2 = 5.;
// style parameters
static const Int_t    C0     = 1.;
static const Int_t    C1     = 810;
static const Int_t    C2     = 890;
static const Int_t    M0     = 1;
static const Int_t    M1     = 4;
static const Int_t    Font   = 42;
static const Int_t    Center = 1;
static const Double_t Label  = 0.02;


void MakePlots() {

  cout << "\n  Making plots..." << endl;


  TFile *fOut = new TFile("pp200r12.plots.r03a02rm1full.d19m7y2017.root", "recreate");
  TFile *fIn  = new TFile("pp200r12.match.r03a02rm1full.d17m7y2017.root");
  if (!fIn) {
    cerr << "PANIC: input file could not be opened!" << endl;
    assert(fIn);
  }

  cout << "    Grabbing histograms..." << endl;

  const Int_t nJ = 3;
  const Int_t nM = 2;
  // evt histograms
  TH1D *hEa;
  TH1D *hEpt;
  TH2D *hRa;
  TH2D *hRpt;
  // jet histograms [0='G', 1='U', 2='M']
  TH1D *hA[nJ];
  TH1D *hPtc[nJ];
  // matching calculation [0='U', 1='M']
  TH1D *hQt[nM];
  TH1D *hDr[nM];
  TH2D *hQtVsDr[nM];
  if (fIn) {
    fIn -> GetObject("hEfficiencyA", hEa);
    fIn -> GetObject("hEfficiencyPt", hEpt);
    fIn -> GetObject("hResponseA", hRa);
    fIn -> GetObject("hResponsePtc", hRpt);
    fIn -> GetObject("GeantJets/hJetAreaG", hA[0]);
    fIn -> GetObject("MuDstJets/hJetAreaU", hA[1]);
    fIn -> GetObject("MatchJets/hJetAreaM", hA[2]);
    fIn -> GetObject("GeantJets/hJetPtCorrG", hPtc[0]);
    fIn -> GetObject("MuDstJets/hJetPtCorrU", hPtc[1]);
    fIn -> GetObject("MatchJets/hJetPtCorrM", hPtc[2]);
    fIn -> GetObject("MuDstJets/hJetQtU", hQt[0]);
    fIn -> GetObject("MatchJets/hJetQtM", hQt[1]);
    fIn -> GetObject("MuDstJets/hJetDrU", hDr[0]);
    fIn -> GetObject("MatchJets/hJetDrM", hDr[1]);
    fIn -> GetObject("MuDstJets/hJetQtVsDrU", hQtVsDr[0]);
    fIn -> GetObject("MatchJets/hJetQtVsDrM", hQtVsDr[1]);
  }
  fOut -> cd();


  cout << "    Setting styles..." << endl;

  // set styles
  hEa     -> SetTitleFont(Font);
  hEa     -> GetXaxis() -> SetRangeUser(A1, A2);
  hEa     -> GetXaxis() -> SetTitleFont(Font);
  hEa     -> GetXaxis() -> CenterTitle(Center);
  hEa     -> GetXaxis() -> SetLabelSize(Label);
  hEa     -> GetYaxis() -> SetTitleFont(Font);
  hEa     -> GetYaxis() -> CenterTitle(Center);
  hEa     -> GetYaxis() -> SetLabelSize(Label);
  hEpt    -> SetTitleFont(Font);
  hEpt    -> GetXaxis() -> SetRangeUser(P1, P2);
  hEpt    -> GetXaxis() -> SetTitleFont(Font);
  hEpt    -> GetXaxis() -> CenterTitle(Center);
  hEpt    -> GetXaxis() -> SetLabelSize(Label);
  hEpt    -> GetYaxis() -> SetTitleFont(Font);
  hEpt    -> GetYaxis() -> CenterTitle(Center);
  hEpt    -> GetYaxis() -> SetLabelSize(Label);
  hRa     -> SetTitleFont(Font);
  hRa     -> GetXaxis() -> SetRangeUser(A1, A2);
  hRa     -> GetXaxis() -> SetTitleFont(Font);
  hRa     -> GetXaxis() -> CenterTitle(Center);
  hRa     -> GetXaxis() -> SetLabelSize(Label);
  hRa     -> GetYaxis() -> SetRangeUser(A1, A2);
  hRa     -> GetYaxis() -> SetTitleFont(Font);
  hRa     -> GetYaxis() -> CenterTitle(Center);
  hRa     -> GetYaxis() -> SetLabelSize(Label);
  hRpt    -> SetTitleFont(Font);
  hRpt    -> GetXaxis() -> SetRangeUser(P1, P2);
  hRpt    -> GetXaxis() -> SetTitleFont(Font);
  hRpt    -> GetXaxis() -> CenterTitle(Center);
  hRpt    -> GetXaxis() -> SetLabelSize(Label);
  hRpt    -> GetYaxis() -> SetRangeUser(P1, P2);
  hRpt    -> GetYaxis() -> SetTitleFont(Font);
  hRpt    -> GetYaxis() -> CenterTitle(Center);
  hRpt    -> GetYaxis() -> SetLabelSize(Label);
  hA[0]   -> SetLineColor(C0);
  hA[0]   -> SetMarkerColor(C0);
  hA[0]   -> SetTitleFont(Font);
  hA[0]   -> GetXaxis() -> SetRangeUser(A1, A2);
  hA[0]   -> GetXaxis() -> SetTitleFont(Font);
  hA[0]   -> GetXaxis() -> CenterTitle(Center);
  hA[0]   -> GetXaxis() -> SetLabelSize(Label);
  hA[0]   -> GetYaxis() -> SetTitleFont(Font);
  hA[0]   -> GetYaxis() -> CenterTitle(Center);
  hA[0]   -> GetYaxis() -> SetLabelSize(Label);
  hPtc[0] -> SetLineColor(C0);
  hPtc[0] -> SetMarkerColor(C0);
  hPtc[0] -> SetTitleFont(Font);
  hPtc[0] -> GetXaxis() -> SetRangeUser(P1, P2);
  hPtc[0] -> GetXaxis() -> SetTitleFont(Font);
  hPtc[0] -> GetXaxis() -> CenterTitle(Center);
  hPtc[0] -> GetXaxis() -> SetLabelSize(Label);
  hPtc[0] -> GetYaxis() -> SetTitleFont(Font);
  hPtc[0] -> GetYaxis() -> CenterTitle(Center);
  hPtc[0] -> GetYaxis() -> SetLabelSize(Label);
  hA[1]   -> SetLineColor(C1);
  hA[1]   -> SetMarkerColor(C1);
  hA[1]   -> SetTitleFont(Font);
  hA[1]   -> GetXaxis() -> SetRangeUser(A1, A2);
  hA[1]   -> GetXaxis() -> SetTitleFont(Font);
  hA[1]   -> GetXaxis() -> CenterTitle(Center);
  hA[1]   -> GetXaxis() -> SetLabelSize(Label);
  hA[1]   -> GetYaxis() -> SetTitleFont(Font);
  hA[1]   -> GetYaxis() -> CenterTitle(Center);
  hA[1]   -> GetYaxis() -> SetLabelSize(Label);
  hPtc[1] -> SetLineColor(C1);
  hPtc[1] -> SetMarkerColor(C1);
  hPtc[1] -> SetTitleFont(Font);
  hPtc[1] -> GetXaxis() -> SetRangeUser(P1, P2);
  hPtc[1] -> GetXaxis() -> SetTitleFont(Font);
  hPtc[1] -> GetXaxis() -> CenterTitle(Center);
  hPtc[1] -> GetXaxis() -> SetLabelSize(Label);
  hPtc[1] -> GetYaxis() -> SetTitleFont(Font);
  hPtc[1] -> GetYaxis() -> CenterTitle(Center);
  hPtc[1] -> GetYaxis() -> SetLabelSize(Label);
  hA[2]   -> SetLineColor(C2);
  hA[2]   -> SetMarkerColor(C2);
  hA[2]   -> SetTitleFont(Font);
  hA[2]   -> GetXaxis() -> SetRangeUser(A1, A2);
  hA[2]   -> GetXaxis() -> SetTitleFont(Font);
  hA[2]   -> GetXaxis() -> CenterTitle(Center);
  hA[2]   -> GetXaxis() -> SetLabelSize(Label);
  hA[2]   -> GetYaxis() -> SetTitleFont(Font);
  hA[2]   -> GetYaxis() -> CenterTitle(Center);
  hA[2]   -> GetYaxis() -> SetLabelSize(Label);
  hPtc[2] -> SetLineColor(C2);
  hPtc[2] -> SetMarkerColor(C2);
  hPtc[2] -> SetTitleFont(Font);
  hPtc[2] -> GetXaxis() -> SetRangeUser(P1, P2);
  hPtc[2] -> GetXaxis() -> SetTitleFont(Font);
  hPtc[2] -> GetXaxis() -> CenterTitle(Center);
  hPtc[2] -> GetXaxis() -> SetLabelSize(Label);
  hPtc[2] -> GetYaxis() -> SetTitleFont(Font);
  hPtc[2] -> GetYaxis() -> CenterTitle(Center);
  hPtc[2] -> GetYaxis() -> SetLabelSize(Label);
  for (Int_t i = 0; i < nM; i++) {
    if (i == 0) {
      hQt[i] -> SetLineColor(C1);
      hQt[i] -> SetMarkerColor(C1);
      hQt[i] -> SetMarkerStyle(M0);
      hDr[i] -> SetLineColor(C1);
      hDr[i] -> SetMarkerColor(C1);
      hDr[i] -> SetMarkerStyle(M0);
    }
    if (i == 1) {
      hQt[i] -> SetLineColor(C2);
      hQt[i] -> SetMarkerColor(C2);
      hQt[i] -> SetMarkerStyle(M1);
      hDr[i] -> SetLineColor(C2);
      hDr[i] -> SetMarkerColor(C2);
      hDr[i] -> SetMarkerStyle(M1);
    }
    hQt[i]     -> SetTitleFont(Font);
    hQt[i]     -> GetXaxis() -> SetRangeUser(Q1, Q2);
    hQt[i]     -> GetXaxis() -> SetTitleFont(Font);
    hQt[i]     -> GetXaxis() -> CenterTitle(Center);
    hQt[i]     -> GetXaxis() -> SetLabelSize(Label);
    hQt[i]     -> GetYaxis() -> SetTitleFont(Font);
    hQt[i]     -> GetYaxis() -> CenterTitle(Center);
    hQt[i]     -> GetYaxis() -> SetLabelSize(Label);
    hDr[i]     -> SetTitleFont(Font);
    hDr[i]     -> GetXaxis() -> SetRangeUser(R1, R2);
    hDr[i]     -> GetXaxis() -> SetTitleFont(Font);
    hDr[i]     -> GetXaxis() -> CenterTitle(Center);
    hDr[i]     -> GetXaxis() -> SetLabelSize(Label);
    hDr[i]     -> GetYaxis() -> SetTitleFont(Font);
    hDr[i]     -> GetYaxis() -> CenterTitle(Center);
    hDr[i]     -> GetYaxis() -> SetLabelSize(Label);
    hQtVsDr[i] -> SetTitleFont(Font);
    hQtVsDr[i] -> GetXaxis() -> SetRangeUser(R1, R2);
    hQtVsDr[i] -> GetXaxis() -> SetTitleFont(Font);
    hQtVsDr[i] -> GetXaxis() -> CenterTitle(Center);
    hQtVsDr[i] -> GetXaxis() -> SetLabelSize(Label);
    hQtVsDr[i] -> GetYaxis() -> SetRangeUser(Q1, Q2);
    hQtVsDr[i] -> GetYaxis() -> SetTitleFont(Font);
    hQtVsDr[i] -> GetYaxis() -> CenterTitle(Center);
    hQtVsDr[i] -> GetYaxis() -> SetLabelSize(Label);
  }


  cout << "    Creating legends and labels..." << endl;

  const Int_t CL = 0;
  const Int_t FL = 0;
  TLegend *lGB = new TLegend(0.1, 0.1, 0.3, 0.3);
  lGB -> SetLineColor(CL);
  lGB -> SetLineStyle(FL);
  lGB -> SetFillColor(CL);
  lGB -> SetFillStyle(FL);
  lGB -> SetTextFont(Font);
  lGB -> AddEntry(hA[0], "particle jets");
  lGB -> AddEntry(hA[1], "all detector jets");
  lGB -> AddEntry(hA[2], "matched detector jets");

  TLegend *lMB = new TLegend(0.1, 0.1, 0.3, 0.3);
  lMB -> SetLineColor(CL);
  lMB -> SetLineStyle(FL);
  lMB -> SetFillColor(CL);
  lMB -> SetFillStyle(FL);
  lMB -> SetTextFont(Font);
  lMB -> AddEntry(hQt[0], "all detector jets");
  lMB -> AddEntry(hQt[1], "matched detector jets");


  // create label
  TString line1;
  TString line2;
  TString line3;
  TString line4;
  if (Rjet == 0.3) {
    line1 = "Anti-k_{T}, R = 0.3";
    line2 = "A_{jet} > 0.2, N_{rm} = ";
    line2 += Nrm;
    if (Ghosts)
      line3 = "Active area, explicit ghosts";
    else
      line3 = "Active area, implicit ghosts";
  }
  else if (Rjet == 0.5) {
    line1 = "Anti-k_{T}, R = 0.5";
    line2 = "A_{jet} > 0.65, N_{rm} = ";
    line2 += Nrm;
    if (Ghosts)
      line3 = "Active area, explicit ghosts";
    else
      line3 = "Active area, implicit ghosts";
  }
  else if (Rjet == 0.7) {
    line1 = "Anti-k_{T}, R = 0.7";
    line2 = "A_{jet} > 1.2, N_{rm} = ";
    line2 += Nrm;
    if (Ghosts)
      line3 = "Active area, explicit ghosts";
    else
      line3 = "Active area, implicit ghosts";
  }
  else {
    cerr << "PANIC: Whoah! Check what Rjet is set to..." << endl;
    assert(0);
  }
  switch (Type) {
    case 0:
      line4 = "Charged jets";
      break;
    case 1:
      line4 = "Full jets";
      break;
    case 2:
      line4 = "Neutral jets";
      break;
    default:
      cerr << "PANIC: Whoah! Check what Type is set to..." << endl;
      assert(0);
      break;
  }

  TPaveText *pt = new TPaveText(0.3, 0.3, 0.5, 0.5, "NDC NB");
  pt -> SetFillColor(CL);
  pt -> SetFillStyle(FL);
  pt -> SetLineColor(CL);
  pt -> SetLineStyle(FL);
  pt -> SetTextFont(Font);
  pt -> AddText(line1.Data());
  pt -> AddText(line2.Data());
  pt -> AddText(line3.Data());
  pt -> AddText(line4.Data());


  // create evt titles
  TString sEa("Reconstruction efficiency, #epsilon(A_{jet}) = N_{matched}(A_{jet}) / N_{pythia}(A_{jet}); A_{jet}; #epsilon(A_{jet})");
  TString sEpt("Reconstruction efficiency, #epsilon(p_{T}^{jet}) = N_{matched}(p_{T}^{jet}) / N_{pythia}(p_{T}^{jet}); p_{T}^{jet}; #epsilon(p_{T}^{jet})");
  TString sRa("Response, jet area; A_{jet}(matched); A_{jet}(particle)");
  TString sRpt("Response, jet p_{T}; p_{T}^{jet}(matched); p_{T}^{jet}(particle)");
  // create jet titles
  TString sAgb("Jet Area; A_{jet}(particle), #color[810]{A_{jet}(detector)}, #color[890]{A_{jet}(matched)}; (1/N_{evt})dN_{jet}/dA_{jet}");
  TString sPgb("Jet p_{T}^{corr}; p_{T}^{corr}(particle), #color[810]{p_{T}^{corr}(detector)}, #color[890]{p_{T}^{corr}(matched)}; (1/N_{evt})dN_{jet}/dp_{T}^{corr}");
  // create matching titles
  TString sQtb("Jet q_{T}; #color[810]{q_{T}(detector)}, #color[890]{q_{T}(matched)}; counts");
  TString sDrb("Jet #Deltar; #color[810]{#Deltar(detector)}, #color[890]{#Deltar(matched)}; counts");
  TString sQtVsDru("Jet q_{T} vs. #Deltar, all detector jets; #Deltar; q_{T}");
  TString sQtVsDrm("Jet q_{T} vs. #Deltar, matched detector jets; #Deltar; q_{T}");


  cout << "    Drawing plots..." << endl;

  // draw evt plots
  const Int_t width  = 800;
  const Int_t height = 800;
  TCanvas *cEa = new TCanvas("cEa", "Eff(Ajet)", width, height);
  cEa -> cd();
  cEa -> SetGrid(0, 0);
  hEa -> SetTitle(sEa);
  hEa -> Draw();
  pt  -> Draw();
  cEa -> Write();
  cEa -> Close();

  TCanvas *cEpt = new TCanvas("cEpt", "Eff(pTcorr)", width, height);
  cEpt -> cd();
  cEpt -> SetGrid(0, 0);
  hEpt -> SetTitle(sEpt);
  hEpt -> Draw();
  pt   -> Draw();
  cEpt -> Write();
  cEpt -> Close();

  TCanvas *cRa = new TCanvas("cRa", "Response(Ajet)", width, height);
  cRa -> cd();
  cRa -> SetGrid(0, 0);
  cRa -> SetLogz(1);
  hRa -> SetTitle(sRa);
  hRa -> Draw("colz");
  pt  -> Draw();
  cRa -> Write();
  cRa -> Close();

  TCanvas *cRpt = new TCanvas("cRpt", "Response(pTcorr)", width, height);
  cRpt -> cd();
  cRpt -> SetGrid(0, 0);
  cRpt -> SetLogz(1);
  hRpt -> SetTitle(sRpt);
  hRpt -> Draw("colz");
  pt   -> Draw();
  cRpt -> Write();
  cRpt -> Close();

  // draw jet plots
  TCanvas *cAgb = new TCanvas("cAgb", "Ajet", width, height);
  cAgb  -> cd();
  cAgb  -> SetGrid(0, 0);
  cAgb  -> SetLogy(1);
  hA[0] -> SetTitle(sAgb);
  hA[0] -> Draw();
  hA[1] -> Draw("same");
  hA[2] -> Draw("same");
  lGB   -> Draw();
  pt    -> Draw();
  cAgb  -> Write();
  cAgb  -> Close();

  TCanvas *cPgb = new TCanvas("cPgb", "pTcorr", width, height);
  cPgb    -> cd();
  cPgb    -> SetGrid(0, 0);
  cPgb    -> SetLogy(1);
  hPtc[0] -> SetTitle(sPgb);
  hPtc[0] -> Draw();
  hPtc[1] -> Draw("same");
  hPtc[2] -> Draw("same");
  lGB     -> Draw();
  pt      -> Draw();
  cPgb    -> Write();
  cPgb    -> Close();

  // draw matching plots
  TCanvas *cQtu = new TCanvas("cQtb", "qT", width, height);
  cQtb   -> cd();
  cQtb   -> SetGrid(0, 0);
  hQt[0] -> SetTitle(sQtb);
  hQt[0] -> Draw();
  hQt[1] -> Draw("same");
  lMB    -> Draw();
  pt     -> Draw();
  cQtb   -> Write();
  cQtb   -> Close();

  TCanvas *cDru = new TCanvas("cDrb", "dR", width, height);
  cDrb   -> cd();
  cDrb   -> SetGrid(0, 0);
  cDrb   -> SetLogy(1);
  hDr[0] -> SetTitle(sDrb);
  hDr[0] -> Draw();
  hDr[1] -> Draw("same");
  lMB    -> Draw();
  pt     -> Draw();
  cDrb   -> Write();
  cDrb   -> Close();

  TCanvas *cQtVsDru = new TCanvas("cQtVsDru", "qT vs. dR, all MuDst jets", width, height);
  cQtVsDru   -> cd();
  cQtVsDru   -> SetGrid(0, 0);
  cQtVsDru   -> SetLogz(1);
  hQtVsDr[0] -> SetTitle(sQtVsDru);
  hQtVsDr[0] -> Draw("colz");
  pt         -> Draw();
  cQtVsDru   -> Write();
  cQtVsDru   -> Close();

  TCanvas *cQtVsDrm = new TCanvas("cQtVsDrm", "qT vs. dR, matched MuDst jets", width, height);
  cQtVsDrm   -> cd();
  cQtVsDrm   -> SetGrid(0, 0);
  cQtVsDrm   -> SetLogz(1);
  hQtVsDr[1] -> SetTitle(sQtVsDrm);
  hQtVsDr[1] -> Draw("colz");
  pt         -> Draw();
  cQtVsDrm   -> Write();
  cQtVsDrm   -> Close();


  cout << "    Closing files..." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();


  cout << "  Plots made!\n" << endl;

}

// End ------------------------------------------------------------------------
