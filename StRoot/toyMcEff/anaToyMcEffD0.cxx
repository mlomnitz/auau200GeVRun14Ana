/* *********************************************************************
 *
 *  Analysis code to read D0 *.toyMc.root files.
 *
 *  Authors:
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
*/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
#include "iostream"
#include <string>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"

#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH3F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"
#endif

#include "d0Nt.h"
#include "anaCuts.h"

using namespace std;

int getD0PtIndex(float const pt)
{
   int bin = -1;
   for (int i = 0; i < anaCuts::nPtBins; i++)
   {
      if ((pt >= anaCuts::PtEdge[i]) && (pt < anaCuts::PtEdge[i + 1]))
         bin = i;
   }
   return bin;
}

int getD0CentIndex(float const cent)
{
   int bin = -1;

   for (int i = 0; i < anaCuts::physNCentralities; ++i)
   {
      if (cent <= anaCuts::physCentralityEdges[i + 1])
      {
         bin = i;
         break;
      }
   }

   return bin;
}

bool isGoodTrack(float const pt, float const eta)
{
   return pt > anaCuts::pt && fabs(eta) < anaCuts::eta;
}

bool isGoodPair(float const pt, float const y, float const cosTheta, float const pDca, float const kDca,
                float const dca12, float const decayLength, float const dcaV0ToPv)
{
   int tmpIndex = getD0PtIndex(pt);

   return fabs(y) < anaCuts::rapidity &&
          cosTheta > anaCuts::cosTheta[tmpIndex] &&
          pDca > anaCuts::pDca[tmpIndex] && kDca > anaCuts::kDca[tmpIndex] &&
          dca12 < anaCuts::dcaDaughters[tmpIndex] &&
          decayLength > anaCuts::decayLength[tmpIndex] &&
          dcaV0ToPv < anaCuts::dcaV0ToPv[tmpIndex];
}

struct Hists
{
  int centrality;
  TH1D* hNoCuts;
  TH1D* hNoCutsPhysBinning;
  TH1D* hTopoCuts;
  TH1D* hHftMatchingOnly;
  TH1D* hTpcOnly;
  TH1D* hTpcHftTopo;

  Hists(int cent)
  {
    centrality = cent;

    int nBins = 60;
    float minPt = 0.;
    float maxPt = 12.;
    hNoCuts = new TH1D(Form("hNoCuts_%i",cent),Form("No Cuts %s",anaCuts::physCentralityName[cent].Data()),nBins,minPt,maxPt);
    hNoCutsPhysBinning = new TH1D(Form("hNoCutsPhysBinning_%i",cent),Form("No Cuts Physics Binning %s",anaCuts::physCentralityName[cent].Data()),anaCuts::physNPtBins,anaCuts::physPtEdge);
    hTopoCuts = new TH1D(Form("hTopoCuts_%i",cent),Form("Topo Cuts %s",anaCuts::physCentralityName[cent].Data()),nBins,minPt,maxPt);
    hHftMatchingOnly = new TH1D(Form("hHftMatchingOnly_%i",cent),Form("HFT Matching Only %s",anaCuts::physCentralityName[cent].Data()),nBins,minPt,maxPt);
    hTpcOnly = new TH1D(Form("hTpcOnly_%i",cent),Form("TPC Only %s",anaCuts::physCentralityName[cent].Data()),nBins,minPt,maxPt);
    hTpcHftTopo = new TH1D(Form("hTpcHftTopo_%i",cent),Form("TPC + HFT + Topo %s",anaCuts::physCentralityName[cent].Data()),nBins,minPt,maxPt);

    hNoCuts->Sumw2();
    hNoCutsPhysBinning->Sumw2();
    hTopoCuts->Sumw2();
    hHftMatchingOnly->Sumw2();
    hTpcOnly->Sumw2();
    hTpcHftTopo->Sumw2();
  }

  void fill(d0Nt const* const t)
  {
    bool passTopologicalCuts = isGoodPair(t->rPt, t->rY, t->cosTheta, t->pRDca, t->kRDca, t->dca12, t->decayLength, t->dcaD0ToPv);
    bool passHft = t->kHft > 0 && t->pHft > 0;
    bool passTpc = t->kTpc>0 && t->pTpc>0;
    float weight = t->pt * t->w;

    hNoCuts->Fill(t->rPt,weight);
    hNoCutsPhysBinning->Fill(t->rPt,weight);
    if(passTopologicalCuts) hTopoCuts->Fill(t->rPt,weight);
    if(passHft) hHftMatchingOnly->Fill(t->rPt,weight);
    if(passTpc) hTpcOnly->Fill(t->rPt,weight);
    if(passTpc && passHft && passTopologicalCuts) hTpcHftTopo->Fill(t->rPt,weight);
  }

  void makeEffciency(TFile* fOut, TH1D* hPass,TH1D* hTotal,TH1D* hTotalPhysBinning=NULL,bool graphAsym=false)
  {
    fOut->cd();
    TString name = hPass->GetName();
    name.Replace(0,1,"");
    TH1D* hEff = (TH1D*)hPass->Clone(Form("hEff%s",name.Data()));
    hEff->SetTitle(Form("%s Eff.",hPass->GetTitle()));
    hEff->Divide(hTotal);
    hEff->Write();

    if(hTotalPhysBinning)
    {
      TH1* hEffPhysBinning = hPass->Rebin(anaCuts::physNPtBins,Form("hEffPhysBinning%s",name.Data()),anaCuts::physPtEdge);
      hEffPhysBinning->SetTitle(Form("%s Eff. - Physics Binning",hPass->GetTitle()));
      hEffPhysBinning->Divide(hTotalPhysBinning);
      hEffPhysBinning->Write();
    }

    if(graphAsym)
    {
      TGraphAsymmErrors* grEff = new TGraphAsymmErrors(hPass,hTotal);
      grEff->SetName(Form("grEff%s",name.Data()));
      grEff->SetTitle(Form("%s Eff.",hPass->GetTitle()));
      grEff->Write();
    }
  }

  void write(TFile* fOut)
  {
    hNoCuts->Write();
    hTopoCuts->Write();
    hHftMatchingOnly->Write();
    hTpcOnly->Write();
    hTpcHftTopo->Write();
    makeEffciency(fOut,hTopoCuts,hNoCuts,hNoCutsPhysBinning);
    makeEffciency(fOut,hHftMatchingOnly,hNoCuts,hNoCutsPhysBinning);
    makeEffciency(fOut,hTpcOnly,hNoCuts,hNoCutsPhysBinning);
    makeEffciency(fOut,hTpcHftTopo,hNoCuts,hNoCutsPhysBinning);
  }
};

struct TopoHists
{
  TH3F* mcPointingAngle;
  TH3F* mcDecayL;
  TH3F* mcDca12;
  TH3F* mcPionDca2Vtx;
  TH3F* mcKaonDca2Vtx;
  TH3F* mcD0Dca2Vtx;
  
  TopoHists()
  {
    mcPointingAngle = new TH3F(Form("%s_se_us_pointingangle", "mc"), "Same Event US pointing angle; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 1000, 0.9, 1.0);
    mcDecayL = new TH3F(Form("%s_se_us_decayL", "mc"), "Same Event US Decay Length; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.1);
    mcDca12 = new TH3F(Form("%s_se_us_dcaDaughters", "mc"), "Same Event US dca daughters; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.05);
    mcPionDca2Vtx = new TH3F(Form("%s_se_us_pionDca", "mc"), "Same Event #pi dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.2);
    mcKaonDca2Vtx = new TH3F(Form("%s_se_us_kaonDca", "mc"), "Same Event US K dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.2);
    mcD0Dca2Vtx = new TH3F(Form("%s_se_us_D0Dca2Vtx", "mc"), "SameEvent US D0 dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.05);

    mcPointingAngle->Sumw2();
    mcDecayL->Sumw2();
    mcDca12->Sumw2();
    mcPionDca2Vtx->Sumw2();
    mcKaonDca2Vtx->Sumw2();
    mcD0Dca2Vtx->Sumw2();
  }

  void fill(d0Nt* t)
  {
    if (t->rM <  anaCuts::massMin || t->rM > anaCuts::massMax) return;
    bool passHft = t->kHft > 0 && t->pHft > 0;
    if (!passHft) return;

    float weight = t->pt * t->w;
    int ptIndex = getD0PtIndex(t->rPt);

    //Cos theta
    if (t->pRDca > anaCuts::pDca[ptIndex] && t->kRDca > anaCuts::kDca[ptIndex] &&
        t->dca12 < anaCuts::dcaDaughters[ptIndex] && t->decayLength > anaCuts::decayLength[ptIndex] &&
        //std::cos(t->pointingAngle()) > anaCuts::cosTheta[ptIndex] &&
        t->dcaD0ToPv < anaCuts::dcaV0ToPv[ptIndex]) mcPointingAngle->Fill(t->rPt, t->cent,  t->cosTheta, weight);

    //DecayL
    if (t->pRDca > anaCuts::pDca[ptIndex] && t->kRDca > anaCuts::kDca[ptIndex] &&
        t->dca12 < anaCuts::dcaDaughters[ptIndex] &&
        //t->decayLength() > anaCuts::decayLength[ptIndex] &&
        t->cosTheta > anaCuts::cosTheta[ptIndex] && t->dcaD0ToPv < anaCuts::dcaV0ToPv[ptIndex]) mcDecayL->Fill(t->rPt, t->cent, t->decayLength/1.e4, weight);

    //DcaDaughter
    if (t->pRDca > anaCuts::pDca[ptIndex] && t->kRDca > anaCuts::kDca[ptIndex] &&
        //t->dca12 < anaCuts::dcaDaughters[ptIndex] &&
        t->decayLength > anaCuts::decayLength[ptIndex] &&
        t->cosTheta > anaCuts::cosTheta[ptIndex] &&
        t->dcaD0ToPv < anaCuts::dcaV0ToPv[ptIndex])
      mcDca12->Fill(t->rPt, t->cent, t->dca12/1.e4, weight);

    //PionDca
    if (//t->pRDca > anaCuts::pDca[ptIndex] &&
        t->kRDca > anaCuts::kDca[ptIndex] &&
        t->dca12 < anaCuts::dcaDaughters[ptIndex] && t->decayLength > anaCuts::decayLength[ptIndex] &&
        t->cosTheta > anaCuts::cosTheta[ptIndex] && t->dcaD0ToPv < anaCuts::dcaV0ToPv[ptIndex]) mcPionDca2Vtx->Fill(t->rPt, t->cent, t->pRDca/1.e4, weight);

    //Kaon Dca
    if (t->pRDca > anaCuts::pDca[ptIndex] &&
        //t->kRDca > anaCuts::kDca[ptIndex] &&
        t->dca12 < anaCuts::dcaDaughters[ptIndex] &&
        t->decayLength > anaCuts::decayLength[ptIndex] && t->cosTheta > anaCuts::cosTheta[ptIndex] &&
        t->dcaD0ToPv < anaCuts::dcaV0ToPv[ptIndex]) mcKaonDca2Vtx->Fill(t->rPt, t->cent,  t->kRDca/1.e4, weight);

    //D0 dca
    if (t->pRDca > anaCuts::pDca[ptIndex] && t->kRDca > anaCuts::kDca[ptIndex] &&
        t->dca12 < anaCuts::dcaDaughters[ptIndex] && t->decayLength > anaCuts::decayLength[ptIndex] &&
        t->cosTheta > anaCuts::cosTheta[ptIndex]) mcD0Dca2Vtx->Fill(t->rPt, t->cent, t->dcaD0ToPv/1.e4, weight);

  }

  void write(TFile* fOut)
  {
    fOut->cd();
    mcPointingAngle->Write();
    mcDecayL->Write();
    mcDca12->Write();
    mcPionDca2Vtx->Write();
    mcKaonDca2Vtx->Write();
    mcD0Dca2Vtx->Write();
  }
};

int main(int argc, char **argv)
{
   d0Nt* t = new d0Nt();

   TFile* fOut = new TFile("eff.root", "recreate");

   Long64_t nEntries = t->GetEntries();
   cout << "nEntries = " << nEntries << endl;

   std::vector<Hists> hists;
   TopoHists          topoHists;

   for(int iCent = 0; iCent < anaCuts::physNCentralities; ++iCent)
   {
     hists.push_back(Hists(iCent));
   }

   for (Long64_t i = 0; i < t->GetEntries(); ++i)
   {
      t->GetEntry(i);

      if (i && i % 1000000 == 0) cout << static_cast<float>(i) / nEntries << endl;

      if (fabs(t->y) > anaCuts::rapidity) continue;
      if (!isGoodTrack(t->kRPt, t->kREta) || !isGoodTrack(t->pRPt, t->pREta)) continue;

      int d0CentBin = getD0CentIndex(t->cent);
      if(d0CentBin<0) continue;
      hists[d0CentBin].fill(t);

      topoHists.fill(t);
   }

   for(int iCent = 0; iCent < anaCuts::physNCentralities; ++iCent)
   {
     hists[iCent].write(fOut);
   }

   topoHists.write(fOut);
   fOut->Close();
}
