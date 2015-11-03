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
#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"

#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#endif

#include "dataDrivenFastSimulator.h"
#include "d0Nt.h"
#include "anaCuts.h"

using namespace std;

TGraphErrors* gHftRatioCorrection = NULL;
TF1*          gf1AuAu010Weight = NULL;

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

bool isGoodTrack(float const ptCut, float const pt, float const eta)
{
   return pt > ptCut && fabs(eta) < anaCuts::eta;
}

bool isGoodPair(float const pt, float const y, float const cosTheta, float const pDca, float const kDca,
                float const dca12, float const decayLength, float const dcaV0ToPv)
{
   int tmpIndex = getD0PtIndex(pt);
    if(tmpIndex < 0) return false;

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
  float minPtCut;
  TH1D* hNoCuts;
  TH1D* hNoCutsPhysBinning;
  TH1D* hTopoCuts;
  TH1D* hHftMatchingOnly;
  TH1D* hTpcOnly;
  TH1D* hTpcHftTopo;
  TH2D* h2MassPt;

  Hists(float ptCut,int cent)
  {
    centrality = cent;
    minPtCut = ptCut;

    int nBins = 300;
    float minPt = 0.;
    float maxPt = 15.;
    hNoCuts = new TH1D(Form("hNoCuts_minPt%i_%i",(int)(minPtCut*1e3),cent),Form("No Cuts %s minPt>%1.2f",anaCuts::physCentralityName[cent].Data(),minPtCut),nBins,minPt,maxPt);
    hNoCutsPhysBinning = new TH1D(Form("hNoCutsPhysBinning_minPt%i_%i",(int)(minPtCut*1e3),cent),Form("No Cuts Physics Binning %s minPt>%1.2f",anaCuts::physCentralityName[cent].Data(),minPtCut),anaCuts::physNPtBins,anaCuts::physPtEdge);
    hTopoCuts = new TH1D(Form("hTopoCuts_minPt%i_%i",(int)(minPtCut*1e3),cent),Form("Topo Cuts %s minPt>%1.2f",anaCuts::physCentralityName[cent].Data(),minPtCut),nBins,minPt,maxPt);
    hHftMatchingOnly = new TH1D(Form("hHftMatchingOnly_minPt%i_%i",(int)(minPtCut*1e3),cent),Form("HFT Matching Only %s minPt>%1.2f",anaCuts::physCentralityName[cent].Data(),minPtCut),nBins,minPt,maxPt);
    hTpcOnly = new TH1D(Form("hTpcOnly_minPt%i_%i",(int)(minPtCut*1e3),cent),Form("TPC Only %s minPt>%1.2f",anaCuts::physCentralityName[cent].Data(),minPtCut),nBins,minPt,maxPt);
    hTpcHftTopo = new TH1D(Form("hTpcHftTopo_minPt%i_%i",(int)(minPtCut*1e3),cent),Form("TPC + HFT + Topo %s minPt>%1.2f",anaCuts::physCentralityName[cent].Data(),minPtCut),nBins,minPt,maxPt);
    h2MassPt = new TH2D(Form("h2MassPt_minPt%i_%i",(int)(minPtCut*1e3),cent),Form("Invariant Mass vs. Pt %s minPt>%1.2f",anaCuts::physCentralityName[cent].Data(),minPtCut),nBins,minPt,maxPt,210,0,2.1);

    hNoCuts->Sumw2();
    hNoCutsPhysBinning->Sumw2();
    hTopoCuts->Sumw2();
    hHftMatchingOnly->Sumw2();
    hTpcOnly->Sumw2();
    hTpcHftTopo->Sumw2();
    h2MassPt->Sumw2();
  }

  void fill(d0Nt const* const t)
  {
    if(t->rPt > 10.) return;
    bool passTopologicalCuts = isGoodPair(t->rPt, t->rY, t->cosTheta, t->pRDca, t->kRDca, t->dca12, t->decayLength, t->dcaD0ToPv);
    // bool passHft = t->kHft > 0 && t->pHft > 0;
    bool passTpc = t->kTpc>0 && t->pTpc>0;
    // float weight = t->pt * t->w;
    float weight = gf1AuAu010Weight->Eval(t->pt);

    float hftRatioWeight = matchHft(0,t->vz,t->cent,t->pRPt,t->pRPhi,t->pREta);
    hftRatioWeight      *= matchHft(1,t->vz,t->cent,t->kRPt,t->kRPhi,t->kREta);

    float hftRatioCorrection = 1.0;

    hNoCuts->Fill(t->rPt,weight);
    hNoCutsPhysBinning->Fill(t->rPt,weight);
    if (!isGoodTrack(minPtCut,t->kRPt, t->kREta) || !isGoodTrack(minPtCut,t->pRPt, t->pREta)) return;

    /*if(t->pRPt>0.6 && t->pRPt<1.99)
    {
      hftRatioCorrection *= gHftRatioCorrection->Eval(t->pRPt);
    }

    if(t->kRPt>0.6 && t->kRPt<1.99)
    {
      hftRatioCorrection *= gHftRatioCorrection->Eval(t->kRPt);
    }
    */

    if(passTopologicalCuts) hTopoCuts->Fill(t->rPt,weight);
    if(passTpc) hTpcOnly->Fill(t->rPt,weight);

    hHftMatchingOnly->Fill(t->rPt,weight * hftRatioCorrection * hftRatioWeight);

    if(passTpc && passTopologicalCuts) 
    {
      hTpcHftTopo->Fill(t->rPt,weight * hftRatioCorrection * hftRatioWeight);
      h2MassPt->Fill(t->rPt,t->rM,weight);
    }
  }

  void makeEffciency(TDirectory* fOut, TH1D* hPass,TH1D* hTotal,TH1D* hTotalPhysBinning=NULL,bool graphAsym=false)
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
      TGraphAsymmErrors* grEff = new TGraphAsymmErrors(hPass,hTotal,"n");
      grEff->SetName(Form("grEff%s",name.Data()));
      grEff->SetTitle(Form("%s Eff.",hPass->GetTitle()));
      grEff->Write();

      if(hTotalPhysBinning)
      {
        TH1* hPassPhysBinning = hPass->Rebin(anaCuts::physNPtBins,Form("hPassPhysBinning%s",name.Data()),anaCuts::physPtEdge);
        TGraphAsymmErrors* grEffPhysBinning = new TGraphAsymmErrors(hPassPhysBinning,hTotalPhysBinning,"n");
        grEffPhysBinning->SetName(Form("grEffPhysBinning%s",name.Data()));
        grEffPhysBinning->SetTitle(Form("%s Eff. Physics Binning",hPass->GetTitle()));
        grEffPhysBinning->Write();
        delete hPassPhysBinning;
      }
    }
  }

  void write(TFile* fOut)
  {
    TDirectory* dir = NULL;
    if(!(dir = (TDirectory*)fOut->Get(Form("minPt%iMeV",(int)(minPtCut*1e3)))))
        dir = (TDirectory*)fOut->mkdir(Form("minPt%iMeV",(int)(minPtCut*1e3)));

    dir->cd();
    hNoCuts->Write();
    hTopoCuts->Write();
    hHftMatchingOnly->Write();
    hTpcOnly->Write();
    hTpcHftTopo->Write();
    h2MassPt->Write();
    makeEffciency(dir,hTopoCuts,hNoCuts,hNoCutsPhysBinning,true);
    makeEffciency(dir,hHftMatchingOnly,hNoCuts,hNoCutsPhysBinning,true);
    makeEffciency(dir,hTpcOnly,hNoCuts,hNoCutsPhysBinning,true);
    makeEffciency(dir,hTpcHftTopo,hNoCuts,hNoCutsPhysBinning,true);
  }
};

struct TopoHists
{
  float minPtCut;
  TH3F* mcPointingAngle;
  TH3F* mcDecayL;
  TH3F* mcDca12;
  TH3F* mcPionDca2Vtx;
  TH3F* mcKaonDca2Vtx;
  TH3F* mcD0Dca2Vtx;
  
  TopoHists(float const minPt)
  {
    float const maxDca = 0.3;

    int const nDcaBins = 300;
    float const maxDca12 = 0.1;
    int const nDca12Bins = 200;

    minPtCut = minPt;
    mcPointingAngle = new TH3F(Form("%s_se_us_pointingangle_minPt%i", "mc",(int)(minPtCut*1.e3)), "Same Event US pointing angle; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 1000, 0.9, 1.0);
    mcDecayL = new TH3F(Form("%s_se_us_decayL_minPt%i", "mc",(int)(minPtCut*1.e3)), "Same Event US Decay Length; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins, 0, maxDca);
    mcDca12 = new TH3F(Form("%s_se_us_dcaDaughters_minPt%i", "mc",(int)(minPtCut*1.e3)), "Same Event US dca daughters; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDca12Bins, 0, maxDca12);
    mcPionDca2Vtx = new TH3F(Form("%s_se_us_pionDca_minPt%i", "mc",(int)(minPtCut*1.e3)), "Same Event #pi dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins, 0, maxDca);
    mcKaonDca2Vtx = new TH3F(Form("%s_se_us_kaonDca_minPt%i", "mc",(int)(minPtCut*1.e3)), "Same Event US K dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins, 0, maxDca);
    mcD0Dca2Vtx = new TH3F(Form("%s_se_us_D0Dca2Vtx_minPt%i", "mc",(int)(minPtCut*1.e3)), "SameEvent US D0 dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.05);

    mcPointingAngle->Sumw2();
    mcDecayL->Sumw2();
    mcDca12->Sumw2();
    mcPionDca2Vtx->Sumw2();
    mcKaonDca2Vtx->Sumw2();
    mcD0Dca2Vtx->Sumw2();
  }

  void fill(d0Nt* t)
  {
    if (!isGoodTrack(minPtCut,t->kRPt, t->kREta) || !isGoodTrack(minPtCut,t->pRPt, t->pREta)) return;
    if (t->rM <  anaCuts::massMin || t->rM > anaCuts::massMax) return;
    bool passHft = t->kHft > 0 && t->pHft > 0;
    bool passTpc = t->kTpc>0 && t->pTpc>0;
    if (!passHft || !passTpc) return;

    float weight = t->pt * t->w;
    int ptIndex = getD0PtIndex(t->rPt);
    if(ptIndex < 0) return;

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
   loadHftRatio();
   TFile* fHftRatioCorrection = new TFile("hftRatioCorrection_v1.root");
   gHftRatioCorrection = (TGraphErrors*)fHftRatioCorrection->Get("Graph");

   TFile* fAuAu010Weight = new TFile("AuAu010_weight.root");
   gf1AuAu010Weight = (TF1*)fAuAu010Weight->Get("f1Levy010");

   std::string file = argv[1];
   d0Nt* t = new d0Nt(file);

   TFile* fOut = new TFile("eff.root", "recreate");

   Long64_t nEntries = t->GetEntries();
   cout << "nEntries = " << nEntries << endl;

   std::vector<Hists> hists200MeV;
   std::vector<Hists> hists600MeV;
   std::vector<Hists> hists700MeV;
   std::vector<Hists> hists800MeV;
   std::vector<Hists> hists900MeV;
   std::vector<Hists> hists1000MeV;
   std::vector<Hists> hists1100MeV;
   std::vector<Hists> hists1200MeV;

   TopoHists          topoHists(0.6);

   for(int iCent = 0; iCent < anaCuts::physNCentralities; ++iCent)
   {
     hists200MeV.push_back(Hists(0.2,iCent));
     hists600MeV.push_back(Hists(0.6,iCent));
     hists700MeV.push_back(Hists(0.7,iCent));
     hists800MeV.push_back(Hists(0.8,iCent));
     hists900MeV.push_back(Hists(0.9,iCent));
     hists1000MeV.push_back(Hists(1.0,iCent));
     hists1100MeV.push_back(Hists(1.1,iCent));
     hists1200MeV.push_back(Hists(1.2,iCent));
   }

   for (Long64_t i = 0; i < t->GetEntries(); ++i)
   {
      t->GetEntry(i);

      if (i && i % 1000000 == 0) cout << static_cast<float>(i) / nEntries << endl;

      if (fabs(t->y) > anaCuts::rapidity) continue;

      int d0CentBin = getD0CentIndex(t->cent);
      if(d0CentBin<0) continue;
      hists200MeV[d0CentBin].fill(t);
      hists600MeV[d0CentBin].fill(t);
      hists700MeV[d0CentBin].fill(t);
      hists800MeV[d0CentBin].fill(t);
      hists900MeV[d0CentBin].fill(t);
      hists1000MeV[d0CentBin].fill(t);
      hists1100MeV[d0CentBin].fill(t);
      hists1200MeV[d0CentBin].fill(t);

      topoHists.fill(t);
   }

   for(int iCent = 0; iCent < anaCuts::physNCentralities; ++iCent)
   {
     hists200MeV[iCent].write(fOut);
     hists600MeV[iCent].write(fOut);
     hists700MeV[iCent].write(fOut);
     hists800MeV[iCent].write(fOut);
     hists900MeV[iCent].write(fOut);
     hists1000MeV[iCent].write(fOut);
     hists1100MeV[iCent].write(fOut);
     hists1200MeV[iCent].write(fOut);
   }

   topoHists.write(fOut);
   fOut->Close();
}
