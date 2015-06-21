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
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
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

int main(int argc, char **argv)
{
   d0Nt* t = new d0Nt();

   TFile* fOut = new TFile("eff.root", "recreate");

   TH1D* noCuts010 = new TH1D("noCuts010", "", 60, 0, 12);
   TH1D* noCuts4080 = new TH1D("noCuts4080", "", 60, 0, 12);

   TH1D* topoOnly010 = new TH1D("topoOnly010", "", 60, 0, 12);
   TH1D* topoOnly4080 = new TH1D("topoOnly4080", "", 60, 0, 12);

   TH1D* hftOnly010 = new TH1D("hftOnly010", "", 60, 0, 12);
   TH1D* hftOnly4080 = new TH1D("hftOnly4080", "", 60, 0, 12);

   TH1D* hftTopo010 = new TH1D("hftTopo010", "", 60, 0, 12);
   TH1D* hftTopo4080 = new TH1D("hftTopo4080", "", 60, 0, 12);

   TH1D* hftTpc010 = new TH1D("hftTpc010", "", 60, 0, 12);
   TH1D* hftTpc4080 = new TH1D("hftTpc4080", "", 60, 0, 12);

   Long64_t nEntries = t->GetEntries();
   cout << "nEntries = " << nEntries << endl;

   for (Long64_t i = 0; i < t->GetEntries(); ++i)
   {
      t->GetEntry(i);

      if (i && i % 1000000 == 0) cout << static_cast<float>(i) / nEntries << endl;

      if (fabs(t->y) > anaCuts::rapidity) continue;

      if (!(t->cent >= 7 || t->cent < 4)) continue;

      // fill denominator histograms
      if (t->cent >= 7) noCuts010->Fill(t->rPt);
      else if (t->cent < 4) noCuts4080->Fill(t->rPt);

      if (!isGoodTrack(t->kRPt, t->kREta) || !isGoodTrack(t->pRPt, t->pREta)) continue;
      bool passTopologicalCuts = isGoodPair(t->rPt, t->rY, t->cosTheta, t->pRDca, t->kRDca, t->dca12, t->decayLength, t->dcaD0ToPv);
      bool passHft = t->kHft > 0 && t->pHft > 0;
      bool passTpc = t->reco > 0;

      if (t->cent >= 7) // 0-10
      {
         if (passHft) hftOnly010->Fill(t->rPt);
         if (passTopologicalCuts) topoOnly010->Fill(t->rPt);
         if (passTopologicalCuts && passHft) hftTopo010->Fill(t->rPt);
         if (passTopologicalCuts && passHft && passTpc) hftTpc010->Fill(t->rPt);
      }
      else if (t->cent < 4) // 40-80
      {
         if (passHft) hftOnly4080->Fill(t->rPt);
         if (passTopologicalCuts) topoOnly4080->Fill(t->rPt);
         if (passTopologicalCuts && passHft) hftTopo4080->Fill(t->rPt);
         if (passTopologicalCuts && passHft && passTpc) hftTpc4080->Fill(t->rPt);
      }
   } // end event looping

   noCuts010->Scale(1 / 2.);
   topoOnly010->Scale(1 / 2.);
   hftOnly010->Scale(1 / 2.);
   hftTopo010->Scale(1 / 2.);
   hftTpc010->Scale(1 / 2.);

   noCuts4080->Scale(1 / 4.);
   topoOnly4080->Scale(1. / 4.);
   hftOnly4080->Scale(1. / 4.);
   hftTopo4080->Scale(1 / 4.);
   hftTpc4080->Scale(1. / 4.);

   noCuts010->Sumw2();
   noCuts4080->Sumw2();
   topoOnly010->Sumw2();
   topoOnly4080->Sumw2();
   hftOnly010->Sumw2();
   hftOnly4080->Sumw2();
   hftTopo010->Sumw2();
   hftTopo4080->Sumw2();
   hftTpc010->Sumw2();
   hftTpc4080->Sumw2();

   fOut->cd();
   noCuts010->Write();
   noCuts4080->Write();
   topoOnly010->Write();
   topoOnly4080->Write();
   hftOnly010->Write();
   hftOnly4080->Write();
   hftTopo010->Write();
   hftTopo4080->Write();
   hftTpc010->Write();
   hftTpc4080->Write();
   fOut->Close();
}
