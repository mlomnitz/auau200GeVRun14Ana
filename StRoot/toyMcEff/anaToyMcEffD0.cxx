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
  for(int i=0; i < anaCuts::nPtBins; i++) 
  {
    if( (pt >= anaCuts::PtEdge[i]) && (pt < anaCuts::PtEdge[i+1]) )
     bin = i;
  }
  return bin;
}

bool isGoodPair(float const pt, float const cosTheta, float const pDca, float const kDca,
                float const dca12, float const decayLength, float const dcaV0ToPv)
{
  int tmpIndex = getD0PtIndex(pt);

  return cosTheta > anaCuts::cosTheta[tmpIndex] &&
    pDca > anaCuts::pDca[tmpIndex] && kDca > anaCuts::kDca[tmpIndex] &&
    dca12 < anaCuts::dcaDaughters[tmpIndex] &&
    decayLength > anaCuts::decayLength[tmpIndex] &&
    dcaV0ToPv < anaCuts::dcaV0ToPv[tmpIndex];
}

int main(int argc, char **argv)
{
   d0Nt* t = new d0Nt();

   TFile* fOut = new TFile("eff.root","recreate");
   TH1D* topoOnly010= new TH1D("topoOnly010","",60,0,12);
   TH1D* topoOnly4080= new TH1D("topoOnly4080","",60,0,12);
   TH1D* hftOnly010= new TH1D("hftOnly010","",60,0,12);
   TH1D* hftOnly4080= new TH1D("hftOnly4080","",60,0,12);
   TH1D* hftTpc010= new TH1D("hftTpc010","",60,0,12);
   TH1D* hftTpc4080= new TH1D("hftTpc4080","",60,0,12);

   for (Long64_t i = 0; i < t->GetEntries(); ++i)
   {
      t->GetEntry(i);

      if(!(t->cent>=7 || t->cent<4)) continue;

      bool passTopologicalCuts = isGoodPair(t->rPt,t->cosTheta,t->pRDca,t->kRDca,t->dca12,t->decayLength,t->dcaD0ToPv);

      if(!passTopologicalCuts) continue;

      if(t->cent>=7) // 0-10
      {
        topoOnly010->Fill(t->rPt);
        if(t->kHft>0 && t->pHft>0) 
        {
          hftOnly010->Fill(t->rPt);
          if(t->reco) hftOnly010->Fill(t->rPt);
        }
      }
      else if(t->cent<4) // 40-80
      {
        topoOnly4080->Fill(t->rPt);
        if(t->kHft>0 && t->pHft>0) 
        {
          hftOnly4080->Fill(t->rPt);
          if(t->reco) hftOnly4080->Fill(t->rPt);
        }
      }
   } // end event looping

   topoOnly010->Scale(1/2.);
   topoOnly4080->Scale(1./4.);
   hftOnly010->Scale(1/2.);
   hftOnly4080->Scale(1./4.);
   hftTpc010->Scale(1/2.);
   hftTpc4080->Scale(1./4.);

   topoOnly010->Sumw2();
   topoOnly4080->Sumw2();
   hftOnly010->Sumw2();
   hftOnly4080->Sumw2();
   hftTpc010->Sumw2();
   hftTpc4080->Sumw2();

   fOut->cd();
   topoOnly010->Write();
   topoOnly4080->Write();
   hftOnly010->Write();
   hftOnly4080->Write();
   hftTpc010->Write();
   hftTpc4080->Write();
   fOut->Close();
}
