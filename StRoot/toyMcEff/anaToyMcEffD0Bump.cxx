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

#include "d0BumpNt.h"
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
   d0BumpNt* t = new d0BumpNt();

   TFile* fOut = new TFile("D0Bump.hists.root", "recreate");

   TH1F* h763 = new TH1F("h763","",160,0.5,2.1);
   TH1F* h785 = new TH1F("h785","",160,0.5,2.1);
   TH1F* h765 = new TH1F("h765","",160,0.5,2.1);
   TH1F* h764 = new TH1F("h764","",160,0.5,2.1);
   TH1F* h719 = new TH1F("h719","",160,0.5,2.1);

   TH1F* h763MisPid = new TH1F("h763MisPid","",160,0.5,2.1);
   TH1F* h785MisPid = new TH1F("h785MisPid","",160,0.5,2.1);
   TH1F* h765MisPid = new TH1F("h765MisPid","",160,0.5,2.1);
   TH1F* h764MisPid = new TH1F("h764MisPid","",160,0.5,2.1);
   TH1F* h719MisPid = new TH1F("h719MisPid","",160,0.5,2.1);

   TH1F* h763x = new TH1F("h763x","",160,0.5,2.1);
   TH1F* h785x = new TH1F("h785x","",160,0.5,2.1);
   TH1F* h765x = new TH1F("h765x","",160,0.5,2.1);
   TH1F* h764x = new TH1F("h764x","",160,0.5,2.1);
   TH1F* h719x = new TH1F("h719x","",160,0.5,2.1);

   TH1F* h763xMisPid = new TH1F("h763xMisPid","",160,0.5,2.1);
   TH1F* h785xMisPid = new TH1F("h785xMisPid","",160,0.5,2.1);
   TH1F* h765xMisPid = new TH1F("h765xMisPid","",160,0.5,2.1);
   TH1F* h764xMisPid = new TH1F("h764xMisPid","",160,0.5,2.1);
   TH1F* h719xMisPid = new TH1F("h719xMisPid","",160,0.5,2.1);

   Long64_t nEntries = t->GetEntries();
   cout << "nEntries = " << nEntries << endl;

   for (Long64_t i = 0; i < t->GetEntries(); ++i)
   {
      t->GetEntry(i);

      if (i && i % 1000000 == 0) cout << static_cast<float>(i) / nEntries << endl;

      if ( t->cent < 4 || t->cent >6 ) continue;
      if (fabs(t->y) > anaCuts::rapidity || fabs(t->rY)> anaCuts::rapidity) continue;

      if (!isGoodTrack(t->kRPt, t->kREta) || !isGoodTrack(t->pRPt, t->pREta)) continue;
      bool passTopologicalCuts = isGoodPair(t->rPt, t->cosTheta, t->pRDca, t->kRDca, t->dca12, t->decayLength, t->dcaD0ToPv);

      switch( static_cast<int>(t->decayChannel) )
      {
        case 763:
          h763->Fill(t->rM, t->w);
          h763MisPid->Fill(t->misPidM, t->w);
          if(passTopologicalCuts)
          {
            h763x->Fill(t->rM, t->w);
            h763xMisPid->Fill(t->misPidM, t->w);
          }
          break;
        case 785:
          h785->Fill(t->rM, t->w);
          h785MisPid->Fill(t->misPidM, t->w);
          if(passTopologicalCuts)
          {
            h785x->Fill(t->rM, t->w);
            h785xMisPid->Fill(t->misPidM, t->w);
          }
          break;
        case 765:
          h765->Fill(t->rM, t->w);
          h765MisPid->Fill(t->misPidM, t->w);
          if(passTopologicalCuts)
          {
            h765x->Fill(t->rM, t->w);
            h765xMisPid->Fill(t->misPidM, t->w);
          }
          break;
        case 764:
          h764->Fill(t->rM, t->w);
          h764MisPid->Fill(t->misPidM, t->w);
          if(passTopologicalCuts)
          {
            h764x->Fill(t->rM, t->w);
            h764xMisPid->Fill(t->misPidM, t->w);
          }
          break;
        case 719:
          h719->Fill(t->rM, t->w);
          h719MisPid->Fill(t->misPidM, t->w);
          if(passTopologicalCuts)
          {
            h719x->Fill(t->rM, t->w);
            h719xMisPid->Fill(t->misPidM, t->w);
          }
          break;
        default:
          cout<<"Uknown decayChannel!"<<endl;
          break;
      }

   } // end event looping

   fOut->cd();
   h763->Write();
   h785->Write();
   h765->Write();
   h764->Write();
   h719->Write();

   h763MisPid->Write();
   h785MisPid->Write();
   h765MisPid->Write();
   h764MisPid->Write();
   h719MisPid->Write();

   h763x->Write();
   h785x->Write();
   h765x->Write();
   h764x->Write();
   h719x->Write();

   h763xMisPid->Write();
   h785xMisPid->Write();
   h765xMisPid->Write();
   h764xMisPid->Write();
   h719xMisPid->Write();

   fOut->Close();
}
