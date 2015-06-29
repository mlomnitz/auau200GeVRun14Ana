/* *********************************************************************
 *
 *  Analysis code to read D0 Bump *.toyMc.root files.
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
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TNtuple.h"
#endif

#include "d0BumpNt.h"
#include "anaCuts.h"


TGraphErrors* grDoubleMisPID = NULL;
TF1*          fDoubleMisPID = NULL;

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

bool isMisPid(float const pt)
{
   float const eff = pt > 0.5 ? fDoubleMisPID->Eval(pt) : grDoubleMisPID->Eval(pt);
   return gRandom->Rndm() <= eff;
}

class hists
{
public:

   hists(int decayChannel, std::string title): h2Mass(NULL), h2MassMisPid(NULL), h2MassX(NULL), h2MassMisPidX(NULL),
                                               h2KDcaVsPt(NULL), h2PiDcaVsPt(NULL), h2CosThetaVsPt(NULL), h2DcaToPvVsPt(NULL)
   {
      h2Mass = new TH2F(Form("h%i", decayChannel), Form("%s without cuts", title.c_str()), 100, 0, 10, 160, 0.5, 2.1);
      h2MassMisPid = new TH2F(Form("h%iMisPid", decayChannel), Form("%s misPid", title.c_str()), 100, 0, 10, 160, 0.5, 2.1);
      h2MassX = new TH2F(Form("h%ix", decayChannel), Form("%s with cuts", title.c_str()), 100, 0, 10, 160, 0.5, 2.1);
      h2MassMisPidX = new TH2F(Form("h%ixMisPid", decayChannel), Form("%s misPid with wuts", title.c_str()), 100, 0, 10, 160, 0.5, 2.1);

      h2KDcaVsPt = new TH2F(Form("hKDcaVsPt%i",decayChannel),"",100,0,10,2000,0,20000);
      h2PiDcaVsPt = new TH2F(Form("hPiDcaVsPt%i",decayChannel),"",100,0,10,2000,0,20000);
      h2CosThetaVsPt = new TH2F(Form("hCosThetaVsPt%i",decayChannel),"",100,0,10,2000,-1.,1.);
      h2DcaToPvVsPt = new TH2F(Form("hDcaToPv%i",decayChannel),"",100,0,10,2000,0,20000);

      h2Mass->Sumw2();
      h2MassMisPid->Sumw2();
      h2MassX->Sumw2();
      h2MassMisPidX->Sumw2();

      h2KDcaVsPt->Sumw2();
      h2PiDcaVsPt->Sumw2();
      h2CosThetaVsPt->Sumw2();
      h2DcaToPvVsPt->Sumw2();
   }

   void fill(d0BumpNt const* const t, bool passTopologicalCuts, bool misPid)
   {
     h2KDcaVsPt->Fill(t->rPt,t->kRDca);
     h2PiDcaVsPt->Fill(t->rPt,t->pRDca);
     h2CosThetaVsPt->Fill(t->rPt,t->cosTheta);
     h2DcaToPvVsPt->Fill(t->rPt,t->dcaD0ToPv);

     h2Mass->Fill(t->rPt, t->rM, t->pt * t->w);
     if (misPid) h2MassMisPid->Fill(t->rPt, t->misPidM, t->pt * t->w);
     if (passTopologicalCuts)
     {
       h2MassX->Fill(t->rPt, t->rM, t->pt * t->w);
       if (misPid) h2MassMisPidX->Fill(t->rPt, t->misPidM, t->pt * t->w);
     }
   }

   void write(TFile* fOut)
   {
      fOut->cd();
      h2Mass->Write();
      h2MassMisPid->Write();
      h2MassX->Write();
      h2MassMisPidX->Write();
      h2KDcaVsPt->Write();
      h2PiDcaVsPt->Write();
      h2CosThetaVsPt->Write();
      h2DcaToPvVsPt->Write();
   }

private:
   TH2F* h2Mass;
   TH2F* h2MassMisPid;
   TH2F* h2MassX; // X means topological cuts applied.
   TH2F* h2MassMisPidX;
   TH2F* h2KDcaVsPt;
   TH2F* h2PiDcaVsPt;
   TH2F* h2CosThetaVsPt;
   TH2F* h2DcaToPvVsPt;
};

int main(int argc, char **argv)
{
   TH1::AddDirectory(false);
   gRandom->SetSeed();

   cout << "Decay channel  = " << 763 << " : D0 --> K- pi+" << endl;
   cout << "Decay channel  = " << 785 << " : D0 --> K- pi+ pi0" << endl;
   cout << "Decay channel  = " << 765 << " : D0 --> K- rho+ --> K- pi+ pi0" << endl;
   cout << "Decay channel  = " << 764 << " : D0 --> K*- pi+  --> K- pi0 pi+" << endl;
   cout << "Decay channel  = " << 786 << " : D0 --> K- pi+ rho0 --> K- pi+ pi+ pi-" << endl;
   cout << "Decay channel  = " << 719 << " : D+ --> K- pi+ pi+" << endl;

   TFile* fMisPid = new TFile("doubleCountingFit.root");
   grDoubleMisPID = (TGraphErrors*)fMisPid->Get("grD0R")->Clone("grDoubleMisPID");
   fDoubleMisPID = (TF1*)fMisPid->Get("f1")->Clone("fDoubleMisPID");

   d0BumpNt* t = new d0BumpNt();

   TFile* fOut = new TFile("D0Bump.hists.root", "recreate");
   hists hists763(763, "D0 --> K- pi+");
   hists hists785(785, "D0 --> K- pi+ pi0");
   hists hists765(765, "D0 --> K- rho+ --> K- pi+ pi0");
   hists hists764(764, "D0 --> K*- pi+  --> K- pi0 pi+");
   hists hists786(786, "D0 --> K- pi+ rho0 --> K- pi+ pi+ pi-");
   hists hists719(719, "D+ --> K- pi+ pi+");

   Long64_t nEntries = t->GetEntries();
   cout << "nEntries = " << nEntries << endl;

   for (Long64_t i = 0; i < t->GetEntries(); ++i)
   {
      t->GetEntry(i);

      if (i && i % 1000000 == 0) cout << static_cast<float>(i) / nEntries << endl;

      if (t->cent < 4 || t->cent > 6) continue;
      if (fabs(t->y) > anaCuts::rapidity || fabs(t->rY) > anaCuts::rapidity) continue;

      if (!isGoodTrack(t->kRPt, t->kREta) || !isGoodTrack(t->pRPt, t->pREta)) continue;
      bool const passTopologicalCuts = isGoodPair(t->rPt, t->cosTheta, t->pRDca, t->kRDca, t->dca12, t->decayLength, t->dcaD0ToPv);

      bool const misPid = isMisPid(t->rPt);

      switch (static_cast<int>(t->decayChannel))
      {
         case 763:
            hists763.fill(t,passTopologicalCuts,misPid);
            break;
         case 785:
            hists785.fill(t,passTopologicalCuts,misPid);
            break;
         case 765:
            hists765.fill(t,passTopologicalCuts,misPid);
            break;
         case 764:
            hists764.fill(t,passTopologicalCuts,misPid);
            break;
         case 786:
            hists786.fill(t,passTopologicalCuts,misPid);
            break;
         case 719:
            hists719.fill(t,passTopologicalCuts,misPid);
            break;
         default:
            cout << "Uknown decayChannel!" << endl;
            break;
      }

   } // end event looping

   hists763.write(fOut);
   hists785.write(fOut);
   hists765.write(fOut);
   hists764.write(fOut);
   hists786.write(fOut);
   hists719.write(fOut);
}
