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
#include "TH2D.h"
#include "TH3F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"
#endif

#include "d0ZeroDecayLength.h"

const Int_t nParticles = 2;
const Int_t nCent = 9;

//----------------
const Int_t nEtasHftRatio = 10;
const Double_t EtaEdgeHftRatio[nEtasHftRatio + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
const Int_t nVzsHftRatio = 6;
const Double_t VzEdgeHftRatio[nVzsHftRatio + 1] = { -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0};
const Int_t nPtBinsHftRatio = 35;
const Double_t ptEdgeHftRatio[nPtBinsHftRatio + 1] = { 0.0, 0.2, 0.4,  0.6,  0.8,
                                       1.0, 1.2, 1.4,  1.6,  1.8,
                                       2.0, 2.2, 2.4,  2.6,  2.8,
                                       3.0, 3.2, 3.4,  3.6,  3.8,
                                       4.0, 4.2, 4.4,  4.6,  4.8,
                                       5.0, 5.4, 5.8,  6.2,  6.6,
                                       7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
const Int_t nPhisHftRatio = 30;
const Double_t PhiEdgeHftRatio[nPhisHftRatio + 1] = { -3.14159, -2.93215, -2.72271, -2.51327, -2.30383, -2.0944, -1.88496, -1.67552, -1.46608, -1.25664, -1.0472, -0.837758, -0.628319, -0.418879, -0.20944, 0.0, 0.20944, 0.418879, 0.628319, 0.837758, 1.0472, 1.25664, 1.46608, 1.67552, 1.88496, 2.0944, 2.30383, 2.51327, 2.72271, 2.93215, 3.14159};

//----------------
int const nVzs = 4;
float const VzEdge[nVzs + 1] = { -6., -3., 0, 3., 6.};

int const nPtBins = 28;
const Double_t ptEdge[nPtBins + 1] = { 0. , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 ,
      1. , 1.2 , 1.4 , 1.6 , 1.8 , 2.  , 2.2 , 2.4 , 2.6 , 2.8 , 
      3. , 3.5 , 4.  , 4.5 , 5. , 6. , 8.0 , 10. , 12.0};

int const nEtas = 5;
float const EtaEdge[nEtas + 1] = { -1.0, -0.6, -0.2, 0.2, 0.6, 1.0};

int const nDcas = 140;
float const DcaEdge[nDcas + 1] =
{ -1 ,-0.97 ,-0.94 ,-0.91 ,-0.88 ,-0.85 ,-0.82 ,-0.79 ,-0.76 ,-0.73 ,-0.7 ,-0.67 ,-0.64 ,-0.61 ,-0.58 ,-0.55 ,-0.52 ,-0.49 ,-0.46 ,-0.43 ,-0.4 ,-0.37 ,-0.34 ,-0.31 ,-0.28 ,-0.25 ,-0.22 ,-0.19 ,-0.16 ,-0.13 ,// 30 //0.03cm/perbin
  -0.1 ,-0.0975 ,-0.095 ,-0.0925 ,-0.09 ,-0.0875 ,-0.085 ,-0.0825 ,-0.08 ,-0.0775 ,-0.075 ,-0.0725 ,-0.07 ,-0.0675 ,-0.065 ,-0.0625 ,-0.06 ,-0.0575 ,-0.055 ,-0.0525 ,-0.05 ,-0.0475 ,-0.045 ,-0.0425 ,-0.04 ,-0.0375 ,-0.035 ,-0.0325 ,-0.03 ,-0.0275 ,-0.025 ,-0.0225 ,-0.02 ,-0.0175 ,-0.015 ,-0.0125 ,-0.01 ,-0.0075 ,-0.005 ,-0.0025 ,0 ,0.0025 ,0.005 ,0.0075 ,0.01 ,0.0125 ,0.015 ,0.0175 ,0.02 ,0.0225 ,0.025 ,0.0275 ,0.03 ,0.0325 ,0.035 ,0.0375 ,0.04 ,0.0425 ,0.045 ,0.0475 ,0.05 ,0.0525 ,0.055 ,0.0575 ,0.06 ,0.0625 ,0.065 ,0.0675 ,0.07 ,0.0725 ,0.075 ,0.0775 ,0.08 ,0.0825 ,0.085 ,0.0875 ,0.09 ,0.0925 ,0.095 ,0.0975 ,0.1 ,//80 //0.0025/perbin
  0.13 ,0.16 ,0.19 ,0.22 ,0.25 ,0.28 ,0.31 ,0.34 ,0.37 ,0.4 ,0.43 ,0.46 ,0.49 ,0.52 ,0.55 ,0.58 ,0.61 ,0.64 ,0.67 ,0.7 ,0.73 ,0.76 ,0.79 ,0.82 ,0.85 ,0.88 ,0.91 ,0.94 ,0.97 ,1 //30 //0.03cm/perbin
};

TH1D* hHftDenominator[nParticles][nEtasHftRatio][nVzsHftRatio][nPhisHftRatio][nCent];
TH1D* hHftRatio[nParticles][nEtasHftRatio][nVzsHftRatio][nPhisHftRatio][nCent];
TH1D* h1Dca[nParticles][nEtas][nVzs][nCent][nPtBins];
TH1D* h1DcaZ[nParticles][nEtas][nVzs][nCent][nPtBins];
TH1D* h1DcaXY[nParticles][nEtas][nVzs][nCent][nPtBins];

void bookHistograms()
{
  for (int iParticle = 0; iParticle < nParticles; ++iParticle)
  {
    for (int iEta = 0; iEta < nEtas; ++iEta)
    {
      for (int iVz = 0; iVz < nVzs; ++iVz)
      {
        for (int ii = 0; ii < nCent; ++ii)
        {
          for (int jj = 0; jj < nPtBins; ++jj)
          {
            TString name = Form("simDcaPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj);
            TString title = Form("simDcaPtCentPartEtaVz_%i_%1.1f_%1.1f_%i_%1.2f_pT_%1.2f", iParticle, EtaEdge[iEta], VzEdge[iVz], ii, ptEdge[jj], ptEdge[jj+1]);
            TString nameZ = Form("simDcaZPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj);
            TString titleZ = Form("simDcaZPtCentPartEtaVz_%i_%1.1f_%1.1f_%i_%1.2f_pT_%1.2f", iParticle, EtaEdge[iEta], VzEdge[iVz], ii, ptEdge[jj], ptEdge[jj+1]);
            TString nameXY = Form("simDcaXyPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj);
            TString titleXY = Form("simDcaXYPtCentPartEtaVz_%i_%1.1f_%1.1f_%i_%1.2f_pT_%1.2f", iParticle, EtaEdge[iEta], VzEdge[iVz], ii, ptEdge[jj], ptEdge[jj+1]);
            h1Dca[iParticle][iEta][iVz][ii][jj] = new TH1D(name.Data(),title.Data(),nDcas,DcaEdge);
            h1DcaZ[iParticle][iEta][iVz][ii][jj] = new TH1D(nameZ.Data(),titleZ.Data(),nDcas,DcaEdge);
            h1DcaXY[iParticle][iEta][iVz][ii][jj] = new TH1D(nameXY.Data(),titleXY.Data(),nDcas,DcaEdge);
          }
        }
      }
    }

    for(int iCent = 0; iCent < nCent; ++iCent)
    {
      for (int iEta = 0; iEta < nEtasHftRatio; ++iEta)
      {
        for (int iVz = 0; iVz < nVzsHftRatio; ++iVz)
        {
          for(int iPhi = 0; iPhi < nPhisHftRatio; ++iPhi)
          {
            TString name = Form("mh1HFT1PtCentPartEtaVzPhiRatio_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iPhi, iCent);
            TString title = Form("mh1HFT1PtCentPartEtaVzPhiRatio_%i_%1.1f_%1.1f_%1.1f_%i", iParticle, EtaEdgeHftRatio[iEta], VzEdge[iVz], PhiEdgeHftRatio[iPhi], iCent);
            hHftRatio[iParticle][iEta][iVz][iPhi][iCent] = new TH1D(name.Data(),title.Data(),nPtBinsHftRatio,ptEdgeHftRatio);

            name = Form("mh1HFT1PtCentPartEtaVzPhiDenominator_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iPhi, iCent);
            title = Form("mh1HFT1PtCentPartEtaVzPhiDenominator_%i_%1.1f_%1.1f_%1.1f_%i", iParticle, EtaEdgeHftRatio[iEta], VzEdge[iVz], PhiEdgeHftRatio[iPhi], iCent);
            hHftDenominator[iParticle][iEta][iVz][iPhi][iCent] = new TH1D(name.Data(),title.Data(),nPtBinsHftRatio,ptEdgeHftRatio);
          }
        }
      }
    }
  }
}

void write(TFile* f)
{
  f->cd();
  for (int iParticle = 0; iParticle < nParticles; ++iParticle)
  {
    for (int iEta = 0; iEta < nEtas; ++iEta)
    {
      for (int iVz = 0; iVz < nVzs; ++iVz)
      {
        for (int ii = 0; ii < nCent; ++ii)
        {
          for (int jj = 0; jj < nPtBins; ++jj)
          {
            h1Dca[iParticle][iEta][iVz][ii][jj]->Write();
            h1DcaZ[iParticle][iEta][iVz][ii][jj]->Write();
            h1DcaXY[iParticle][iEta][iVz][ii][jj]->Write();
          }
        }
      }
    }

    for(int iCent = 0; iCent < nCent; ++iCent)
    {
      for (int iEta = 0; iEta < nEtasHftRatio; ++iEta)
      {
        for (int iVz = 0; iVz < nVzsHftRatio; ++iVz)
        {
          for(int iPhi = 0; iPhi < nPhisHftRatio; ++iPhi)
          {
            hHftRatio[iParticle][iEta][iVz][iPhi][iCent]->Divide(hHftDenominator[iParticle][iEta][iVz][iPhi][iCent]);
            hHftRatio[iParticle][iEta][iVz][iPhi][iCent]->Write();
            hHftDenominator[iParticle][iEta][iVz][iPhi][iCent]->Write();
          }
        }
      }
    }
  }
}

using namespace std;

int getPtIndex(double pT)
{
   for (int i = 0; i < nPtBins; i++)
   {
      if ((pT >= ptEdge[i]) && (pT < ptEdge[i + 1]))
         return i;
   }
   return nPtBins - 1 ;
}

int getEtaIndex(double Eta)
{
   for (int i = 0; i < nEtas; i++)
   {
      if ((Eta >= EtaEdge[i]) && (Eta < EtaEdge[i + 1]))
         return i;
   }
   return nEtas - 1 ;
}

int getVzIndex(double Vz)
{
   for (int i = 0; i < nVzs; i++)
   {
      if ((Vz >= VzEdge[i]) && (Vz < VzEdge[i + 1]))
         return i;
   }
   return nVzs - 1 ;
}

int getEtaIndexHftRatio(double Eta)
{
   for (int i = 0; i < nEtasHftRatio; i++)
   {
      if ((Eta >= EtaEdgeHftRatio[i]) && (Eta < EtaEdgeHftRatio[i + 1]))
         return i;
   }
   return nEtasHftRatio - 1 ;
}

int getVzIndexHftRatio(double Vz)
{
   for (int i = 0; i < nVzsHftRatio; i++)
   {
      if ((Vz >= VzEdgeHftRatio[i]) && (Vz < VzEdgeHftRatio[i + 1]))
         return i;
   }
   return nVzsHftRatio - 1 ;
}

int getPhiIndexHftRatio(double Phi)
{
   for (int i = 0; i < nPhisHftRatio; i++)
   {
      if ((Phi >= PhiEdgeHftRatio[i]) && (Phi < PhiEdgeHftRatio[i + 1]))
         return i;
   }
   return nPhisHftRatio - 1 ;
}

bool isGoodTrack(float const eta)
{
   return fabs(eta) < 1;
}

int main(int argc, char **argv)
{
   d0Nt* t = new d0Nt();

   TFile* fOut = new TFile("fastSim.3DcaBinning.2DDcaSampling.Dca.root", "recreate");

   // TH1::AddDirectory(0);
   bookHistograms();

   Long64_t nEntries = t->GetEntries();
   cout << "nEntries = " << nEntries << endl;

   for (Long64_t i = 0; i < t->GetEntries(); ++i)
   {
      t->GetEntry(i);

      if (i && i % 1000000 == 0) cout << static_cast<float>(i) / nEntries << endl;

      int const vzIdx = getVzIndex(t->vz/1.e4);
      int const vzIdxHftRatio = getVzIndexHftRatio(t->vz/1.e4);

      if(isGoodTrack(t->pREta))
      {
        h1Dca[0][getEtaIndex(t->pREta)][vzIdx][(int)t->cent][getPtIndex(t->pRPt)]->Fill(t->pRSDca/1.e4);
        h1DcaZ[0][getEtaIndex(t->pREta)][vzIdx][(int)t->cent][getPtIndex(t->pRPt)]->Fill(t->pRDcaZ/1.e4);
        h1DcaXY[0][getEtaIndex(t->pREta)][vzIdx][(int)t->cent][getPtIndex(t->pRPt)]->Fill(t->pRDcaXY/1.e4);
        hHftDenominator[0][getEtaIndexHftRatio(t->pREta)][vzIdxHftRatio][getPhiIndexHftRatio(t->pRPhi)][(int)t->cent]->Fill(t->pRPt);
        if(t->pHft > 1) hHftRatio[0][getEtaIndexHftRatio(t->pREta)][vzIdxHftRatio][getPhiIndexHftRatio(t->pRPhi)][(int)t->cent]->Fill(t->pRPt);
      }

      if(isGoodTrack(t->kREta))
      {
        h1Dca[1][getEtaIndex(t->kREta)][vzIdx][(int)t->cent][getPtIndex(t->kRPt)]->Fill(t->kRSDca/1.e4);
        h1DcaZ[1][getEtaIndex(t->kREta)][vzIdx][(int)t->cent][getPtIndex(t->kRPt)]->Fill(t->kRDcaZ/1.e4);
        h1DcaXY[1][getEtaIndex(t->kREta)][vzIdx][(int)t->cent][getPtIndex(t->kRPt)]->Fill(t->kRDcaXY/1.e4);
        hHftDenominator[1][getEtaIndexHftRatio(t->kREta)][vzIdxHftRatio][getPhiIndexHftRatio(t->kRPhi)][(int)t->cent]->Fill(t->kRPt);
        if(t->kHft > 0) hHftRatio[1][getEtaIndexHftRatio(t->kREta)][vzIdxHftRatio][getPhiIndexHftRatio(t->kRPhi)][(int)t->cent]->Fill(t->kRPt);
      }
   }

   write(fOut);
   fOut->Close();
}
