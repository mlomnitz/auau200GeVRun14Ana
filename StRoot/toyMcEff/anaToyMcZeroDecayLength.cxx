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
// const Int_t nEtas = 10;
// const Double_t EtaEdge[nEtas + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
int const nEtas = 5;
float const EtaEdge[nEtas + 1] = { -1.0, -0.6, -0.2, 0.2, 0.6, 1.0};
// const Int_t nVzs = 6;
// const Double_t VzEdge[nVzs + 1] = { -6.0e4, -4.0e4, -2.0e4, 0.0, 2.0e4, 4.0e4, 6.0e4};
int const nVzs = 4;
float const VzEdge[nVzs + 1] = { -6., -3., 0, 3., 6.};
// const Int_t nPtBins = 35;
// const Double_t ptEdge[nPtBins + 1] = { 0.0, 0.2, 0.4,  0.6,  0.8,
//                                        1.0, 1.2, 1.4,  1.6,  1.8,
//                                        2.0, 2.2, 2.4,  2.6,  2.8,
//                                        3.0, 3.2, 3.4,  3.6,  3.8,
//                                        4.0, 4.2, 4.4,  4.6,  4.8,
//                                        5.0, 5.4, 5.8,  6.2,  6.6,
//                                        7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

int const nPtBins = 28;
const Double_t ptEdge[nPtBins + 1] = { 0. , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 ,
      1. , 1.2 , 1.4 , 1.6 , 1.8 , 2.  , 2.2 , 2.4 , 2.6 , 2.8 , 
      3. , 3.5 , 4.  , 4.5 , 5. , 6. , 8.0 , 10. , 12.0};

// int const nDcas=380;
// float const DcaEdge[nDcas+1]= {
// -1.0 ,-0.99 ,-0.98 ,-0.97 ,-0.96 ,-0.95 ,-0.94 ,-0.93 ,-0.92 ,-0.91 ,-0.9 ,-0.89 ,-0.88 ,-0.87 ,-0.86 ,-0.85 ,-0.84 ,-0.83 ,-0.82 ,-0.81 ,-0.8 ,-0.79 ,-0.78 ,-0.77 ,-0.76 ,-0.75 ,-0.74 ,-0.73 ,-0.72 ,-0.71 ,-0.7 ,-0.69 ,-0.68 ,-0.67 ,-0.66 ,-0.65 ,-0.64 ,-0.63 ,-0.62 ,-0.61 ,-0.6 ,-0.59 ,-0.58 ,-0.57 ,-0.56 ,-0.55 ,-0.54 ,-0.53 ,-0.52 ,-0.51 ,-0.5 ,-0.49 ,-0.48 ,-0.47 ,-0.46 ,-0.45 ,-0.44 ,-0.43 ,-0.42 ,-0.41 ,-0.4 ,-0.39 ,-0.38 ,-0.37 ,-0.36 ,-0.35 ,-0.34 ,-0.33 ,-0.32 ,-0.31 ,-0.3 ,-0.29 ,-0.28 ,-0.27 ,-0.26 ,-0.25 ,-0.24 ,-0.23 ,-0.22 ,-0.21 ,-0.2 ,-0.19 ,-0.18 ,-0.17 ,-0.16 ,-0.15 ,-0.14 ,-0.13 ,-0.12 ,-0.11 ,// 90
// -0.1 ,-0.099 ,-0.098 ,-0.097 ,-0.096 ,-0.095 ,-0.094 ,-0.093 ,-0.092 ,-0.091 ,-0.09 ,-0.089 ,-0.088 ,-0.087 ,-0.086 ,-0.085 ,-0.084 ,-0.083 ,-0.082 ,-0.081 ,-0.08 ,-0.079 ,-0.078 ,-0.077 ,-0.076 ,-0.075 ,-0.074 ,-0.073 ,-0.072 ,-0.071 ,-0.07 ,-0.069 ,-0.068 ,-0.067 ,-0.066 ,-0.065 ,-0.064 ,-0.063 ,-0.062 ,-0.061 ,-0.06 ,-0.059 ,-0.058 ,-0.057 ,-0.056 ,-0.055 ,-0.054 ,-0.053 ,-0.052 ,-0.051 ,-0.05 ,-0.049 ,-0.048 ,-0.047 ,-0.046 ,-0.045 ,-0.044 ,-0.043 ,-0.042 ,-0.041 ,-0.04 ,-0.039 ,-0.038 ,-0.037 ,-0.036 ,-0.035 ,-0.034 ,-0.033 ,-0.032 ,-0.031 ,-0.03 ,-0.029 ,-0.028 ,-0.027 ,-0.026 ,-0.025 ,-0.024 ,-0.023 ,-0.022 ,-0.021 ,-0.02 ,-0.019 ,-0.018 ,-0.017 ,-0.016 ,-0.015 ,-0.014 ,-0.013 ,-0.012 ,-0.011 ,-0.01 ,-0.009 ,-0.008 ,-0.007 ,-0.006 ,-0.005 ,-0.004 ,-0.003 ,-0.002 ,-0.001 ,0.0 ,0.001 ,0.002 ,0.003 ,0.004 ,0.005 ,0.006 ,0.007 ,0.008 ,0.009 ,0.01 ,0.011 ,0.012 ,0.013 ,0.014 ,0.015 ,0.016 ,0.017 ,0.018 ,0.019 ,0.02 ,0.021 ,0.022 ,0.023 ,0.024 ,0.025 ,0.026 ,0.027 ,0.028 ,0.029 ,0.03 ,0.031 ,0.032 ,0.033 ,0.034 ,0.035 ,0.036 ,0.037 ,0.038 ,0.039 ,0.04 ,0.041 ,0.042 ,0.043 ,0.044 ,0.045 ,0.046 ,0.047 ,0.048 ,0.049 ,0.05 ,0.051 ,0.052 ,0.053 ,0.054 ,0.055 ,0.056 ,0.057 ,0.058 ,0.059 ,0.06 ,0.061 ,0.062 ,0.063 ,0.064 ,0.065 ,0.066 ,0.067 ,0.068 ,0.069 ,0.07 ,0.071 ,0.072 ,0.073 ,0.074 ,0.075 ,0.076 ,0.077 ,0.078 ,0.079 ,0.08 ,0.081 ,0.082 ,0.083 ,0.084 ,0.085 ,0.086 ,0.087 ,0.088 ,0.089 ,0.09 ,0.091 ,0.092 ,0.093 ,0.094 ,0.095 ,0.096 ,0.097 ,0.098 ,0.099 ,0.1 ,//201
// 0.11 ,0.12 ,0.13 ,0.14 ,0.15 ,0.16 ,0.17 ,0.18 ,0.19 ,0.2 ,0.21 ,0.22 ,0.23 ,0.24 ,0.25 ,0.26 ,0.27 ,0.28 ,0.29 ,0.3 ,0.31 ,0.32 ,0.33 ,0.34 ,0.35 ,0.36 ,0.37 ,0.38 ,0.39 ,0.4 ,0.41 ,0.42 ,0.43 ,0.44 ,0.45 ,0.46 ,0.47 ,0.48 ,0.49 ,0.5 ,0.51 ,0.52 ,0.53 ,0.54 ,0.55 ,0.56 ,0.57 ,0.58 ,0.59 ,0.6 ,0.61 ,0.62 ,0.63 ,0.64 ,0.65 ,0.66 ,0.67 ,0.68 ,0.69 ,0.7 ,0.71 ,0.72 ,0.73 ,0.74 ,0.75 ,0.76 ,0.77 ,0.78 ,0.79 ,0.8 ,0.81 ,0.82 ,0.83 ,0.84 ,0.85 ,0.86 ,0.87 ,0.88 ,0.89 ,0.9 ,0.91 ,0.92 ,0.93 ,0.94 ,0.95 ,0.96 ,0.97 ,0.98 ,0.99 ,1.0 //90
                                // };

int const nDcas = 140;
float const DcaEdge[nDcas + 1] =
{ -1 ,-0.97 ,-0.94 ,-0.91 ,-0.88 ,-0.85 ,-0.82 ,-0.79 ,-0.76 ,-0.73 ,-0.7 ,-0.67 ,-0.64 ,-0.61 ,-0.58 ,-0.55 ,-0.52 ,-0.49 ,-0.46 ,-0.43 ,-0.4 ,-0.37 ,-0.34 ,-0.31 ,-0.28 ,-0.25 ,-0.22 ,-0.19 ,-0.16 ,-0.13 ,// 30 //0.03cm/perbin
  -0.1 ,-0.0975 ,-0.095 ,-0.0925 ,-0.09 ,-0.0875 ,-0.085 ,-0.0825 ,-0.08 ,-0.0775 ,-0.075 ,-0.0725 ,-0.07 ,-0.0675 ,-0.065 ,-0.0625 ,-0.06 ,-0.0575 ,-0.055 ,-0.0525 ,-0.05 ,-0.0475 ,-0.045 ,-0.0425 ,-0.04 ,-0.0375 ,-0.035 ,-0.0325 ,-0.03 ,-0.0275 ,-0.025 ,-0.0225 ,-0.02 ,-0.0175 ,-0.015 ,-0.0125 ,-0.01 ,-0.0075 ,-0.005 ,-0.0025 ,0 ,0.0025 ,0.005 ,0.0075 ,0.01 ,0.0125 ,0.015 ,0.0175 ,0.02 ,0.0225 ,0.025 ,0.0275 ,0.03 ,0.0325 ,0.035 ,0.0375 ,0.04 ,0.0425 ,0.045 ,0.0475 ,0.05 ,0.0525 ,0.055 ,0.0575 ,0.06 ,0.0625 ,0.065 ,0.0675 ,0.07 ,0.0725 ,0.075 ,0.0775 ,0.08 ,0.0825 ,0.085 ,0.0875 ,0.09 ,0.0925 ,0.095 ,0.0975 ,0.1 ,//80 //0.0025/perbin
  0.13 ,0.16 ,0.19 ,0.22 ,0.25 ,0.28 ,0.31 ,0.34 ,0.37 ,0.4 ,0.43 ,0.46 ,0.49 ,0.52 ,0.55 ,0.58 ,0.61 ,0.64 ,0.67 ,0.7 ,0.73 ,0.76 ,0.79 ,0.82 ,0.85 ,0.88 ,0.91 ,0.94 ,0.97 ,1 //30 //0.03cm/perbin
};

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

bool isGoodTrack(float const eta)
{
   return fabs(eta) < 1;
}

int main(int argc, char **argv)
{
   d0Nt* t = new d0Nt();

   TFile* fOut = new TFile("fastSim.3DcaBinning.Dca.root", "recreate");

   // TH1::AddDirectory(0);
   bookHistograms();

   Long64_t nEntries = t->GetEntries();
   cout << "nEntries = " << nEntries << endl;

   for (Long64_t i = 0; i < t->GetEntries(); ++i)
   {
      t->GetEntry(i);

      if (i && i % 1000000 == 0) cout << static_cast<float>(i) / nEntries << endl;

      int const vzIdx = getVzIndex(t->vz);

      if(isGoodTrack(t->pREta) && t->pHft>0)
      {
        h1Dca[0][getEtaIndex(t->pREta)][vzIdx][(int)t->cent][getPtIndex(t->pRPt)]->Fill(t->pRSDca);
        h1DcaZ[0][getEtaIndex(t->pREta)][vzIdx][(int)t->cent][getPtIndex(t->pRPt)]->Fill(t->pRDcaZ);
        h1DcaXY[0][getEtaIndex(t->pREta)][vzIdx][(int)t->cent][getPtIndex(t->pRPt)]->Fill(t->pRDcaXY);
      }

      if(isGoodTrack(t->kREta) && t->kHft>0)
      {
        h1Dca[1][getEtaIndex(t->kREta)][vzIdx][(int)t->cent][getPtIndex(t->kRPt)]->Fill(t->kRSDca);
        h1DcaZ[1][getEtaIndex(t->kREta)][vzIdx][(int)t->cent][getPtIndex(t->kRPt)]->Fill(t->kRDcaZ);
        h1DcaXY[1][getEtaIndex(t->kREta)][vzIdx][(int)t->cent][getPtIndex(t->kRPt)]->Fill(t->kRDcaXY);
      }
   }

   write(fOut);
   fOut->Close();
}
