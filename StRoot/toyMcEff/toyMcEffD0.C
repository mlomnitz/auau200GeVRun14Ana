/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for D0 decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for D0 --> Kpi
 *
 *  Authors:
 *            Guannan Xie (guannanxie@lbl.gov)
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu (hqiu@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
*/

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TClonesArray.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TRandom3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"

using namespace std;

void setDecayChannels(int const mdme);
void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TVector3 v00);
void getKinematics(TLorentzVector& b, double const mass);
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos);
TVector3 smearPosData(int iParticleIndex, double vz, int cent, TLorentzVector const& rMom, TVector3 const& pos);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
TVector3 getVertex(int centrality);
bool matchHft(int iParticleIndex, double vz, int cent, TLorentzVector const& mom);
bool reconstructD0(int const centrality, TLorentzVector const& mom);
void bookObjects();
void write();
int getPtIndex(double);
int getEtaIndex(double);
int getVzIndex(double);
int getPhiIndex(double);

TPythia6Decayer* pydecay;
TNtuple* nt;
TFile* result;

TF1* fKaonMomResolution = NULL;
TF1* fPionMomResolution = NULL;
TF1* fWeightFunction = NULL;
TGraph* grEff[3];
const Int_t nParticles = 2;
const Int_t nEtas = 10;
const Int_t nVzs = 6;
const Int_t nPhis = 30;
const Int_t nCent = 9;
const Int_t nPtBins = 35;
const Double_t EtaEdge[nEtas + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
const Double_t VzEdge[nVzs + 1] = { -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0};
const Double_t PhiEdge[nPhis + 1] = { -3.14159, -2.93215, -2.72271, -2.51327, -2.30383, -2.0944, -1.88496, -1.67552, -1.46608, -1.25664, -1.0472, -0.837758, -0.628319, -0.418879, -0.20944, 0.0, 0.20944, 0.418879, 0.628319, 0.837758, 1.0472, 1.25664, 1.46608, 1.67552, 1.88496, 2.0944, 2.30383, 2.51327, 2.72271, 2.93215, 3.14159};
const Double_t ptEdge[nPtBins + 1] = { 0.0, 0.2, 0.4,  0.6,  0.8,
                                       1.0, 1.2, 1.4,  1.6,  1.8,
                                       2.0, 2.2, 2.4,  2.6,  2.8,
                                       3.0, 3.2, 3.4,  3.6,  3.8,
                                       4.0, 4.2, 4.4,  4.6,  4.8,
                                       5.0, 5.4, 5.8,  6.2,  6.6,
                                       7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

TH1D* hHftRatio[nCent];
TH1D* h1DcaZ[nCent][nPtBins];
TH1D* h1DcaXY[nCent][nPtBins];
TH1D* h1Vz[nCent];

TH1D* hHftRatio1[nParticles][nEtas][nVzs][nCent];
TH1D* h1DcaZ1[nParticles][nEtas][nVzs][nCent][nPtBins];
TH1D* h1DcaXY1[nParticles][nEtas][nVzs][nCent][nPtBins];

string outFileName = "D0.toyMc.root";
std::pair<int, int> const decayChannels(747, 807);
std::pair<float, float> const momentumRange(0.3, 12);

float const acceptanceRapidity = 1.0;
float const sigmaPos0 = 15.2;
float const pxlLayer1Thickness = 0.00486;
float const sigmaVertexCent[nCent] = {31., 18.1, 12.8, 9.3, 7.2, 5.9, 5., 4.6, 4.};

//============== main  program ==================
void toyMcEffD0(int npart = 100)
{
   gRandom->SetSeed();
   bookObjects();

   pydecay = TPythia6Decayer::Instance();
   pydecay->Init();

   setDecayChannels(763); // D0 --> Kpi

   TLorentzVector* b_d = new TLorentzVector;
   TClonesArray ptl("TParticle", 10);
   for (int ipart = 0; ipart < npart; ipart++)
   {
      if (!(ipart % 100000))
         cout << "____________ ipart = " << ipart << " ________________" << endl;

      getKinematics(*b_d, M_D_0);

      decayAndFill(421, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
      decayAndFill(-421, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
   }

   write();
}

void setDecayChannels(int const mdme)
{
   for (int idc = decayChannels.first; idc < decayChannels.second + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
   TPythia6::Instance()->SetMDME(mdme, 1, 1);
}

void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters)
{
   pydecay->Decay(kf, b);
   pydecay->ImportParticles(&daughters);

   TLorentzVector kMom;
   TLorentzVector pMom;
   TVector3 v00;

   int nTrk = daughters.GetEntriesFast();
   for (int iTrk = 0; iTrk < nTrk; ++iTrk)
   {
      TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

      switch (abs(ptl0->GetPdgCode()))
      {
         case 321:
            ptl0->Momentum(kMom);
            v00.SetXYZ(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.); // converted to μm
            break;
         case 211:
            ptl0->Momentum(pMom);
            break;
         default:
            break;
      }
   }
   daughters.Clear();

   fill(kf, b, weight, kMom, pMom, v00);
}

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& pMom, TVector3 v00)
{
   int const centrality = floor(nCent * gRandom->Rndm());

   TVector3 const vertex = getVertex(centrality);
   // smear primary vertex
   // float const sigmaVertex = sigmaVertexCent[cent];
   // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));

   v00 += vertex;

   // smear momentum
   TLorentzVector const kRMom = smearMom(kMom, fKaonMomResolution);
   TLorentzVector const pRMom = smearMom(pMom, fPionMomResolution);

   // smear position
   TVector3 const kRPos = smearPosData(1, vertex.z(), centrality, kRMom, v00);
   TVector3 const pRPos = smearPosData(0, vertex.z(), centrality, pRMom, v00);
   // TVector3 const kRPos = smearPos(kMom, kRMom, v00);
   // TVector3 const pRPos = smearPos(pMom, pRMom, v00);

   // reconstruct
   TLorentzVector const rMom = kRMom + pRMom;
   float const kDca = dca(kMom.Vect(), v00, vertex);
   float const pDca = dca(pMom.Vect(), v00, vertex);
   float const kRDca = dca(kRMom.Vect(), kRPos, vertex);
   float const pRDca = dca(pRMom.Vect(), pRPos, vertex);

   TVector3 v0;
   float const dca12 = dca1To2(kRMom.Vect(), kRPos, pRMom.Vect(), pRPos, v0);
   float const decayLength = (v0 - vertex).Mag();
   float const dcaD0ToPv = dca(rMom.Vect(), v0, vertex);
   float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
   float const angle12 = kRMom.Vect().Angle(pRMom.Vect());

   TLorentzVector kRMomRest = kRMom;
   TVector3 beta;
   beta.SetMagThetaPhi(rMom.Beta(), rMom.Theta(), rMom.Phi());
   kRMomRest.Boost(-beta);
   float const cosThetaStar = rMom.Vect().Unit().Dot(kRMomRest.Vect().Unit());

   // save
   float arr[100];
   int iArr = 0;
   arr[iArr++] = centrality;
   arr[iArr++] = vertex.X();
   arr[iArr++] = vertex.Y();
   arr[iArr++] = vertex.Z();

   arr[iArr++] = kf;
   arr[iArr++] = weight;
   arr[iArr++] = b->M();
   arr[iArr++] = b->Perp();
   arr[iArr++] = b->PseudoRapidity();
   arr[iArr++] = b->Rapidity();
   arr[iArr++] = b->Phi();
   arr[iArr++] = v00.X();
   arr[iArr++] = v00.Y();
   arr[iArr++] = v00.Z();

   arr[iArr++] = rMom.M();
   arr[iArr++] = rMom.Perp();
   arr[iArr++] = rMom.PseudoRapidity();
   arr[iArr++] = rMom.Rapidity();
   arr[iArr++] = rMom.Phi();
   arr[iArr++] = v0.X();
   arr[iArr++] = v0.Y();
   arr[iArr++] = v0.Z();
   arr[iArr++] = reconstructD0(centrality, rMom);

   arr[iArr++] = dca12;
   arr[iArr++] = decayLength;
   arr[iArr++] = dcaD0ToPv;
   arr[iArr++] = cosTheta;
   arr[iArr++] = angle12;
   arr[iArr++] = cosThetaStar;

   arr[iArr++] = kMom.M();
   arr[iArr++] = kMom.Perp();
   arr[iArr++] = kMom.PseudoRapidity();
   arr[iArr++] = kMom.Rapidity();
   arr[iArr++] = kMom.Phi();
   arr[iArr++] = kDca;

   arr[iArr++] = kRMom.M();
   arr[iArr++] = kRMom.Perp();
   arr[iArr++] = kRMom.PseudoRapidity();
   arr[iArr++] = kRMom.Rapidity();
   arr[iArr++] = kRMom.Phi();
   arr[iArr++] = kRPos.X();
   arr[iArr++] = kRPos.Y();
   arr[iArr++] = kRPos.Z();
   arr[iArr++] = kRDca;

   arr[iArr++] = pMom.M();
   arr[iArr++] = pMom.Perp();
   arr[iArr++] = pMom.PseudoRapidity();
   arr[iArr++] = pMom.Rapidity();
   arr[iArr++] = pMom.Phi();
   arr[iArr++] = pDca;

   arr[iArr++] = pRMom.M();
   arr[iArr++] = pRMom.Perp();
   arr[iArr++] = pRMom.PseudoRapidity();
   arr[iArr++] = pRMom.Rapidity();
   arr[iArr++] = pRMom.Phi();
   arr[iArr++] = pRPos.X();
   arr[iArr++] = pRPos.Y();
   arr[iArr++] = pRPos.Z();
   arr[iArr++] = pRDca;

   arr[iArr++] = matchHft(1, vertex.z(), centrality, kRMom);
   arr[iArr++] = matchHft(0, vertex.z(), centrality, pRMom);

   nt->Fill(arr);
}

void getKinematics(TLorentzVector& b, double const mass)
{
   float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   float const phi = TMath::TwoPi() * gRandom->Rndm();

   float const mT = sqrt(mass * mass + pt * pt);
   float const pz = mT * sinh(y);
   float const E = mT * cosh(y);

   b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
   TVector3 posDiff = pos - vertex;
   return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
}

float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0)
{
   TVector3 posDiff = pos2 - pos1;
   TVector3 pu1 = p1.Unit();
   TVector3 pu2 = p2.Unit();
   double pu1Pu2 = pu1.Dot(pu2);
   double g = posDiff.Dot(pu1);
   double k = posDiff.Dot(pu2);
   double s2 = (k - pu1Pu2 * g) / (pu1Pu2 * pu1Pu2 - 1.);
   double s1 = g + s2 * pu1Pu2;
   TVector3 posDca1 = pos1 + pu1 * s1;
   TVector3 posDca2 = pos2 + pu2 * s2;
   v0 = 0.5 * (posDca1 + posDca2);
   return (posDca1 - posDca2).Mag();
}

TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution)
{
   float const pt = b.Perp();
   float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt));

   TLorentzVector sMom;
   sMom.SetXYZM(sPt * cos(b.Phi()), sPt * sin(b.Phi()), sPt * sinh(b.PseudoRapidity()), b.M());
   return sMom;
}

TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos)
{
   float thetaMCS = 13.6 / mom.Beta() / rMom.P() / 1000 * sqrt(pxlLayer1Thickness / fabs(sin(mom.Theta())));
   float sigmaMCS = thetaMCS * 28000 / fabs(sin(mom.Theta()));
   float sigmaPos = sqrt(pow(sigmaMCS, 2) + pow(sigmaPos0, 2));

   return TVector3(gRandom->Gaus(pos.X(), sigmaPos), gRandom->Gaus(pos.Y(), sigmaPos), gRandom->Gaus(pos.Z(), sigmaPos));
}

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

int getPhiIndex(double Phi)
{
   for (int i = 0; i < nPtBins; i++)
   {
      if ((Phi >= PhiEdge[i]) && (Phi < PhiEdge[i + 1]))
         return i;
   }
   return nPhis - 1 ;
}

TVector3 smearPosData(int const iParticleIndex, double const vz, int const cent, TLorentzVector const& rMom, TVector3 const& pos)
{
   int const iEtaIndex = getEtaIndex(rMom.PseudoRapidity());
   int const iVzIndex = getVzIndex(vz);
   int const iPtIndex = getPtIndex(rMom.Perp());

   float sigmaPosZ = 0;
   float sigmaPosXY = 0;

   if (h1DcaZ1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetEntries())
     sigmaPosZ = h1DcaZ1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom() * 1e4;

   if (h1DcaXY1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetEntries())
     sigmaPosXY = h1DcaXY1[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom() * 1e4;

   TVector3 newPos(pos);
   newPos.SetZ(0);
   TVector3 momPerp(-rMom.Vect().Y(), rMom.Vect().X(), 0.0);
   newPos += momPerp.Unit() * sigmaPosXY;

   return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

TVector3 getVertex(int const centrality)
{
   double rdmVz;
   if (h1Vz[centrality]->GetEntries() == 0) rdmVz = 0.;
   else rdmVz = h1Vz[centrality]->GetRandom() * 1e4;
   return TVector3(0., 0., rdmVz);
}

bool reconstructD0(int const centrality, TLorentzVector const& mom)
{
   TGraph* gr = NULL;

   if (centrality < 4) gr = grEff[0];
   else if (centrality < 7) gr = grEff[1];
   else gr = grEff[2];

   return gRandom->Rndm() < gr->Eval(mom.Perp());
}

bool matchHft(int const iParticleIndex, double const vz, int const cent, TLorentzVector const& mom)
{
   int const iEtaIndex = getEtaIndex(mom.PseudoRapidity());
   int const iVzIndex = getVzIndex(vz);
   // int const iPhiIndex = getPhiIndex(mom.Phi());
   int const bin = hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][cent]->FindBin(mom.Perp());
   return gRandom->Rndm() < hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][cent]->GetBinContent(bin);
}
//___________
void bookObjects()
{
   result = new TFile(outFileName.c_str(), "recreate");
   result->cd();

   TH1::AddDirectory(false);
   nt = new TNtuple("nt", "", "cent:vx:vy:vz:"
                    "pid:w:m:pt:eta:y:phi:v0x:v0y:v0z:" // MC D0
                    "rM:rPt:rEta:rY:rPhi:rV0x:rV0y:rV0z:reco:" // Rc D0
                    "dca12:decayLength:dcaD0ToPv:cosTheta:angle12:cosThetaStar:" // Rc pair
                    "kM:kPt:kEta:kY:kPhi:kDca:" // MC Kaon
                    "kRM:kRPt:kREta:kRY:kRPhi:kRVx:kRVy:kRVz:kRDca:" // Rc Kaon
                    "pM:pPt:pEta:pY:pPhi:pDca:" // MC Pion1
                    "pRM:pRPt:pREta:pRY:pRPhi:pRVx:pRVy:pRVz:pRDca:" // Rc Pion1
                    "kHft:pHft");

   TFile f("momentum_resolution.root");
   fPionMomResolution = (TF1*)f.Get("fPion")->Clone("fPionMomResolution");
   fKaonMomResolution = (TF1*)f.Get("fKaon")->Clone("fKaonMomResolution");
   f.Close();

   TFile fPP("pp200_spectra.root");
   fWeightFunction = (TF1*)fPP.Get("run12/f1Levy")->Clone("fWeightFunction");
   fPP.Close();

   TFile fHftRatio("HFT_Ratio_VsPt_Centrality.root");
   TFile fDca("Dca_VsPt_Centrality.root");
   TFile fVertex("Run14_After107_Vz_Cent.root");

   for (int ii = 0; ii < nCent; ++ii)
   {
      hHftRatio[ii] = (TH1D*)(fHftRatio.Get(Form("mh1HFTRatio1_%i", ii))->Clone(Form("mh1HFTRatio1_%i", ii)));
      h1Vz[ii]      = (TH1D*)(fVertex.Get(Form("mh1Vz_%i", ii))->Clone(Form("mh1Vz_%i", ii)));

      for (int jj = 0; jj < nPtBins; ++jj)
      {
         h1DcaXY[ii][jj] = (TH1D*)((fDca.Get(Form("mh3DcaXy_Cent%i_Pt%i", ii, jj)))->Clone(Form("mh3DcaXy_Cent%i_Pt%i", ii, jj)));
         h1DcaZ[ii][jj]  = (TH1D*)((fDca.Get(Form("mh3DcaZ_Cent%i_Pt%i", ii, jj)))->Clone(Form("mh3DcaZ_Cent%i_Pt%i", ii, jj)));
      }
   }

   fHftRatio.Close();
   fDca.Close();
   fVertex.Close();

   TFile fHftRatio1("HFT_Ratio_VsPt_Centrality_Eta_Phi_Vz_Zdcx.root");
   TFile fDca1("Dca_VsPt_Centrality_Eta_Phi_Vz_Zdcx.root");
   for (int iParticle = 0; iParticle < nParticles; ++iParticle)
   {
      for (int iEta = 0; iEta < nEtas; ++iEta)
      {
         for (int iVz = 0; iVz < nVzs; ++iVz)
         {
            for (int ii = 0; ii < nCent; ++ii)
            {
               hHftRatio1[iParticle][iEta][iVz][ii] = (TH1D*)(fHftRatio1.Get(Form("mh1HFT1PtCentPartEtaVzRatio_%i_%i_%i_%i", iParticle, iEta, iVz, ii))->Clone(Form("mh1HFT1PtCentPartEtaVzRatio_%i_%i_%i_%i", iParticle, iEta, iVz, ii)));
               for (int jj = 0; jj < nPtBins; ++jj)
               {
                  h1DcaXY1[iParticle][iEta][iVz][ii][jj] = (TH1D*)((fDca1.Get(Form("mh1DcaXyPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj)))->Clone(Form("mh1DcaXyPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj)));
                  h1DcaZ1[iParticle][iEta][iVz][ii][jj] = (TH1D*)((fDca1.Get(Form("mh1DcaZPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj)))->Clone(Form("mh1DcaZPtCentPartEtaVz_%i_%i_%i_%i_%i", iParticle, iEta, iVz, ii, jj)));
               }
            }
         }
      }
   }

   fHftRatio1.Close();
   fDca1.Close();

   grEff[0] = new TGraph("eff_4080.csv", "%lg %lg", ",");
   grEff[1] = new TGraph("eff_1040.csv", "%lg %lg", ",");
   grEff[2] = new TGraph("eff_010.csv", "%lg %lg", ",");
}
//___________
void write()
{
   result->cd();
   nt->Write();
   result->Close();
}
