/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for studying the D0 we see in data
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
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
#include "TRandom.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"

using namespace std;

void setDecayChannels(int const mdme);
void decayAndFill(int const kf, int const decayChennel, TLorentzVector* b, double const weight, TClonesArray& daughters);
void fill(int const kf, int const decayChennel, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& pMom, TVector3 const& v00);
void getKinematics(TLorentzVector& b, double const mass);
TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos);
TVector3 smearPosData(int const cent, TLorentzVector const& rMom, TVector3 const& pos);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
bool matchHft(int const centrality, TLorentzVector const& mom);
bool reconstructD0(int const centrality, TLorentzVector const& mom);
void bookObjects();
void write();
int getptIndex(double);

TPythia6Decayer* pydecay;
TNtuple* nt;
TFile* result;

TF1* fKaonMomResolution = NULL;
TF1* fPionMomResolution = NULL;
TF1* fWeightFunction = NULL;
TGraph* grEff[3];
const Int_t nCent = 9;
const Int_t nPtBins = 35;
const Double_t ptEdge[nPtBins + 1] = { 0.0, 0.2, 0.4,  0.6,  0.8,
                                       1.0, 1.2, 1.4,  1.6,  1.8,
                                       2.0, 2.2, 2.4,  2.6,  2.8,
                                       3.0, 3.2, 3.4,  3.6,  3.8,
                                       4.0, 4.2, 4.4,  4.6,  4.8,
                                       5.0, 5.4, 5.8,  6.2,  6.6,
                                       7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
                                     };

TH1D* hHftRatio[nCent];
TH1D* h1DcaZ[nCent][nPtBins];
TH1D* h1DcaXY[nCent][nPtBins];

string outFileName = "D0Bump.toyMc.root";
std::pair<int, int> const decayChannels(673, 807);
std::pair<float, float> const momentumRange(0.3, 12);

float const acceptanceRapidity = 1.0;
float const sigmaPos0 = 15.2;
float const pxlLayer1Thickness = 0.00486;
float const sigmaVertexCent[nCent] = {31., 18.1, 12.8, 9.3, 7.2, 5.9, 5., 4.6, 4.};

//============== main  program ==================
void toyMcEffD0Bump(int npart = 100)
{
   gRandom->SetSeed();
   bookObjects();

   pydecay = TPythia6Decayer::Instance();
   pydecay->Init();

   TPythia6::Instance()->SetMDME(554, 1, 0); // π0 --> γγ
   TPythia6::Instance()->SetMDME(555, 1, 0); // π0 --> γe-e+

   TPythia6::Instance()->SetMDME(581, 1, 1); // ρ+ --> π+ π0
   TPythia6::Instance()->SetMDME(582, 1, 0); // ρ+ --> π+ gamma

   TPythia6::Instance()->SetMDME(556, 1, 1); // ρ0 --> π+ π-
   TPythia6::Instance()->SetMDME(557, 1, 0); // ρ0 --> π0 γ
   TPythia6::Instance()->SetMDME(558, 1, 0); // ρ0 --> η  γ                                                           
   TPythia6::Instance()->SetMDME(559, 1, 0); // ρ0 --> μ- μ+                                                             
   TPythia6::Instance()->SetMDME(560, 1, 0); // ρ0 --> e- e+                                                              

   TPythia6::Instance()->SetMDME(636, 1, 0); // K*- --> K0 pi-
   TPythia6::Instance()->SetMDME(637, 1, 1); // K*- --> K- π0
   TPythia6::Instance()->SetMDME(638, 1, 0); // K*- --> K- γ

   cout << "Decay channel  = " << 763 << " : D0 --> K- pi+" << endl;
   cout << "Decay channel  = " << 785 << " : D0 --> K- pi+ pi0" << endl;
   cout << "Decay channel  = " << 765 << " : D0 --> K- rho+ --> K- pi+ pi0" << endl;
   cout << "Decay channel  = " << 764 << " : D0 --> K*- pi+  --> K- pi0 pi+" << endl;
   cout << "Decay channel  = " << 786 << " : D0 --> K- pi+ rho0 --> K- pi+ pi+ pi-" <<endl;
   cout << "Decay channel  = " << 719 << " : D+ --> K- pi+ pi+" << endl;

   TLorentzVector* b_d = new TLorentzVector;
   TClonesArray ptl("TParticle", 10);
   for (int ipart = 0; ipart < npart; ipart++)
   {
      if (!(ipart % 100000))
         cout << "____________ ipart = " << ipart << " ________________" << endl;

      getKinematics(*b_d, M_D_0);

      setDecayChannels(763);
      decayAndFill(421, 763, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
      decayAndFill(-421, 763, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);

      setDecayChannels(785);
      decayAndFill(421, 785, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
      decayAndFill(-421, 785, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);

      setDecayChannels(765);
      decayAndFill(421, 765, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
      decayAndFill(-421, 765, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);

      setDecayChannels(764);
      decayAndFill(421, 764, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
      decayAndFill(-421, 764, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);

      setDecayChannels(786);
      decayAndFill(421, 786, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
      decayAndFill(-421, 786, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);

      // D+ 
      b_d->SetVectM(b_d->Vect(),M_D_PLUS);

      setDecayChannels(719);
      decayAndFill(411, 719, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
      decayAndFill(-411, 719, b_d, fWeightFunction->Eval(b_d->Perp()), ptl);
   }

   write();
}

void setDecayChannels(int const mdme)
{
   for (int idc = decayChannels.first; idc < decayChannels.second + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
   TPythia6::Instance()->SetMDME(mdme, 1, 1);
}

void decayAndFill(int const kf, int const decayChannel, TLorentzVector* b, double const weight, TClonesArray& daughters)
{
   pydecay->Decay(kf, b);
   pydecay->ImportParticles(&daughters);

   TLorentzVector kMom;
   TLorentzVector pi1Mom;
   TLorentzVector pi2Mom;
   TVector3 v00;

   int nTrk = daughters.GetEntriesFast();
   for (int iTrk = 0; iTrk < nTrk; ++iTrk)
   {
      TParticle* ptl0 = (TParticle*)daughters.At(iTrk);
      TLorentzVector tmp;

      switch (abs(ptl0->GetPdgCode()))
      {
         case 321:
            ptl0->Momentum(kMom);
            v00.SetXYZ(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.); // converted to μm
            break;
         case 211:
            if(!pi1Mom.P()) ptl0->Momentum(pi1Mom);
            else ptl0->Momentum(pi2Mom);
            break;
         default:
            break;
      }
   }
   daughters.Clear();

   fill(kf,decayChannel, b,weight,kMom,pi1Mom,v00);
   if((decayChannel==719 || decayChannel==786)) fill(kf,decayChannel,b,weight,kMom,pi2Mom,v00);
}

void fill(int const kf, int const decayChannel, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& pMom, TVector3 const& v00)
{
   // smear momentum
   TLorentzVector const kRMom = smearMom(kMom, fKaonMomResolution);
   TLorentzVector const pRMom = smearMom(pMom, fPionMomResolution);

   int const cent = floor(nCent * gRandom->Rndm());
   // smear position
   TVector3 const kRPos = smearPosData(cent, kRMom, v00);
   TVector3 const pRPos = smearPosData(cent, pRMom, v00);
   // TVector3 const kRPos = smearPos(kMom, kRMom, v00);
   // TVector3 const pRPos = smearPos(kMom, pRMom, v00);

   // smear primary vertex
   // float const sigmaVertex = sigmaVertexCent[cent];
   // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));
   TVector3 const vertex(0., 0., 0.);

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

   // misPID
   TLorentzVector kMisPidMom = kRMom;
   TLorentzVector pMisPidMom = pRMom;
   kMisPidMom.SetVectM(kRMom.Vect(), M_PION_PLUS);
   pMisPidMom.SetVectM(pRMom.Vect(), M_KAON_PLUS);
   TLorentzVector const rMisPidMom = kMisPidMom + pMisPidMom;

   // save
   float arr[100];
   int iArr = 0;
   arr[iArr++] = decayChannel;
   arr[iArr++] = cent;
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
   arr[iArr++] = rMisPidMom.M();
   arr[iArr++] = rMom.Perp();
   arr[iArr++] = rMom.PseudoRapidity();
   arr[iArr++] = rMom.Rapidity();
   arr[iArr++] = rMom.Phi();
   arr[iArr++] = reconstructD0(cent, rMom);

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

   arr[iArr++] = matchHft(cent, kRMom);
   arr[iArr++] = matchHft(cent, pRMom);

   nt->Fill(arr);
}

void getKinematics(TLorentzVector& b, double const mass)
{
   float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   float const phi = TMath::TwoPi() * gRandom->Rndm();

   float const mT = sqrt(mass*mass+pt*pt);
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
   float const sPt = gRandom->Gaus(pt, 1.2 * pt * fMomResolution->Eval(pt));

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

int getptIndex(double pT)
{
   for (int i = 0; i < nPtBins; i++)
   {
      if ((pT >= ptEdge[i]) && (pT < ptEdge[i + 1]))
         return i;
   }
}

TVector3 smearPosData(int const cent, TLorentzVector const& rMom, TVector3 const& pos)
{
   int ptIndex = getptIndex(rMom.Perp());
   float sigmaPosZ = h1DcaZ[cent][ptIndex]->GetRandom() * 1e4;
   float sigmaPosXY = h1DcaXY[cent][ptIndex]->GetRandom() * 1e4;

   TVector3 newPos(pos);
   newPos.SetZ(0);
   TVector3 momPerp(-rMom.Vect().Y(), rMom.Vect().X(), 0.0);
   newPos += momPerp.Unit() * sigmaPosXY;

   return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

bool reconstructD0(int const centrality, TLorentzVector const& mom)
{
   TGraph* gr = NULL;

   if (centrality < 4) gr = grEff[0];
   else if (centrality < 7) gr = grEff[1];
   else gr = grEff[2];

   return gRandom->Rndm() < gr->Eval(mom.Perp());
}

bool matchHft(int const cent, TLorentzVector const& mom)
{
   int const bin = hHftRatio[cent]->FindBin(mom.Perp());
   return gRandom->Rndm() < hHftRatio[cent]->GetBinContent(bin);
}

//___________
void bookObjects()
{
   result = new TFile(outFileName.c_str(), "recreate");
   result->cd();

   TH1::AddDirectory(false);
   nt = new TNtuple("nt", "", "decayChannel:cent:vx:vy:vz:"
                    "pid:w:m:pt:eta:y:phi:v0x:v0y:v0z:" // MC D0
                    "rM:misPidM:rPt:rEta:rY:rPhi:reco:" // Rc D0
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

   for (int ii = 0; ii < 9; ++ii)
   {
      hHftRatio[ii] = (TH1D*)(fHftRatio.Get(Form("mh1HFTRatio1_%i", ii))->Clone(Form("mh1HFTRatio1_%i", ii)));
      result->cd();
      hHftRatio[ii]->Write();
   }

   result->cd();
   for (int ii = 0; ii < nCent; ++ii)
   {
      for (int jj = 0; jj < nPtBins; ++jj)
      {
         h1DcaXY[ii][jj] = (TH1D*)((fDca.Get(Form("mh3DcaXy_Cent%i_Pt%i", ii, jj)))->Clone(Form("mh3DcaXy_Cent%i_Pt%i", ii, jj)));
         h1DcaZ[ii][jj]  = (TH1D*)((fDca.Get(Form("mh3DcaZ_Cent%i_Pt%i", ii, jj)))->Clone(Form("mh3DcaZ_Cent%i_Pt%i", ii, jj)));
      }
   }

   fHftRatio.Close();
   fDca.Close();

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
