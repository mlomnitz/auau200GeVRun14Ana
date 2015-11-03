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
#include "dataDrivenFastSimulator.h"

using namespace std;

void setDecayChannels(int const mdme);
void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TVector3 v00);
void getKinematics(TLorentzVector& b, double const mass);
bool reconstructD0(int const centrality, TLorentzVector const& mom);
void bookObjects();
void write();

TPythia6Decayer* pydecay;
TNtuple* nt;
TFile* result;

TGraph* grEff[3];
TF1* fWeightFunction = NULL;

string outFileName = "D0.toyMc.root";
std::pair<int, int> const decayChannels(747, 807);
std::pair<float, float> const momentumRange(0, 12);

float const acceptanceRapidity = 1.0;

//============== main  program ==================
void toyMcEffD0(int npart = 100)
{
   loadAllDistributions();

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
   TLorentzVector const kRMom = smearMom(1,kMom);
   TLorentzVector const pRMom = smearMom(0,pMom);

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

   int const charge = kf > 0? 1: -1;

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
   arr[iArr++] = tpcReconstructed(1,-1*charge,centrality,kRMom);

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
   arr[iArr++] = tpcReconstructed(0,charge,centrality,pRMom);

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

bool reconstructD0(int const centrality, TLorentzVector const& mom)
{
   TGraph* gr = NULL;

   if (centrality < 4) gr = grEff[0];
   else if (centrality < 7) gr = grEff[1];
   else gr = grEff[2];

   return gRandom->Rndm() < gr->Eval(mom.Perp());
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
                    "kRM:kRPt:kREta:kRY:kRPhi:kRVx:kRVy:kRVz:kRDca:kTpc:" // Rc Kaon
                    "pM:pPt:pEta:pY:pPhi:pDca:" // MC Pion1
                    "pRM:pRPt:pREta:pRY:pRPhi:pRVx:pRVy:pRVz:pRDca:pTpc:" // Rc Pion1
                    "kHft:pHft");

   cout << "Loading input spectra ..." << endl;
   TFile fPP("pp200_spectra.root");
   fWeightFunction = (TF1*)fPP.Get("run12/f1Levy")->Clone("f1Levy");
   fPP.Close();

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
