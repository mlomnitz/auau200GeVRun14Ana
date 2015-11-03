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
#include "dataDrivenFastSimulator.h"

using namespace std;

void setDecayChannels(int const mdme);
void decayAndFill(int const kf, int const decayChennel, TLorentzVector* b, double const weight, TClonesArray& daughters);
void fill(int const kf, int const decayChennel, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& pMom, TVector3 const& v00);
void getKinematics(TLorentzVector& b, double const mass);

void bookObjects();
void write();
TPythia6Decayer* pydecay;
TNtuple* nt;
TFile* result;

TF1* fWeightFunction = NULL;

string outFileName = "D0Bump.toyMc.root";
std::pair<int, int> const decayChannels(673, 807);
std::pair<float, float> const momentumRange(0.3, 12);

float const acceptanceRapidity = 1.0;

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
//___________
void bookObjects()
{
   result = new TFile(outFileName.c_str(), "recreate");
   result->cd();

   int BufSize = (int)pow(2., 16.);
   nt = new TNtuple("nt", "", "decayChannel:cent:vx:vy:vz:"
                    "pid:w:m:pt:eta:y:phi:v0x:v0y:v0z:" // MC D0
                    "rM:misPidM:rPt:rEta:rY:rPhi:reco:" // Rc D0
                    "dca12:decayLength:dcaD0ToPv:cosTheta:angle12:cosThetaStar:" // Rc pair
                    "kM:kPt:kEta:kY:kPhi:kDca:" // MC Kaon
                    "kRM:kRPt:kREta:kRY:kRPhi:kRVx:kRVy:kRVz:kRDca:" // Rc Kaon
                    "pM:pPt:pEta:pY:pPhi:pDca:" // MC Pion1
                    "pRM:pRPt:pREta:pRY:pRPhi:pRVx:pRVy:pRVz:pRDca:" // Rc Pion1
                    "kHft:pHft",BufSize);

   TFile fPP("pp200_spectra.root");
   fWeightFunction = (TF1*)fPP.Get("run12/f1Levy")->Clone("fWeightFunction");
   fPP.Close();
}
//___________
void write()
{
   result->cd();
   nt->Write();
   result->Close();
}
