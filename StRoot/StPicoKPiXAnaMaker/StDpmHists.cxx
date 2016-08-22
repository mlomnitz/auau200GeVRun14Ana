#include <cmath>

#include "TMath.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TProfile.h"
#include "TVector2.h"

#include "StDpmHists.h"
#include "StPicoCharmContainers/StPicoKPiX.h"
#include "StEventPlane/StEventPlane.h"

StDpmHists::StDpmHists(std::string fileBaseName = "", int harmonic)
{   
  mOutfile = new TFile(Form("%s.dpmv%i.root",fileBaseName.c_str(), harmonic),"RECREATE");
  mOutfile->cd();
   hCentVzPsi = new TH3F("hCentVzPsi", "hCentVzPsi", 9, 0, 9, 10, -6, 6, 10, 0, (2.0/harmonic)*TMath::Pi());
   hCentVzPsiNoWeight = new TH3F("hCentVzPsiNoWeight", "hCentVzPsiNoWeight", 9, 0, 9, 10, -6, 6, 10, 0, (2.0/harmonic)*TMath::Pi());
   //
   // EtaSub Histograms
   const int nDimEtaSub = 4;
   int nBinsEtaSub[nDimEtaSub] = {9, 100, 50, 10};//cent, pt, m, dPhi
   double xMinEtaSub[nDimEtaSub] = {0, 0, 1.6, 0};
   double xMaxEtaSub[nDimEtaSub] = {9, 10, 2.1, (2.0/harmonic)*TMath::Pi()};
   std::string etaSubTitle = "hEtaSubDpmCentPtMDphi";
   // Eta Gap Histograms
   const int nDim2 = 5;
   int nBins2[nDim2] = {9, 100, 50, 10, 8};//cent, pt, m, dPhi, etaGap
   double xMin2[nDim2] = {0, 0, 1.6, 0, 0};
   double xMax2[nDim2] = {9, 10, 2.1, (2.0/harmonic)*TMath::Pi(), 0.8};
   std::string etaGapTitle = "hDpmCentPtMDphiEtaGap";
   
   //Eta sub
   hDpmEtaSubCentPtMDphi = new THnF(Form("%s",etaSubTitle.c_str()), Form("%s",etaSubTitle.c_str()), nDimEtaSub, nBinsEtaSub, xMinEtaSub, xMaxEtaSub);
   hDpmEtaSubCentPtMDphiLikeSign = new THnF(Form("%sLikeSign",etaSubTitle.c_str()), Form("%sLikeSign",etaSubTitle.c_str()), nDimEtaSub, nBinsEtaSub, xMinEtaSub, xMaxEtaSub);
   //Eta gap
   hDpmCentPtMDphiEtaGap = new THnF(Form("%s",etaGapTitle.c_str()), Form("%s",etaSubTitle.c_str()), nDim2, nBins2, xMin2, xMax2);
   hDpmCentPtMDphiEtaGapLikeSign = new THnF(Form("%sLikeSign",etaGapTitle.c_str()), Form("%sLikeSign",etaSubTitle.c_str()), nDim2, nBins2, xMin2, xMax2);
   mHarmonic = harmonic;
}
StDpmHists::~StDpmHists()
{
}
void StDpmHists::closeFile()
{
  mOutfile->cd();
  hCentVzPsi->Write();
  hCentVzPsiNoWeight->Write();
  hDpmEtaSubCentPtMDphi->Write();
  hDpmEtaSubCentPtMDphiLikeSign->Write();
  //
  hDpmCentPtMDphiEtaGap->Write();
  hDpmCentPtMDphiEtaGapLikeSign->Write();
  
  mOutfile->Close();
  mOutfile->Delete();
}
//-----------------------------------------------------------------------------
int StDpmHists::getDpmPtIndex(StPicoKPiX const* triplet, std::vector<float> const& edges) const
{
  StThreeVectorF threeMom = triplet->threeMom();
  for (int i = 0; i < edges.size(); i++)
   {
     if ((threeMom.perp() >= edges[i]) && (threeMom.perp() < edges[i + 1]))
       return i;
   }
  return edges.size() - 1;
}
//-----------------------------------------------------------------------------
void StDpmHists::addKaPiX(bool unlike, double const mass_hypo, StPicoKPiX const* triplet, float const pionEta, float const kaonEta, float const xaonEta, StEventPlane* mEventPlane, int const centrality, const double reweight)
{
  StLorentzVectorF fourMom = triplet->fourMom(mass_hypo);
  //   ------------  Event sub
  TVector2 QEtaSub;
  if( fourMom.pseudoRapidity() > 0)
    QEtaSub = mEventPlane->QEtaMinusGap005();
  else
    QEtaSub = mEventPlane->QEtaPlusGap005();
  //   ------------  Remove daughters
  if (pionEta*fourMom.pseudoRapidity() <0 && std::fabs(pionEta) > 0.05 )
    QEtaSub -= mEventPlane->q(triplet->pionIdx());
  if (kaonEta*fourMom.pseudoRapidity() <0 && std::fabs(kaonEta) > 0.05 )
    QEtaSub -= mEventPlane->q(triplet->kaonIdx());
  if (xaonEta*fourMom.pseudoRapidity() <0 && std::fabs(xaonEta) > 0.05 )
    QEtaSub -= mEventPlane->q(triplet->xaonIdx());
  
  double dPhiEtaSub = fourMom.phi() - QEtaSub.Phi() / mHarmonic;
  while (dPhiEtaSub < 0) dPhiEtaSub += (2.0/mHarmonic)*TMath::Pi();
  while (dPhiEtaSub >= (2.0/mHarmonic)*TMath::Pi()) dPhiEtaSub -= (2.0/mHarmonic)*TMath::Pi();
  double toFillEtaSub[5] = { 1.0*centrality+0.5, fourMom.perp(), fourMom.m(), dPhiEtaSub };
 
  if( unlike )
    hDpmEtaSubCentPtMDphi->Fill(toFillEtaSub);
  else
    hDpmEtaSubCentPtMDphiLikeSign->Fill(toFillEtaSub);
  //   ------------  Eta gap
  TVector2 mQEta[20];
  for(int iEta = 0; iEta < 20; ++iEta)
    mQEta[iEta] = mEventPlane->QEta(iEta);
  
  int iEta = (int)(fourMom.pseudoRapidity() * 10 + 10);
  for (int nEtaGaps = 0; nEtaGaps < 8; ++nEtaGaps){
    //Corresponding q-vector
    TVector2 QEtaGap(0,0);
    int iEta_ = iEta;
    if (iEta_ < nEtaGaps) iEta_ = nEtaGaps - 1;
    if (iEta_ > 20 - nEtaGaps) iEta_ = 20 - nEtaGaps;
    for( int ii = 0; ii<20; ++ii){
      if(fabs(ii-iEta) >= nEtaGaps)
	QEtaGap+=mQEta[ii];
    }
    //
    int iEtaPion = (int)(pionEta * 10 + 10);
    if (fabs(iEtaPion - iEta_) >= nEtaGaps)
      QEtaGap -= mEventPlane->q(triplet->pionIdx());
    int iEtaKaon = (int)(kaonEta * 10 + 10);
    if (fabs(iEtaKaon - iEta_) >= nEtaGaps)
      QEtaGap -= mEventPlane->q(triplet->kaonIdx());
    int iEtaXaon = (int)(xaonEta * 10 + 10);
    if (fabs(iEtaXaon - iEta_) >= nEtaGaps)
      QEtaGap -= mEventPlane->q(triplet->xaonIdx());
    //
    
    if (QEtaGap.Mod() == 0){
      cout << "QEtaGap.Mod()==0  nEtaGaps: " << nEtaGaps << endl;
      continue;
    }
    float dPhiEtaGap = fourMom.phi() - QEtaGap.Phi() / mHarmonic;
    while (dPhiEtaGap < 0) dPhiEtaGap += (2.0/mHarmonic)*TMath::Pi();
    while (dPhiEtaGap >= (2.0/mHarmonic)*TMath::Pi()) dPhiEtaGap -= (2.0/mHarmonic)*TMath::Pi();
    double toFillEtaGap[5] = {1.*centrality + 0.5, fourMom.perp(), fourMom.m(), dPhiEtaGap, 0.1 * nEtaGaps + 0.05};

    if( unlike )
      hDpmCentPtMDphiEtaGap->Fill(toFillEtaGap);
    else
      hDpmCentPtMDphiEtaGapLikeSign->Fill(toFillEtaGap);
  }
}
void StDpmHists::addEvent(int const centrality, float const Vz, float const Psi, float const weight)
{
  hCentVzPsiNoWeight->Fill(centrality, Vz, Psi);
  hCentVzPsi->Fill(centrality, Vz, Psi, weight);
}
