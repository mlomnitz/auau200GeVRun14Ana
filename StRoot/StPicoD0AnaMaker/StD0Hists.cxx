#include <cmath>

#include "TMath.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TProfile.h"
#include "TVector2.h"

#include "StD0Hists.h"
#include "StPicoCharmContainers/StKaonPion.h"
#include "StEventPlane/StEventPlane.h"

StD0Hists::StD0Hists(std::string fileBaseName = "", int harmonic)
{   
  mOutfile = new TFile(Form("%s.d0v%i.root",fileBaseName.c_str(), harmonic),"RECREATE");
  mOutfile->cd();
   hCentVzPsi = new TH3F("hCentVzPsi", "hCentVzPsi", 9, 0, 9, 10, -6, 6, 10, 0, (2.0/harmonic)*TMath::Pi());
   hCentVzPsiNoWeight = new TH3F("hCentVzPsiNoWeight", "hCentVzPsiNoWeight", 9, 0, 9, 10, -6, 6, 10, 0, (2.0/harmonic)*TMath::Pi());
   //
   // EtaSub Histograms
   const int nDimEtaSub = 4;
   int nBinsEtaSub[nDimEtaSub] = {9, 100, 50, 10};//cent, pt, m, dPhi
   double xMinEtaSub[nDimEtaSub] = {0, 0, 1.6, 0};
   double xMaxEtaSub[nDimEtaSub] = {9, 10, 2.1, (2.0/harmonic)*TMath::Pi()};
   std::string etaSubTitle = "hEtaSubD0CentPtMDphi";
   // Eta Gap Histograms
   const int nDim2 = 5;
   int nBins2[nDim2] = {9, 100, 50, 10, 8};//cent, pt, m, dPhi, etaGap
   double xMin2[nDim2] = {0, 0, 1.6, 0, 0};
   double xMax2[nDim2] = {9, 10, 2.1, (2.0/harmonic)*TMath::Pi(), 0.8};
   std::string etaGapTitle = "hD0CentPtMDphiEtaGap";
   
   for(int ii = 0 ; ii< topoCuts::nCutsSets ; ++ii){
     //Eta sub
     hD0EtaSubCentPtMDphi[ii] = new THnF(Form("%s_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), Form("%s_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), nDimEtaSub, nBinsEtaSub, xMinEtaSub, xMaxEtaSub);
     hD0EtaSubCentPtMDphiLikeSign[ii] = new THnF(Form("%sLikeSign_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), Form("%sLikeSign_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), nDimEtaSub, nBinsEtaSub, xMinEtaSub, xMaxEtaSub);
     //Eta gap
     hD0CentPtMDphiEtaGap[ii] = new THnF(Form("%s_%s",etaGapTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), Form("%s_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), nDim2, nBins2, xMin2, xMax2);
     hD0CentPtMDphiEtaGapLikeSign[ii] = new THnF(Form("%sLikeSign_%s",etaGapTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), Form("%sLikeSign_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), nDim2, nBins2, xMin2, xMax2);
   }
   mHarmonic = harmonic;
}
StD0Hists::~StD0Hists()
{
}
void StD0Hists::closeFile()
{
  mOutfile->cd();
  hCentVzPsi->Write();
  hCentVzPsiNoWeight->Write();
  for( int ii = 0 ; ii < topoCuts::nCutsSets; ++ii){
    hD0EtaSubCentPtMDphi[ii]->Write();
    hD0EtaSubCentPtMDphiLikeSign[ii]->Write();
    //
    hD0CentPtMDphiEtaGap[ii]->Write();
    hD0CentPtMDphiEtaGapLikeSign[ii]->Write();
  }
  mOutfile->Close();
  mOutfile->Delete();
}
//----------------------------------------------------------------------------- 
bool StD0Hists::isGoodPair(StKaonPion const* pair, topoCuts::TopologicalCuts const& cuts) const
{
  int ptIndex = getD0PtIndex(pair,cuts.PtEdge);
   return (pair->m() > topoCuts::massMin && pair->m() < topoCuts::massMax &&
           std::abs(pair->lorentzVector().rapidity()) < cuts.RapidityCut &&
           pair->pionDca() > cuts.pDca[ptIndex] && pair->kaonDca() > cuts.kDca[ptIndex] &&
           pair->dcaDaughters() < cuts.dcaDaughters[ptIndex] &&
           pair->decayLength() > cuts.decayLength[ptIndex] &&
           std::cos(pair->pointingAngle()) > cuts.cosTheta[ptIndex] &&
           ((pair->decayLength()) * sin(pair->pointingAngle())) < cuts.dcaV0ToPv[ptIndex]);
}
//-----------------------------------------------------------------------------
int StD0Hists::getD0PtIndex(StKaonPion const* pair, std::vector<float> const& edges) const
{
  for (int i = 0; i < edges.size(); i++)
   {
      if ((pair->pt() >= edges[i]) && (pair->pt() < edges[i + 1]))
         return i;
   }
  return edges.size() - 1;
}
//-----------------------------------------------------------------------------
void StD0Hists::addKaonPion(bool unlike, StKaonPion const* pair, float const pionEta, float const kaonEta, StEventPlane* mEventPlane, int const centrality, const double reweight)
{
  bool pass_cut_set[topoCuts::nCutsSets] = { isGoodPair(pair,topoCuts::D0Cuts),
					     isGoodPair(pair,topoCuts::D0Cuts_50eff),
					     isGoodPair(pair,topoCuts::D0Cuts_150eff)};
  bool drop_pair = false;
  for( int iCuts = 0; iCuts<topoCuts::nCutsSets; ++iCuts){
    if(pass_cut_set[iCuts])
      drop_pair = true;
  }
  if( drop_pair == false)
    return;
  //   ------------  Event sub
  TVector2 QEtaSub;
  if( pair->eta() > 0)
    QEtaSub = mEventPlane->QEtaMinusGap005();
  else
    QEtaSub = mEventPlane->QEtaPlusGap005();
  //   ------------  Remove daughters
  if (pionEta*pair->eta() <0 && std::fabs(pionEta) > 0.05 )
    QEtaSub -= mEventPlane->q(pair->pionIdx());
  if (kaonEta*pair->eta() <0 && std::fabs(kaonEta) > 0.05 )
    QEtaSub -= mEventPlane->q(pair->kaonIdx());
  
  double dPhiEtaSub = pair->phi() - QEtaSub.Phi() / mHarmonic;
  while (dPhiEtaSub < 0) dPhiEtaSub += (2.0/mHarmonic)*TMath::Pi();
  while (dPhiEtaSub >= (2.0/mHarmonic)*TMath::Pi()) dPhiEtaSub -= (2.0/mHarmonic)*TMath::Pi();
  double toFillEtaSub[5] = { 1.0*centrality+0.5, pair->pt(), pair->m(), dPhiEtaSub };
  for( int ii = 0; ii<topoCuts::nCutsSets; ++ii){
    if( !pass_cut_set[ii] ) continue;
    if( unlike )
      hD0EtaSubCentPtMDphi[ii]->Fill(toFillEtaSub);
    else
      hD0EtaSubCentPtMDphiLikeSign[ii]->Fill(toFillEtaSub);
  }
  //   ------------  Eta gap
  TVector2 mQEta[20];
  for(int iEta = 0; iEta < 20; ++iEta)
    mQEta[iEta] = mEventPlane->QEta(iEta);
  
  int iEta = (int)(pair->eta() * 10 + 10);
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
      QEtaGap -= mEventPlane->q(pair->pionIdx());
    int iEtaKaon = (int)(kaonEta * 10 + 10);
    if (fabs(iEtaKaon - iEta_) >= nEtaGaps)
      QEtaGap -= mEventPlane->q(pair->kaonIdx());
    //
    
    if (QEtaGap.Mod() == 0){
      cout << "QEtaGap.Mod()==0  nEtaGaps: " << nEtaGaps << endl;
      continue;
    }
    float dPhiEtaGap = pair->phi() - QEtaGap.Phi() / mHarmonic;
    while (dPhiEtaGap < 0) dPhiEtaGap += (2.0/mHarmonic)*TMath::Pi();
    while (dPhiEtaGap >= (2.0/mHarmonic)*TMath::Pi()) dPhiEtaGap -= (2.0/mHarmonic)*TMath::Pi();
    double toFillEtaGap[5] = {1.*centrality + 0.5, pair->pt(), pair->m(), dPhiEtaGap, 0.1 * nEtaGaps + 0.05};
    for(int ii = 0; ii < topoCuts::nCutsSets; ++ii){
      if( !pass_cut_set[ii] ) continue;
      if( unlike )
	hD0CentPtMDphiEtaGap[ii]->Fill(toFillEtaGap);
      else
	hD0CentPtMDphiEtaGapLikeSign[ii]->Fill(toFillEtaGap);
    }
  }
}
void StD0Hists::addEvent(int const centrality, float const Vz, float const Psi, float const weight)
{
  hCentVzPsiNoWeight->Fill(centrality, Vz, Psi);
  hCentVzPsi->Fill(centrality, Vz, Psi, weight);
}
