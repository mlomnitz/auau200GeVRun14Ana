#include <cmath>

#include "TMath.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TProfile.h"

#include "StD0Hists.h"
#include "StPicoCharmContainers/StKaonPion.h"

StD0Hists::StD0Hists(std::string fileBaseName = "", int harmonic)
{
   hCentVzPsi = new TH3F("hCentVzPsi", "hCentVzPsi", 9, 0, 9, 10, -6, 6, 10, 0, (2.0/harmonic)*TMath::Pi());
   hCentVzPsiSameEventNoWeight = new TH3F("hCentVzPsiSameEventNoWeight", "hCentVzPsiSameEventNoWeight", 9, 0, 9, 10, -6, 6, 10, 0, (2.0/harmonic)*TMath::Pi());
   hCentVzPsiSameEvent = new TH3F("hCentVzPsiSameEvent", "hCentVzPsiSameEvent", 9, 0, 9, 10, -6, 6, 10, 0, (2.0/harmonic)*TMath::Pi());

   //D0 histograms
   const int nDimDaug = 5;
   int nBinsDaug[nDimDaug] = {9, 100, 10, 250, 10};//cent, pt, daughterpt1, m, daughterpt2
   double xMinDaug[nDimDaug] = {0, 0, 0.6, 0, 0.6};
   double xMaxDaug[nDimDaug] = {9, 10, 1.6, 2.5, 1.6};
   
   hCentVzPsi = new TH3F("hCentVzPsi", "hCentVzPsi", 9, 0, 9, 10, -6, 6, 10, 0, (2.0/harmonic)*TMath::Pi());
   hCentVzPsiSameEventNoWeight = new TH3F("hCentVzPsiSameEventNoWeight", "hCentVzPsiSameEventNoWeight", 9, 0, 9, 10, -6, 6, 10, 0, (2.0/harmonic)*TMath::Pi());
   hCentVzPsiSameEvent = new TH3F("hCentVzPsiSameEvent", "hCentVzPsiSameEvent", 9, 0, 9, 10, -6, 6, 10, 0, (2.0/harmonic)*TMath::Pi());
   //
   // EtaSub Histograms
   const int nDimEtaSub = 4;
   int nBinsEtaSub[nDimEtaSub] = {9, 100, 50, 10};//cent, pt, m, dPhi
   double xMinEtaSub[nDimEtaSub] = {0, 0, 1.6, 0};
   double xMaxEtaSub[nDimEtaSub] = {9, 10, 2.1, (2.0/harmonic)*TMath::Pi()};
   std::string etaSubTitle = "hD0EtaSubCentPtMDphi";
   // Eta Gap Histograms
   const int nDim2 = 5;
   int nBins2[nDim2] = {9, 100, 50, 10, 8};//cent, pt, m, dPhi, etaGap
   double xMin2[nDim2] = {0, 0, 1.6, 0, 0};
   double xMax2[nDim2] = {9, 10, 2.1, (2.0/harmonic)*TMath::Pi(), 0.8};
   std::string etaGapTitle = "hD0EtaSubCentPtMDphiMixed";
   
   for(int ii = 0 ; ii< topoCuts::nCutsSets ; ++ii){
     //Eta sub
     hD0EtaSubCentPtMDphi[ii] = new THnF(Form("%s_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), Form("%s_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), nDimEtaSub, nBinsEtaSub, xMinEtaSub, xMaxEtaSub);
     hD0EtaSubCentPtMDphiLikeSign[ii] = new THnF(Form("%sLikeSign_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), Form("%sLikeSign_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), nDimEtaSub, nBinsEtaSub, xMinEtaSub, xMaxEtaSub);
     //Eta gap
     hD0CentPtMDphiEtaGap[ii] = new THnF(Form("%s_%s",etaGapTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), Form("%s_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), nDim2, nBins2, xMin2, xMax2);
     hD0CentPtMDphiEtaGapLikeSign[ii] = new THnF(Form("%sLikeSign_%s",etaGapTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), Form("%sLikeSign_%s",etaSubTitle.c_str(),topoCuts::cutsSetName[ii].c_str()), nDim2, nBins2, xMin2, xMax2);
   }

}
StD0Hists::~StD0Hists()
{
}
void StD0Hists::closeFile()
{

   return;
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
