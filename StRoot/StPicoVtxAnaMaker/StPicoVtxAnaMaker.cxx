#include <vector>
#include <cmath>
#include <algorithm>
// #include <random>

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "../StPicoDstMaker/StPicoDst.h"
#include "../StPicoDstMaker/StPicoDstMaker.h"
#include "../StPicoDstMaker/StPicoEvent.h"
#include "../StPicoDstMaker/StPicoTrack.h"
#include "StCuts.h"

#include "StPicoVtxAnaMaker.h"

ClassImp(StPicoVtxAnaMaker)

StPicoVtxAnaMaker::StPicoVtxAnaMaker(char const* makerName, StPicoDstMaker* picoMaker, char const* fileBaseName)
   : StMaker(makerName), mPicoDstMaker(picoMaker), mPicoEvent(NULL), mKfVertexFitter(), mVtxEvent(fileBaseName)
{
}

StPicoVtxAnaMaker::~StPicoVtxAnaMaker()
{
}

Int_t StPicoVtxAnaMaker::Init()
{
   return kStOK;
}

Int_t StPicoVtxAnaMaker::Finish()
{
   mVtxEvent.closeFile();
   return kStOK;
}

void StPicoVtxAnaMaker::Clear(Option_t *opt)
{
}

Int_t StPicoVtxAnaMaker::Make()
{
   if (!mPicoDstMaker)
   {
      LOG_WARN << " No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const * picoDst = mPicoDstMaker->picoDst();
   if (!picoDst)
   {
      LOG_WARN << " No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   mPicoEvent = picoDst->event();

   StThreeVectorF kfVtx(-999.,-999.,-999.);
   StThreeVectorF kfHftVtx(-999.,-999.,-999.);
   StThreeVectorF kfTopVtx(-999.,-999.,-999.);
   StThreeVectorF kfBottomVtx(-999.,-999.,-999.);
   StThreeVectorF kfRightVtx(-999.,-999.,-999.);
   StThreeVectorF kfLeftVtx(-999.,-999.,-999.);
   StThreeVectorF kfVtxSubEvt1(-999.,-999.,-999.);
   StThreeVectorF kfVtxSubEvt2(-999.,-999.,-999.);
   StThreeVectorF kfHftVtxSubEvt1(-999.,-999.,-999.);
   StThreeVectorF kfHftVtxSubEvt2(-999.,-999.,-999.);
   StThreeVectorF kfHftPxlSecVtx[10];
   StThreeVectorF kfHftHighPtPxlSecVtx[10];

   int nTrks = 0;
   int nTrksHft = 0;
   int nTrksTop = 0;
   int nTrksBottom = 0;
   int nTrksRight = 0;
   int nTrksLeft = 0;
   int nTrksSubEvt1 = 0;
   int nTrksSubEvt2 = 0;
   int nTrksHftSubEvt1 = 0;
   int nTrksHftSubEvt2 = 0;
   int nTrksHftPxlSec[10];
   int nTrksHftHighPtPxlSec[10];
   
   for(int i=0; i<10; ++i)
   { 
     kfHftPxlSecVtx[i].set(-999.,-999.,-999.);
     kfHftHighPtPxlSecVtx[i].set(-999.,-999.,-999.);
     nTrksHftPxlSec[i] = 0;
     nTrksHftHighPtPxlSec[i] = 0;
   }

   if (isGoodEvent())
   {
      UInt_t nTracks = picoDst->numberOfTracks();

      std::vector<int> allTracksForVtxFit;
      std::vector<int> hftTracksForVtxFit;
      std::vector<int> topTracksForVtxFit;
      std::vector<int> bottomTracksForVtxFit;
      std::vector<int> rightTracksForVtxFit;
      std::vector<int> leftTracksForVtxFit;
      std::vector<int> tracksPerPxlSecForVtxFit[10];
      std::vector<int> highPtTracksPerPxlSecForVtxFit[10];

      StThreeVectorF const pVtx = mPicoEvent->primaryVertex();

      for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
      {
         StPicoTrack* trk = picoDst->track(iTrack);
         if(!trk) continue;

         if(!isGoodForVertexFit(trk,pVtx)) continue;

         allTracksForVtxFit.push_back(iTrack);
         if(trk->isHFTTrack()) hftTracksForVtxFit.push_back(iTrack);

         float phi = trk->gMom(pVtx,mPicoEvent->bField()).phi();

         if(phi > 0 && phi < TMath::Pi()) topTracksForVtxFit.push_back(iTrack);
         else bottomTracksForVtxFit.push_back(iTrack);

         if(phi > -TMath::Pi()/2. && phi < TMath::Pi()/2.) rightTracksForVtxFit.push_back(iTrack);
         else leftTracksForVtxFit.push_back(iTrack);

         if(!trk->isHFTTrack()) continue;
         if(int nSec = cuts::getPxlSectorNumber(phi))
         {
           tracksPerPxlSecForVtxFit[nSec-1].push_back(iTrack);
           if(trk->gPt() > cuts::pxlVtxFitPtCut) highPtTracksPerPxlSecForVtxFit[nSec-1].push_back(iTrack);
         }
      } // .. end tracks loop

      if(!allTracksForVtxFit.empty())
      {
        // std::random_device rd;
        // std::mt19937 g(rd());
        std::random_shuffle(allTracksForVtxFit.begin(),allTracksForVtxFit.end());

        int middleElement = static_cast<int>(allTracksForVtxFit.size()/2);
        std::vector<int> tracksForVtxFitSub1(middleElement);
        std::vector<int> tracksForVtxFitSub2(allTracksForVtxFit.size()-middleElement);

        std::copy(allTracksForVtxFit.begin(),allTracksForVtxFit.begin()+middleElement,tracksForVtxFitSub1.begin());
        std::copy(allTracksForVtxFit.begin()+middleElement,allTracksForVtxFit.end(),tracksForVtxFitSub2.begin());

        nTrks = allTracksForVtxFit.size();
        nTrksSubEvt1 = tracksForVtxFitSub1.size();
        nTrksSubEvt2 = tracksForVtxFitSub2.size();

        kfVtx = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,allTracksForVtxFit);
        kfVtxSubEvt1 = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,tracksForVtxFitSub1);
        kfVtxSubEvt2 = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,tracksForVtxFitSub2);
      }

      if(!hftTracksForVtxFit.empty())
      {
        // std::random_device rd;
        // std::mt19937 g(rd());
        // std::random_shuffle(hftTracksForVtxFit.begin(),hftTracksForVtxFit.end(),g);
        std::random_shuffle(hftTracksForVtxFit.begin(),hftTracksForVtxFit.end());

        int middleElement = static_cast<int>(hftTracksForVtxFit.size()/2);
        std::vector<int> hftTracksForVtxFitSub1(middleElement);
        std::vector<int> hftTracksForVtxFitSub2(hftTracksForVtxFit.size()-middleElement);

        std::copy(hftTracksForVtxFit.begin(),hftTracksForVtxFit.begin()+middleElement,hftTracksForVtxFitSub1.begin());
        std::copy(hftTracksForVtxFit.begin()+middleElement,hftTracksForVtxFit.end(),hftTracksForVtxFitSub2.begin());

        nTrksHft = hftTracksForVtxFit.size();
        nTrksHftSubEvt1 = hftTracksForVtxFitSub1.size();
        nTrksHftSubEvt2 = hftTracksForVtxFitSub2.size();

        kfHftVtx = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,hftTracksForVtxFit);
        kfHftVtxSubEvt1 = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,hftTracksForVtxFitSub1);
        kfHftVtxSubEvt2 = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,hftTracksForVtxFitSub2);
      }

      if(!topTracksForVtxFit.empty())
      {
        nTrksTop = topTracksForVtxFit.size();
        kfTopVtx = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,topTracksForVtxFit);
      }

      if(!bottomTracksForVtxFit.empty())
      {
        nTrksBottom = bottomTracksForVtxFit.size();
        kfBottomVtx = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,bottomTracksForVtxFit);
      }

      if(!rightTracksForVtxFit.empty())
      {
        nTrksRight = rightTracksForVtxFit.size();
        kfRightVtx = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,rightTracksForVtxFit);
      }

      if(!leftTracksForVtxFit.empty())
      {
        nTrksLeft = leftTracksForVtxFit.size();
        kfLeftVtx = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,leftTracksForVtxFit);
      }

      for(int iSec = 0; iSec < 10; ++iSec)
      {
        if(!tracksPerPxlSecForVtxFit[iSec].empty())
        {
          nTrksHftPxlSec[iSec] = tracksPerPxlSecForVtxFit[iSec].size();
          kfHftPxlSecVtx[iSec] = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,tracksPerPxlSecForVtxFit[iSec]);
        }

        if(!highPtTracksPerPxlSecForVtxFit[iSec].empty())
        {
          nTrksHftHighPtPxlSec[iSec] = highPtTracksPerPxlSecForVtxFit[iSec].size();
          kfHftHighPtPxlSecVtx[iSec] = mKfVertexFitter.primaryVertexRefitUsingTracks(picoDst,highPtTracksPerPxlSecForVtxFit[iSec]);
        }
      }
   } //.. end of good event fill

   mVtxEvent.addEvent(*mPicoEvent,
                      kfVtx,
                      kfHftVtx,
                      kfTopVtx,
                      kfBottomVtx,
                      kfRightVtx,
                      kfLeftVtx,
                      kfVtxSubEvt1,
                      kfVtxSubEvt2,
                      kfHftVtxSubEvt1,
                      kfHftVtxSubEvt2,
                      kfHftPxlSecVtx,
                      kfHftHighPtPxlSecVtx,
                      nTrks,
                      nTrksHft,
                      nTrksTop,
                      nTrksBottom,
                      nTrksRight,
                      nTrksLeft,
                      nTrksSubEvt1,
                      nTrksSubEvt2,
                      nTrksHftSubEvt1,
                      nTrksHftSubEvt2,
                      nTrksHftPxlSec,
                      nTrksHftHighPtPxlSec);
   return kStOK;
}

bool StPicoVtxAnaMaker::isGoodEvent()
{
   return (mPicoEvent->triggerWord() & cuts::triggerWord) &&
          fabs(mPicoEvent->primaryVertex().z()) < cuts::vz &&
          fabs(mPicoEvent->primaryVertex().z() - mPicoEvent->vzVpd()) < cuts::vzVpdVz;
}

bool StPicoVtxAnaMaker::isGoodForVertexFit(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
  StPhysicalHelixD helix = trk->dcaGeometry().helix();
  float dca = (helix.at(helix.pathLength(vtx)) - vtx).mag();

  if (dca > cuts::vtxDca) return false;

  size_t numberOfFitPoints = popcount(trk->map0() >> 1); // drop first bit, primary vertex point
  numberOfFitPoints += popcount(trk->map1() & 0x1FFFFF); // only the first 21 bits are important, see StTrackTopologyMap.cxx

  return numberOfFitPoints >= cuts::vtxNumberOfFitPoints;
}
