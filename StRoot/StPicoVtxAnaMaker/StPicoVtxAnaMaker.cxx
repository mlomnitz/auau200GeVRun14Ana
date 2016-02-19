#include <vector>
#include <cmath>
#include <algorithm>

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "../StPicoDstMaker/StPicoDst.h"
#include "../StPicoDstMaker/StPicoDstMaker.h"
#include "../StPicoDstMaker/StPicoEvent.h"
#include "../StPicoDstMaker/StPicoTrack.h"
#include "StCuts.h"

#include "StPicoVtxAnaMaker.h"

ClassImp(StPicoVtxAnaMaker)

StPicoVtxAnaMaker::StPicoVtxAnaMaker(char const* makerName, StPicoDstMaker* picoMaker, char const* fileBaseName)
   : StMaker(makerName), mPicoDstMaker(picoMaker), mPicoEvent(NULL), mKfVertexFitter(), mVtxEvent(fileBaseName),
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
   mKfVertexEvent.closeFile();
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

   StThreeVectorF kfVertex(-999.,-999.,-999.);
   StThreeVectorF kfVertexSubEvt1(-999.,-999.,-999.);
   StThreeVectorF kfVertexSubEvt2(-999.,-999.,-999.);
   int nTracksFullEvt = 0;
   int nTracksSubEvt1 = 0;
   int nTracksSubEvt2 = 0;
   
   if (isGoodEvent())
   {
      UInt_t nTracks = picoDst->numberOfTracks();

      std::vector<int> idxTracksToRejectFromVtx;
      std::vector<int> allTracksForVtxFit;

      StThreeVectorF const pVtx = mPicoEvent->primaryVertex();

      for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
      {
         StPicoTrack* trk = picoDst->track(iTrack);
         if(!trk) continue;

         if(!isGoodForVertexFit(trk,pVtx)) idxTracksToRejectFromVtx.push_back(iTrack);
         else allTracksForVtxFit.push_back(iTrack);
      } // .. end tracks loop

      if(allTracksForVtxFit.size())
      {
        std::random_shuffle(allTracksForVtxFit.begin(),allTracksForVtxFit.end());

        int middleElement = static_cast<int>(allTracksForVtxFit.size()/2);
        std::vector<int> idxTracksToRejectFromVtxSub1(allTracksForVtxFit.size()-middleElement);
        std::vector<int> idxTracksToRejectFromVtxSub2(middleElement);

        std::copy(allTracksForVtxFit.begin(),allTracksForVtxFit.begin()+middleElement,idxTracksToRejectFromVtxSub2.begin());
        std::copy(allTracksForVtxFit.begin()+middleElement,allTracksForVtxFit.end(),idxTracksToRejectFromVtxSub1.begin());

        nTracksFullEvt = allTracksForVtxFit.size();
        nTracksSubEvt1 = allTracksForVtxFit.size() - idxTracksToRejectFromVtxSub1.size();
        nTracksSubEvt2 = allTracksForVtxFit.size() - idxTracksToRejectFromVtxSub2.size();

        for(size_t ij=0;ij<idxTracksToRejectFromVtx.size();++ij)
        {
          idxTracksToRejectFromVtxSub1.push_back(idxTracksToRejectFromVtx[ij]);
          idxTracksToRejectFromVtxSub2.push_back(idxTracksToRejectFromVtx[ij]);
        }

        kfVertex = mKfVertexFitter.primaryVertexRefit(picoDst,idxTracksToRejectFromVtx);
        kfVertexSubEvt1 = mKfVertexFitter.primaryVertexRefit(picoDst,idxTracksToRejectFromVtxSub1);
        kfVertexSubEvt2 = mKfVertexFitter.primaryVertexRefit(picoDst,idxTracksToRejectFromVtxSub2);
      }
   } //.. end of good event fill

   mVtxEvent.addEvent(*mPicoEvent,&kfVertex,&kfVertexSubEvt1,&kfVertexSubEvt2,
                            nTracksFullEvt,nTracksSubEvt1,nTracksSubEvt2);
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
