#ifndef StPicoVtxAnaEvent__h
#define StPicoVtxAnaEvent__h

/* **************************************************
 *  A class to save event information with KF vertex.
 *
 *  Authors:  **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */

#include "StarClassLibrary/StThreeVectorF.hh"

class TTree;
class TFile;
class StPicoEvent;

class StPicoVtxAnaEvent
{
public:
   StPicoVtxAnaEvent(char const*);
   void addEvent(StPicoEvent const&,
                 StThreeVectorF const& kfVtx, 
                 StThreeVectorF const& kfHftVtx,
                 StThreeVectorF const& kfTopVtx,
                 StThreeVectorF const& kfBottomVtx,
                 StThreeVectorF const& kfRightVtx,
                 StThreeVectorF const& kfLeftVtx,
                 StThreeVectorF const& kfSubEvt1,
                 StThreeVectorF const& kfSubEvt2,
                 StThreeVectorF const& kfHftSubEvt1,
                 StThreeVectorF const& kfHftSubEvt2,
                 StThreeVectorF const* kfHftPxlSecVtx,
                 StThreeVectorF const* kfHftHighPtPxlSecVtx,
                 int nTrks=0,
                 int nTrksHft=0,
                 int nTrksTop=0,
                 int nTrksBottom=0,
                 int nTrksRight=0,
                 int nTrksLeft=0,
                 int nTrksSubEvt1=0,
                 int nTrksSubEvt2=0,
                 int nTrksHftSubEvt1=0,
                 int nTrksHftSubEvt2=0,
                 int* nTrksHftPxlSec = NULL,
                 int* nTrksHftHighPtPxlSec = NULL);
   void closeFile();

private:

   int   mRunId;
   int   mEventId;
   int   mRefMult;
   int   mNTracks;
   int   mNTracksHft;
   int   mNTracksTop;
   int   mNTracksBottom;
   int   mNTracksRight;
   int   mNTracksLeft;
   int   mNTracksSubEvt1;
   int   mNTracksSubEvt2;
   int   mNTracksHftSubEvt1;
   int   mNTracksHftSubEvt2;
   int   mNTracksHftPxlSec[10];
   int   mNTracksHftHighPtPxlSec[10];
   int   mGRefMult;

   float mVx;
   float mVy;
   float mVz;

   float mKfVx;
   float mKfVy;
   float mKfVz;

   float mKfHftVx;
   float mKfHftVy;
   float mKfHftVz;

   float mKfTopVx;
   float mKfTopVy;
   float mKfTopVz;

   float mKfBottomVx;
   float mKfBottomVy;
   float mKfBottomVz;

   float mKfRightVx;
   float mKfRightVy;
   float mKfRightVz;

   float mKfLeftVx;
   float mKfLeftVy;
   float mKfLeftVz;

   float mKfSubEvt1Vx;
   float mKfSubEvt1Vy;
   float mKfSubEvt1Vz;

   float mKfSubEvt2Vx;
   float mKfSubEvt2Vy;
   float mKfSubEvt2Vz;

   float mKfHftSubEvt1Vx;
   float mKfHftSubEvt1Vy;
   float mKfHftSubEvt1Vz;

   float mKfHftSubEvt2Vx;
   float mKfHftSubEvt2Vy;
   float mKfHftSubEvt2Vz;

   float mKfHftPxlSecVx[10];
   float mKfHftPxlSecVy[10];
   float mKfHftPxlSecVz[10];

   float mKfHftHighPtPxlSecVx[10];
   float mKfHftHighPtPxlSecVy[10];
   float mKfHftHighPtPxlSecVz[10];

   TFile* mOutputFile;
   TTree* mTree;
};

#endif
