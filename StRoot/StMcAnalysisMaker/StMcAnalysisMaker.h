#ifndef ST_MCANALYSISMAKER_H
#define ST_MCANALYSISMAKER_H

#include <vector>
#include "TString.h"

#include "StChain/StMaker.h"

class TFile;
class TNtuple;

class StMcTrack;
class StTrack;
class StGlobalTrack;
class StAssociationMaker;
class StMcEvent;
class StEvent;
class StDedxPidTraits;

class StMcAnalysisMaker : public StMaker
{
private:
   TString mOutfileName;
   std::vector<float> firedTriggersIndices;
   double mField; //.. magnetic field

   TFile* mFile;
   TNtuple* mTracks;
   TNtuple* mEventCount; //.. For counting purposes

   StMcEvent* mMcEvent;
   StEvent* mEvent;
   StAssociationMaker* mAssoc;

   StTrack const* findPartner(StMcTrack*, int&) const;
   StMcTrack const* findPartner(StGlobalTrack*, int&) const;
   StDedxPidTraits const* findDedxPidTraits(StTrack const*) const;
   int getNHitsDedx(StTrack const*) const;

   bool passTrigger();
   int  fillEventCounts(float nRTracks = -1, float nMcTracks = -1);
   int  fillTracks(int& nRTracks, int& nMcTracks);
   void fillMcTrack(float* array,int& idx,StMcTrack const*);
   void fillRcTrack(float* array,int& idx,StMcTrack const*,StTrack const*,int const ncom);
   unsigned int  getHftTruth(StMcTrack const*,StTrack const*) const;
   void getDca(StTrack const*,float& dca, float& dcaXY, float& dcaZ) const;

   bool isGoodMcTrack(StMcTrack const*) const;

public:
   StMcAnalysisMaker(TString name);

   int Init();
   int Make();
   int Finish();

   ClassDef(StMcAnalysisMaker, 0)
};
#endif
