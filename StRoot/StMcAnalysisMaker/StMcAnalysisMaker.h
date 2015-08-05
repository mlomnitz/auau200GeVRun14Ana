#ifndef ST_MCANALYSISMAKER_H
#define ST_MCANALYSISMAKER_H

#include <vector>
#include "TString.h"

#include "StChain/StMaker.h"

class TFile;
class TNtuple;
class TH3F;

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
   bool mFillTpcHitsNtuple;

   TFile* mFile;
   TNtuple* mTracks;
   TNtuple* mEventCount; //.. For counting purposes
   TNtuple* mTpcNtuple;

   TH3F* hTpcHitsDiffXVsPadrowVsSector;                                                                                                                                                                  
   TH3F* hTpcHitsDiffYVsPadrowVsSector;
   TH3F* hTpcHitsDiffZVsPadrowVsSector;

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

   void fillTpcNtuple(StMcTrack const* const,StTrack const* const);

public:
   StMcAnalysisMaker(TString name);

   int Init();
   int Make();
   int Finish();

   void fillTpcHitsNtuple(bool t=true);
   ClassDef(StMcAnalysisMaker, 0)
};

inline void StMcAnalysisMaker::fillTpcHitsNtuple(bool t){ mFillTpcHitsNtuple=t;}
#endif
