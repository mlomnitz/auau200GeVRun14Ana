#include "TFile.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TSystem.h"

#include "StarClassLibrary/StParticleDefinition.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StarClassLibrary/StParticleTypes.hh"
#include "StBFChain/StBFChain.h"

#include "StEvent/StEventTypes.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTrackGeometry.h"
#include "StEvent/StTrackNode.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StTrackTopologyMap.h"
#include "StEvent/StEventSummary.h"
#include "StEvent/StBTofCollection.h"
#include "StEvent/StBTofHeader.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StTpcDedxPidAlgorithm.h"
#include "StEventUtilities/StuRefMult.hh"

#include "StMcEvent/StMcEventTypes.hh"
#include "StMcEvent/StMcContainers.hh"
#include "StMcEvent/StMcVertex.hh"
#include "StMcEvent/StMcEvent.hh"

#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"

#include "StMcAnaCuts.h"
#include "StMcAnalysisMaker.h"

ClassImp(StMcAnalysisMaker);

StMcAnalysisMaker::StMcAnalysisMaker(TString name): StMaker(name.Data()), mFile(NULL),
   mTracks(NULL), mEventCount(NULL), mMcEvent(NULL), mEvent(NULL), mAssoc(NULL)
{
   LOG_INFO << "StMcAnalysisMaker() - DONE" << endm;
}
//__________________________________
int StMcAnalysisMaker::Init()
{
   StBFChain *bfChain = (StBFChain *) StMaker::GetChain();

   if (!bfChain) return kStFatal;
   TString fileName(gSystem->BaseName(bfChain->GetFileOut().Data()));
   fileName = fileName.ReplaceAll(".event.root", "");
   fileName = fileName.ReplaceAll(".geant.root", "");
   fileName = fileName.ReplaceAll(".MuDst.root", "");

   if (!fileName.Length()) fileName = "mcAnalysis";
   fileName = fileName.ReplaceAll(".root", "");

   mFile = new TFile(Form("%s.McAna.root", fileName.Data()), "recreate");
   assert(mFile && !mFile->IsZombie());

   mAssoc = (StAssociationMaker*)GetMaker("StAssociationMaker");
   if (!mAssoc)
   {
      LOG_ERROR << "Could not get StAssociationMaker" << endm;
      return kStErr;
   }

   mEventCount = new TNtuple("eventCount", "eventCount", "runId:eventId:mcVx:mcVy:mcVz:vx:vy:vz:vzVpd:"
                             "posRefMult:negRefMult:zdc:bbc:nMcTracks:nRTracks:magField:t0:t1:t2:t3:t4:t5");

   mTracks = new TNtuple("tracks", "", "pt:p:eta:y:phi:geantId:eventGenLabel:startVtxX:startVtxY:startVtxZ:stopVtxX:stopVtxY:stopVtxZ:" // MC
                         "rPt:rEta:rPhi:nFit:nMax:nCom:nDedx:nDedx2:dedx:nSigPi:nSigK:dca:dcaXY:dcaZ:hftTopo:hftTruth"); // global

   LOG_INFO << "Init() - DONE" << endm;

   return kStOk;
}

//__________________________________
int StMcAnalysisMaker::Make()
{
   StMcEvent* mcEvent = (StMcEvent*)GetDataSet("StMcEvent");

   if (!mcEvent)
   {
      LOG_WARN << "No StMcEvent" << endm;
      return kStWarn;
   }

   StEvent* mEvent = (StEvent*)GetDataSet("StEvent");
   if (!mEvent)
   {
      LOG_WARN << "No StEvent" << endm;
      return kStWarn;
   }

   mField = (float)mEvent->summary()->magneticField();

   // Fill
   int nRTracks = -1;
   int nMcTracks = -1;

   int fillTracksStatus = kStOk;
   if (passTrigger())
   {
      fillTracksStatus = fillTracks(nRTracks, nMcTracks);
   }
   else
   {
      LOG_INFO << "No interesting triggers. Counting event then skipping." << endm;
   }

   int fillEventCountStatus = fillEventCounts((float)nRTracks, (float)nMcTracks);

   return fillTracksStatus && fillEventCountStatus;
}

int StMcAnalysisMaker::fillTracks(int& nRTracks, int& nMcTracks)
{
   nRTracks = 0;
   nMcTracks = 0;

   LOG_INFO << "Filling " << mMcEvent->tracks().size() << " mcTracks..." << "\n";

   for (unsigned int iTrk = 0;  iTrk < mMcEvent->tracks().size(); ++iTrk)
   {
      StMcTrack* const mcTrack = mMcEvent->tracks()[iTrk];

      if (!isGoodMcTrack(mcTrack)) continue;
      ++nMcTracks;

      int nCommonHits = 0;
      StTrack const* const rcTrack = findPartner(mcTrack, nCommonHits);

      float array[220];
      for (int ii = 0; ii < 220; ++ii) array[ii] = -999;
      int idx = 0;

      fillMcTrack(array, idx, mcTrack);
      if (rcTrack)
      {
         ++nRTracks;
         fillRcTrack(array, idx, mcTrack,rcTrack, nCommonHits);
      }

      mTracks->Fill(array);
   }

   LOG_INFO << endm;

   return kStOk;
}

void StMcAnalysisMaker::fillMcTrack(float* array, int& idx, StMcTrack const* const mcTrk)
{
   array[idx++] = mcTrk->pt();
   array[idx++] = mcTrk->momentum().mag();
   array[idx++] = mcTrk->pseudoRapidity();
   array[idx++] = mcTrk->rapidity();
   array[idx++] = mcTrk->momentum().phi();
   array[idx++] = mcTrk->geantId();
   array[idx++] = mcTrk->eventGenLabel();
   array[idx++] = mcTrk->startVertex()->position().x();
   array[idx++] = mcTrk->startVertex()->position().y();
   array[idx++] = mcTrk->startVertex()->position().z();
   array[idx++] = mcTrk->stopVertex()->position().x();
   array[idx++] = mcTrk->stopVertex()->position().y();
   array[idx++] = mcTrk->stopVertex()->position().z();
}

void StMcAnalysisMaker::fillRcTrack(float* array, int& idx, StMcTrack const* const mcTrack,StTrack const* const rcTrack, int const ncom)
{
   if (!rcTrack) return;

   array[idx++] = rcTrack->geometry()->momentum().perp();
   array[idx++] = rcTrack->geometry()->momentum().pseudoRapidity();
   array[idx++] = rcTrack->geometry()->momentum().phi();
   array[idx++] = rcTrack->fitTraits().numberOfFitPoints(kTpcId);
   array[idx++] = rcTrack->numberOfPossiblePoints(kTpcId);
   array[idx++] = ncom;

   // dedx info
   float nDedxPts = -9999;
   float dedx = -9999;
   float nSigPi = -9999;
   float nSigK = -9999;
   static StTpcDedxPidAlgorithm aplus(McAnaCuts::dedxMethod);
   static StPionPlus* Pion = StPionPlus::instance();
   static StKaonPlus* Kaon = StKaonPlus::instance();
   StParticleDefinition const* prtcl = rcTrack->pidTraits(aplus);

   if (prtcl)
   {
      nDedxPts = aplus.traits()->numberOfPoints();
      dedx = aplus.traits()->mean();
      nSigPi = aplus.numberOfSigma(Pion);
      nSigK = aplus.numberOfSigma(Kaon);
   }

   array[idx++] = nDedxPts;
   array[idx++] = getNHitsDedx(rcTrack);
   array[idx++] = dedx;
   array[idx++] = nSigPi;
   array[idx++] = nSigK;

   unsigned int hftTruth = 999;
   float dca = -999.;
   float dcaXY = -999.;
   float dcaZ = -999.;

   hftTruth = getHftTruth(mcTrack, rcTrack);
   getDca(rcTrack, dca, dcaXY, dcaZ);

   array[idx++] = dca;
   array[idx++] = dcaXY;
   array[idx++] = dcaZ;
   array[idx++] = rcTrack->topologyMap().data(0) >> 1 & 0x7F;
   array[idx++] = hftTruth;
}

bool StMcAnalysisMaker::isGoodMcTrack(StMcTrack const* const mcTrack) const
{
   return mcTrack->geantId() == McAnaCuts::geantId && mcTrack->startVertex()->position().perp() < McAnaCuts::mcTrackStartVtxR;
}

int StMcAnalysisMaker::fillEventCounts(float nRTracks, float nMcTracks)
{
   float vars[30];

   float vpdVz = -999;
   StBTofHeader* tofheader = 0;
   if (mEvent->btofCollection())  tofheader = mEvent->btofCollection()->tofHeader();
   if (tofheader) vpdVz = tofheader->vpdVz();

   int iVar = 0;
   vars[iVar++] = (float)mEvent->runId();
   vars[iVar++] = (float)mEvent->id();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().x();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().y();
   vars[iVar++] = (float)mMcEvent->primaryVertex()->position().z();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().x();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().y();
   vars[iVar++] = (float)mEvent->primaryVertex()->position().z();
   vars[iVar++] = vpdVz;
   vars[iVar++] = (float)uncorrectedNumberOfPositivePrimaries(*mEvent);
   vars[iVar++] = (float)uncorrectedNumberOfNegativePrimaries(*mEvent);
   vars[iVar++] = (float)mEvent->runInfo()->zdcCoincidenceRate();
   vars[iVar++] = (float)mEvent->runInfo()->bbcCoincidenceRate();
   vars[iVar++] = nMcTracks;
   vars[iVar++] = nRTracks;
   vars[iVar++] = (float)mEvent->summary()->magneticField() / 10;
   vars[iVar++] = firedTriggersIndices.at(0);
   vars[iVar++] = firedTriggersIndices.at(1);
   vars[iVar++] = firedTriggersIndices.at(2);
   vars[iVar++] = firedTriggersIndices.at(3);
   vars[iVar++] = firedTriggersIndices.at(4);
   vars[iVar++] = firedTriggersIndices.at(5);

   mEventCount->Fill(vars);

   return kStOk;
}

bool StMcAnalysisMaker::passTrigger()
{
   LOG_INFO << "Checking triggers..." << endm;
   bool interesting_event = false;

   const StTriggerId* st_trgid = mEvent->triggerIdCollection()->nominal();

   for (unsigned int ii = 0; ii < firedTriggersIndices.size(); ++ii)
   {
      firedTriggersIndices[ii] = -1;
   }

   // Fill interesting triggers
   LOG_INFO << "Interesting fired triggers: " << "\n";

   for (unsigned int ii = 0; ii < st_trgid->maxTriggerIds(); ++ii)
   {
      unsigned int id = st_trgid->triggerId(ii);

      int trgIndex = -1;

      for (unsigned int jj = 0; jj < McAnaCuts::interesting_triggers.size(); ++jj)
      {
         if (McAnaCuts::interesting_triggers[jj] == id)
         {
            trgIndex = jj;
            interesting_event = true;
            LOG_INFO << id << " ";
            break;
         }
      }

      if (trgIndex != -1) firedTriggersIndices[trgIndex] = 1.0;
   }

   LOG_INFO << endm;

   return interesting_event;
}

StTrack const* StMcAnalysisMaker::findPartner(StMcTrack* mcTrack, int& maxCommonTpcHits) const
{
   //..StMcTrack find partner from the StTracks
   pair<mcTrackMapIter, mcTrackMapIter> p = mAssoc->mcTrackMap()->equal_range(mcTrack);

   const StTrack* maxTrack = 0;
   maxCommonTpcHits = 0;
   for (mcTrackMapIter k = p.first; k != p.second; ++k)
   {
      int commonTpcHits = k->second->commonTpcHits();
      const StTrack* track = k->second->partnerTrack()->node()->track(global);//should be global
      if (track && commonTpcHits > maxCommonTpcHits)
      {
         maxTrack = track;
         maxCommonTpcHits = commonTpcHits;
      }
   }
   return maxTrack;
}

void StMcAnalysisMaker::getDca(StTrack const* const rcTrack, float& dca, float& dcaXY, float& dcaZ) const
{
   StPhysicalHelixD helix = rcTrack->geometry()->helix();
   dca = helix.distance(mEvent->primaryVertex()->position());
   dcaXY = helix.geometricSignedDistance(mEvent->primaryVertex()->position().x(), mEvent->primaryVertex()->position().y());

   StThreeVectorF dcaPoint = helix.at(helix.pathLength(mEvent->primaryVertex()->position()));
   dcaZ = dcaPoint.z() - mEvent->primaryVertex()->position().z();
}

unsigned int StMcAnalysisMaker::getHftTruth(StMcTrack const* const mcTrack, StTrack const* const rcTrack) const
{
   bool istTruth = true;
   bool pxlTruth1 = true;
   bool pxlTruth2 = true;

   StPtrVecHit rcIstHits = rcTrack->detectorInfo()->hits(kIstId);
   int nRcIstHits = (int)rcIstHits.size();
   for (int iHit = 0; iHit < nRcIstHits; ++iHit)
   {
      if (rcIstHits[iHit]->idTruth() != mcTrack->key())
      {
         istTruth = false;
         break;
      }
   }

   StPtrVecHit rcPxlHits = rcTrack->detectorInfo()->hits(kPxlId);
   int nRcPxlHits = (int)rcPxlHits.size();
   for (int iHit = 0; iHit < nRcPxlHits; ++iHit)
   {
      if (rcPxlHits[iHit]->idTruth() != mcTrack->key())
      {
         StThreeVectorF pos = rcPxlHits[iHit]->position();

         float const R = pow(pos.x(), 2.0) + pow(pos.y(), 2.0);
         if (R > 3.5 * 3.5) pxlTruth2 = false;
         else pxlTruth1 = false;
      }
   }

   unsigned int hftTruth = 0;
   if (pxlTruth1) hftTruth |= (1 << 0);
   if (pxlTruth2) hftTruth |= (1 << 1);
   if (istTruth)  hftTruth |= (1 << 2);

   return hftTruth;
}

StMcTrack const* StMcAnalysisMaker::findPartner(StGlobalTrack* rcTrack, int& maxCommonTpcHits) const
{
   //.. StGlobalTracks find partner from StMcTracks.
   //.. See example from StRoot/StMcAnalysisMaker
   pair<rcTrackMapIter, rcTrackMapIter> p = mAssoc->rcTrackMap()->equal_range(rcTrack);

   const StMcTrack* maxTrack = 0;
   maxCommonTpcHits = 0;
   for (rcTrackMapIter k = p.first; k != p.second; ++k)
   {
      int commonTpcHits = k->second->commonTpcHits();

      const StMcTrack* track = k->second->partnerMcTrack();

      if (track && commonTpcHits > maxCommonTpcHits)
      {
         maxTrack = track;
         maxCommonTpcHits = commonTpcHits;
      }
   }
   return maxTrack;
}

int StMcAnalysisMaker::Finish()
{
   mFile->cd();
   mFile->Write();
   mFile->Close();
   return kStOk;
}

StDedxPidTraits const* StMcAnalysisMaker::findDedxPidTraits(StTrack const* const rcTrack) const
{
   StDedxPidTraits* pid = 0;
   StPtrVecTrackPidTraits traits = rcTrack->pidTraits(kTpcId);

   for (unsigned int ii = 0; ii < traits.size(); ++ii)
   {
      pid = dynamic_cast<StDedxPidTraits*>(traits[ii]);
      if (pid && pid->method() == McAnaCuts::dedxMethod) break;
   }

   return pid;
}

int StMcAnalysisMaker::getNHitsDedx(StTrack const* const t) const
{
   int ndedx = -1;
   StPtrVecTrackPidTraits pidTraits = t->pidTraits(kTpcId);

   if (pidTraits.size())
   {
      StDedxPidTraits* pid;
      for (unsigned int ii = 0; ii < pidTraits.size(); ++ii)
      {
         pid = dynamic_cast<StDedxPidTraits*>(pidTraits[ii]);

         if (pid && (pid->method() == McAnaCuts::dedxMethod))
         {
            ndedx = pid->numberOfPoints();            //number of dedx hits
            break;
         }
      }
   }

   return ndedx;
}
