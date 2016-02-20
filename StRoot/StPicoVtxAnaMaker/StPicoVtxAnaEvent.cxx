#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoVtxAnaEvent.h"

StPicoVtxAnaEvent::StPicoVtxAnaEvent(char const* fileBaseName): mOutputFile(NULL), mTree(NULL)
{
   TString baseName(fileBaseName);
   mOutputFile = new TFile(Form("%s.kfVertex.root", fileBaseName), "RECREATE");
   mOutputFile->SetCompressionLevel(1);
   int BufSize = (int)pow(2., 16.);
   mTree = new TTree("kfEvent", "event information and kfVertex", BufSize);
   mTree->SetAutoSave(1000000); // autosave every 1 Mbytes

   mTree->Branch("mRunId   ", &mRunId   , "mRunId/I");
   mTree->Branch("mEventId ", &mEventId , "mEventId/I");
   mTree->Branch("mRefMult ", &mRefMult , "mRefMult/I");
   mTree->Branch("mGRefMult", &mGRefMult, "mGRefMult/I");
   mTree->Branch("mNTracks" , &mNTracks,  "mNTracks/I");
   mTree->Branch("mNTracksHft" , &mNTracksHft,  "mNTracksHft/I");
   mTree->Branch("mNTracksTop" , &mNTracksTop,  "mNTracksTop/I");
   mTree->Branch("mNTracksBottom" , &mNTracksBottom,  "mNTracksBottom/I");
   mTree->Branch("mNTracksRight" , &mNTracksRight,  "mNTracksRight/I");
   mTree->Branch("mNTracksLeft" , &mNTracksLeft,  "mNTracksLeft/I");
   mTree->Branch("mNTracksSubEvt1" , &mNTracksSubEvt1,  "mNTracksSubEvt1/I");
   mTree->Branch("mNTracksSubEvt2" , &mNTracksSubEvt2,  "mNTracksSubEvt2/I");
   mTree->Branch("mNTracksHftSubEvt1" , &mNTracksHftSubEvt1,  "mNTracksHftSubEvt1/I");
   mTree->Branch("mNTracksHftSubEvt2" , &mNTracksHftSubEvt2,  "mNTracksHftSubEvt2/I");
   
   mTree->Branch("mVx      ", &mVx      , "mVx/F");
   mTree->Branch("mVy      ", &mVy      , "mVy/F");
   mTree->Branch("mVz      ", &mVz      , "mVz/F");

   mTree->Branch("mKfVx    ", &mKfVx    , "mKfVx/F");
   mTree->Branch("mKfVy    ", &mKfVy    , "mKfVy/F");
   mTree->Branch("mKfVz    ", &mKfVz    , "mKfVz/F");

   mTree->Branch("mKfHftVx    ", &mKfHftVx    , "mKfHftVx/F");
   mTree->Branch("mKfHftVy    ", &mKfHftVy    , "mKfHftVy/F");
   mTree->Branch("mKfHftVz    ", &mKfHftVz    , "mKfHftVz/F");

   mTree->Branch("mKfTopVx    ", &mKfTopVx    , "mKfTopVx/F");
   mTree->Branch("mKfTopVy    ", &mKfTopVy    , "mKfTopVy/F");
   mTree->Branch("mKfTopVz    ", &mKfTopVz    , "mKfTopVz/F");

   mTree->Branch("mKfBottomVx    ", &mKfBottomVx    , "mKfBottomVx/F");
   mTree->Branch("mKfBottomVy    ", &mKfBottomVy    , "mKfBottomVy/F");
   mTree->Branch("mKfBottomVz    ", &mKfBottomVz    , "mKfBottomVz/F");

   mTree->Branch("mKfRightVx    ", &mKfRightVx    , "mKfRightVx/F");
   mTree->Branch("mKfRightVy    ", &mKfRightVy    , "mKfRightVy/F");
   mTree->Branch("mKfRightVz    ", &mKfRightVz    , "mKfRightVz/F");

   mTree->Branch("mKfLeftVx    ", &mKfLeftVx    , "mKfLeftVx/F");
   mTree->Branch("mKfLeftVy    ", &mKfLeftVy    , "mKfLeftVy/F");
   mTree->Branch("mKfLeftVz    ", &mKfLeftVz    , "mKfLeftVz/F");

   mTree->Branch("mKfSubEvt1Vx    ", &mKfSubEvt1Vx    , "mKfSubEvt1Vx/F");
   mTree->Branch("mKfSubEvt1Vy    ", &mKfSubEvt1Vy    , "mKfSubEvt1Vy/F");
   mTree->Branch("mKfSubEvt1Vz    ", &mKfSubEvt1Vz    , "mKfSubEvt1Vz/F");

   mTree->Branch("mKfSubEvt2Vx    ", &mKfSubEvt2Vx    , "mKfSubEvt2Vx/F");
   mTree->Branch("mKfSubEvt2Vy    ", &mKfSubEvt2Vy    , "mKfSubEvt2Vy/F");
   mTree->Branch("mKfSubEvt2Vz    ", &mKfSubEvt2Vz    , "mKfSubEvt2Vz/F");

   mTree->Branch("mKfHftSubEvt1Vx    ", &mKfHftSubEvt1Vx    , "mKfHftSubEvt1Vx/F");
   mTree->Branch("mKfHftSubEvt1Vy    ", &mKfHftSubEvt1Vy    , "mKfHftSubEvt1Vy/F");
   mTree->Branch("mKfHftSubEvt1Vz    ", &mKfHftSubEvt1Vz    , "mKfHftSubEvt1Vz/F");

   mTree->Branch("mKfHftSubEvt2Vx    ", &mKfHftSubEvt2Vx    , "mKfHftSubEvt2Vx/F");
   mTree->Branch("mKfHftSubEvt2Vy    ", &mKfHftSubEvt2Vy    , "mKfHftSubEvt2Vy/F");
   mTree->Branch("mKfHftSubEvt2Vz    ", &mKfHftSubEvt2Vz    , "mKfHftSubEvt2Vz/F");
}

void StPicoVtxAnaEvent::closeFile()
{
  mOutputFile->cd();
  mOutputFile->Write();
  mOutputFile->Close();
}

void StPicoVtxAnaEvent::addEvent(StPicoEvent const& picoEvent,
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
                 int nTrks,
                 int nTrksHft,
                 int nTrksTop,
                 int nTrksBottom,
                 int nTrksRight,
                 int nTrksLeft,
                 int nTrksSubEvt1,
                 int nTrksSubEvt2,
                 int nTrksHftSubEvt1,
                 int nTrksHftSubEvt2)
{
  mRunId    = picoEvent.runId();
  mEventId  = picoEvent.eventId();
  mRefMult  = picoEvent.refMult();
  mGRefMult = picoEvent.grefMult();

  mNTracks        = nTrks;
  mNTracksHft     = nTrksHft;
  mNTracksTop     = nTrksTop;
  mNTracksBottom  = nTrksBottom;
  mNTracksRight   = nTrksRight;
  mNTracksLeft    = nTrksLeft;
  mNTracksSubEvt1 = nTrksSubEvt1;
  mNTracksSubEvt2 = nTrksSubEvt2;
  mNTracksHftSubEvt1 = nTrksHftSubEvt1;
  mNTracksHftSubEvt2 = nTrksHftSubEvt2;

  mVx       = picoEvent.primaryVertex().x();
  mVy       = picoEvent.primaryVertex().y();
  mVz       = picoEvent.primaryVertex().z();

  mKfVx     = kfVtx.x();
  mKfVy     = kfVtx.y();
  mKfVz     = kfVtx.z();

  mKfHftVx     = kfHftVtx.x();
  mKfHftVy     = kfHftVtx.y();
  mKfHftVz     = kfHftVtx.z();

  mKfTopVx     = kfTopVtx.x();
  mKfTopVy     = kfTopVtx.y();
  mKfTopVz     = kfTopVtx.z();

  mKfBottomVx  = kfBottomVtx.x();
  mKfBottomVy  = kfBottomVtx.y();
  mKfBottomVz  = kfBottomVtx.z();

  mKfRightVx   = kfRightVtx.x();
  mKfRightVy   = kfRightVtx.y();
  mKfRightVz   = kfRightVtx.z();

  mKfLeftVx    = kfLeftVtx.x();
  mKfLeftVy    = kfLeftVtx.y();
  mKfLeftVz    = kfLeftVtx.z();

  mKfSubEvt1Vx     = kfSubEvt1.x();
  mKfSubEvt1Vy     = kfSubEvt1.y();
  mKfSubEvt1Vz     = kfSubEvt1.z();

  mKfSubEvt2Vx     = kfSubEvt2.x();
  mKfSubEvt2Vy     = kfSubEvt2.y();
  mKfSubEvt2Vz     = kfSubEvt2.z();

  mKfHftSubEvt1Vx     = kfHftSubEvt1.x();
  mKfHftSubEvt1Vy     = kfHftSubEvt1.y();
  mKfHftSubEvt1Vz     = kfHftSubEvt1.z();

  mKfHftSubEvt2Vx     = kfHftSubEvt2.x();
  mKfHftSubEvt2Vy     = kfHftSubEvt2.y();
  mKfHftSubEvt2Vz     = kfHftSubEvt2.z();

  mTree->Fill();
}
