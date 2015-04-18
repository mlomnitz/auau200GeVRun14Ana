#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StRoot/StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoD0AnaMaker.h"
#include "StPicoD0AnaHists.h"
#include "StAnaCuts.h"

ClassImp(StPicoD0AnaMaker)

StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name, TString const inputFilesList,
                                   TString const outFileBaseName, StPicoDstMaker* picoDstMaker):
   StMaker(name), mPicoDstMaker(picoDstMaker), mPicoD0Event(NULL),
   mInputFilesList(inputFilesList), mOutFileBaseName(outFileBaseName), mChain(NULL), mEventCounter(0),
   mHists(NULL)
{}

Int_t StPicoD0AnaMaker::Init()
{
   mPicoD0Event = new StPicoD0Event();

   mChain = new TChain("T");
   std::ifstream listOfFiles(mInputFilesList.Data());
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoD0AnaMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
      return kStErr;
   }

   mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   // -------------- USER VARIABLES -------------------------
   mHists = new StPicoD0AnaHists(mOutFileBaseName);

   return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
   /*  */
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
   mHists->closeFile();
   return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
   readNextEvent();

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();

   if (!picoDst)
   {
      LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if (mPicoD0Event->runId() != picoDst->event()->runId() ||
         mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
      LOG_ERROR << " StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!" << endm;
      LOG_ERROR << " StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync." << endm;
      exit(1);
   }

   // -------------- USER ANALYSIS -------------------------

   if (isGoodEvent(picoDst->event()))
   {
      mHists->addEvent(picoDst->event());

      TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();

      StThreeVectorF const pVtx = picoDst->event()->primaryVertex();

      for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
      {
         // this is an example of how to get the kaonPion pairs and their corresponsing tracks
         StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);

         if (!isGoodPair(kp)) continue;

         StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
         StPicoTrack const* pion = picoDst->track(kp->pionIdx());

         if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;
         if (!isTpcPion(pion)) continue;

         bool tpcKaon = isTpcKaon(kaon);
         float kBeta = getTofBeta(kaon,*pVtx);
         bool tofAvailable = kBeta>0;
         bool tofKaon = tofAvailable && isTofKaon(kaon,kBeta);

         if (tpcKaon || tofKaon)
         {
            bool unlike = kaon->charge() * pion->charge() < 0 ? true : false;
            mHists->addKaonPion(kp, unlike, tpcKaon, ((tofAvailable && tofKaon) || (!tofAvailable && tpcKaon)));
         }

      } // end of kaonPion loop
   } // end of isGoodEvent

   return kStOK;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodEvent(StPicoEvent const * const picoEvent) const
{
   return picoEvent->triggerWord() & anaCuts::triggerWord;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
   return trk->nHitsFit() >= anaCuts::nHitsFit;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const * const trk) const
{
   return fabs(trk->nSigmaPion()) < anaCuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const * const trk) const
{
   return fabs(trk->nSigmaKaon()) < anaCuts::nSigmaKaon;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodPair(StKaonPion const* const kp) const
{
   return cos(kp->pointingAngle()) > anaCuts::cosTheta &&
          kp->pionDca() > anaCuts::pDca && kp->kaonDca() > anaCuts::kDca &&
          kp->dcaDaughters() < anaCuts::dcaDaughters;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
   bool tofKaon = false;

   if(beta>0)
   {
     double ptot = trk->dcaGeometry().momentum().mag();
     float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
     tofKaon = fabs(1/beta - 1/beta_k) < anaCuts::kTofBetaDiff ? true : false;
   }
   
   return tofKaon;
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
   StDcaGeometry dcaG = trk->dcaGeometry();
   double ptot = dcaG.momentum().mag();
   StPhysicalHelixD helix = dcaG.helix();

   int index2tof = trk->bTofPidTraitsIndex();

   float beta = std::numeric_limits<float>::quiet_NaN();

   if(index2tof >= 0)
   {
      StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

      if(tofPid)
      {
         beta = tofPid->btofBeta();

         if (beta < 1e-4)
         {
            StThreeVectorF const btofHitPos = tofPid->btofHitPos();

            float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
         }
         else
         {
           beta = std::numeric_limits<float>::quiet_NaN();
         }
      }
   }

   return beta;
}
