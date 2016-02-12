#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
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
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "kfEvent.h"

#include "StMemStat.h"

ClassImp(StPicoD0AnaMaker)

StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name, TString const inputFilesList, char const * kfFileList,
                                   TString const outFileBaseName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil):
   StMaker(name), mPicoDstMaker(picoDstMaker), mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil), mKfEvent(NULL),
   mInputFilesList(inputFilesList), mKfFileList(kfFileList), mOutFileBaseName(outFileBaseName),
   mChain(NULL), mKfChain(NULL), mEventCounter(0), mFillQaHists(false), mHists(NULL)
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

   // -------------Next is include KfVertex tree
   // mKfEvent = new kfEvent();
   mKfChain = new TChain("kfEvent");
   std::ifstream listOfKfFiles;
   listOfKfFiles.open(mKfFileList);
   if (listOfKfFiles.is_open())
   {
      std::string kffile;
      while (getline(listOfKfFiles, kffile))
      {
         LOG_INFO << "StPicoD0AnaMaker - Adding kfVertex tree:" << kffile << endm;
         mKfChain->Add(kffile.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoD0AnaMaker - Could not open list of corresponding kfEvent files. ABORT!" << endm;
      return kStErr;
   }
   mKfEvent = new kfEvent(mKfChain);

   // -------------- USER VARIABLES -------------------------
   mHists = new StPicoD0AnaHists(mOutFileBaseName,mFillQaHists);

   return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
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
   StMemStat mem;
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

   if (mPicoD0Event->runId() != picoDst->event()->runId() ||  mPicoD0Event->runId() != mKfEvent->mRunId  ||
         mPicoD0Event->eventId() != picoDst->event()->eventId() || mPicoD0Event->eventId() != mKfEvent->mEventId)
   {
      LOG_ERROR << " StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!" << "\n";
      LOG_ERROR << " StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event and KfEvent are not in sync." << endm;
      exit(1);
   }

   // -------------- USER ANALYSIS -------------------------

   mGRefMultCorrUtil->init(picoDst->event()->runId());

   if (!mGRefMultCorrUtil)
   {
     LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
     return kStWarn;
   }

   if (mGRefMultCorrUtil->isBadRun(picoDst->event()->runId()))
   {
     //cout<<"This is a bad run from mGRefMultCorrUtil! Skip! " << endl;
     return kStOK;
   }

   mHists->addEventBeforeCut(picoDst->event());

   if (picoDst->event()->primaryVertex().x() != mKfEvent->mVx)
   {
      LOG_ERROR << " StPicoMixedEventMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!" << "\n";
      LOG_ERROR << " StPicoMixedEventMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoDst and KfEvent vertex are not in sync." << endm;
      exit(1);
   }
   StThreeVectorF const kfVtx(mKfEvent->mKfVx, mKfEvent->mKfVy, mKfEvent->mKfVz);

   if (isGoodEvent(picoDst->event(),kfVtx))
   {
      TClonesArray const* aKaonPion = mPicoD0Event->kaonPionArray();
      if (aKaonPion->GetEntries()) mHists->addEvent(picoDst->event());

      mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(), kfVtx.z(), picoDst->event()->ZDCx()) ;

      int centrality  = mGRefMultCorrUtil->getCentralityBin9();
      const double reweight = mGRefMultCorrUtil->getWeight();
      const double refmultCor = mGRefMultCorrUtil->getRefMultCorr();
      mHists->addCent(refmultCor, centrality, reweight, kfVtx.z());

      //Basiclly add some QA plots
      UInt_t nTracks = picoDst->numberOfTracks();

      if(mFillQaHists)
      {
        for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
        {
          StPicoTrack const* trk = picoDst->track(iTrack);
          if (!trk) continue;
          StPhysicalHelixD helix = trk->helix();
          float dca = float(helix.geometricSignedDistance(kfVtx));
          StThreeVectorF momentum = trk->gMom(kfVtx, picoDst->event()->bField());

          if (!isGoodQaTrack(trk, momentum, dca)) continue;

          StThreeVectorF dcaPoint = helix.at(helix.pathLength(kfVtx.x(), kfVtx.y()));
          float dcaZ = dcaPoint.z() - kfVtx.z();
          double dcaXy = helix.geometricSignedDistance(kfVtx.x(), kfVtx.y());

         bool tpcPion = isTpcPion(trk);
         bool tpcKaon = isTpcKaon(trk);
         float pBeta = getTofBeta(trk, kfVtx);
         float kBeta = pBeta;
         bool pTofAvailable = !isnan(pBeta) && pBeta > 0;
         bool kTofAvailable = !isnan(kBeta) && kBeta > 0;
         bool tofPion = isTofPion(trk, pBeta, kfVtx);
         bool tofKaon = isTofKaon(trk, kBeta, kfVtx);

         bool goodPion = (pTofAvailable && tofPion && tpcPion) || (!pTofAvailable && tpcPion);//Always require TPC
         bool goodKaon = (kTofAvailable && tofKaon && tpcKaon) || (!kTofAvailable && tpcKaon);
         // bool goodKaon = (momentum.perp() <= 1.6 && kTofAvailable && tofKaon && tpcKaon) || (momentum.perp() > 1.6 && tpcKaon);//strict Kaon pid

         if (trk  && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon ))
         {
            mHists->addDcaPtCent(dca, dcaXy, dcaZ, goodPion, goodKaon, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), kfVtx.z(), picoDst->event()->ZDCx() / 1000.); //add Dca distribution
         }
         if (trk  && fabs(dca) < 1.5 && (goodPion || goodKaon ))
         {
            mHists->addTpcDenom1(goodPion, goodKaon, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), kfVtx.z(), picoDst->event()->ZDCx() / 1000.); //Dca cut on 1.5cm, add Tpc Denominator
         }
         if (trk && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon ) && fabs(dcaXy) < 1. && fabs(dcaZ) < 1.)
         {
            mHists->addHFTNumer1(goodPion, goodKaon, momentum.perp(), centrality,  momentum.pseudoRapidity(), momentum.phi(), kfVtx.z(), picoDst->event()->ZDCx() / 1000.); //Dca cut on 1.5cm, add HFT Numerator
         }
        } // .. end tracks loop
      }

      for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
      {
         StKaonPion const* kp = (StKaonPion*)aKaonPion->UncheckedAt(idx);

         if (!isGoodPair(kp)) continue;

         StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
         StPicoTrack const* pion = picoDst->track(kp->pionIdx());

         if (!isGoodTrack(kaon, kfVtx) || !isGoodTrack(pion, kfVtx)) continue;

         // PID
         if(!isTpcPion(pion) || !isTpcKaon(kaon)) continue;
         float pBeta = getTofBeta(pion, &kfVtx);
         float kBeta = getTofBeta(kaon, &kfVtx);
         bool pTofAvailable = !isnan(pBeta) && pBeta > 0;
         bool kTofAvailable = !isnan(kBeta) && kBeta > 0;
         bool tofPion = pTofAvailable ? isTofPion(pion, pBeta, kfVtx) : true;//this is bybrid pid, not always require tof
         bool tofKaon = kTofAvailable ? isTofKaon(kaon, kBeta, kfVtx) : true;//this is bybrid pid, not always require tof
         bool tof = tofPion && tofKaon;

         bool unlike = kaon->charge() * pion->charge() < 0 ? true : false;
         bool isd0 = (unlike && pion->charge() > 0 && kaon->charge() < 0) ? true : false;//d0--> pi+ k-
         bool isPosLike = (!unlike && pion->charge() > 0 && kaon->charge() > 0 ) ? true : false;//positive like sign  pi+ k+

         int pionRegion = trkHalf(pion, kfVtx);// -1 -> left half ,  +1 -> right half side
         int kaonRegion = trkHalf(kaon, kfVtx);// -1 -> left half ,  +1 -> right half side
         bool isLeft  = (pionRegion+kaonRegion) < 0  ? true : false;
         bool isRight = (pionRegion+kaonRegion) > 0  ? true : false;
         bool isMixed = (pionRegion+kaonRegion) == 0 ? true : false;

         mHists->addKaonPion(kp, unlike, isd0, isPosLike, isLeft, isRight, isMixed, true, tof, centrality, reweight);

      } // end of kaonPion loop
   } // end of isGoodEvent

   return kStOK;
}
//-----------------------------------------------------------------------------
int StPicoD0AnaMaker::getD0PtIndex(StKaonPion const* const kp) const
{
   for (int i = 0; i < anaCuts::nPtBins; i++)
   {
      if ((kp->pt() >= anaCuts::PtBinsEdge[i]) && (kp->pt() < anaCuts::PtBinsEdge[i + 1]))
         return i;
   }
   return anaCuts::nPtBins - 1;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodEvent(StPicoEvent const* const picoEvent, StThreeVectorF const& kfVtx) const
{
   return (picoEvent->triggerWord() & anaCuts::triggerWord) &&
          fabs(kfVtx.z()) < anaCuts::vz &&
          fabs(kfVtx.z() - picoEvent->vzVpd()) < anaCuts::vzVpdVz &&
          !(fabs(kfVtx.x()) < anaCuts::Verror && fabs(kfVtx.y()) < anaCuts::Verror && fabs(kfVtx.z()) < anaCuts::Verror) &&
          sqrt(TMath::Power(kfVtx.x(), 2) + TMath::Power(kfVtx.y(), 2)) <=  anaCuts::Vrcut;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodQaTrack(StPicoTrack const* const trk, StThreeVectorF const& momentum, double const dca) const
{
   return trk->gPt() > anaCuts::qaGPt && trk->nHitsFit() >= anaCuts::qaNHitsFit && fabs(momentum.pseudoRapidity()) <= anaCuts::Eta;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
   StThreeVectorF mom = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField());


   return mom.perp() > anaCuts::minPt &&
          trk->nHitsFit() >= anaCuts::nHitsFit &&
          fabs(mom.pseudoRapidity()) <= anaCuts::Eta;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const* const trk) const
{
   return fabs(trk->nSigmaPion()) < anaCuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const* const trk) const
{
   return fabs(trk->nSigmaKaon()) < anaCuts::nSigmaKaon;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodPair(StKaonPion const* const kp) const
{
   int tmpIndex = getD0PtIndex(kp);
   return cos(kp->pointingAngle()) > anaCuts::cosTheta[tmpIndex] &&
          kp->pionDca() > anaCuts::pDca[tmpIndex] && kp->kaonDca() > anaCuts::kDca[tmpIndex] &&
          kp->dcaDaughters() < anaCuts::dcaDaughters[tmpIndex] &&
          kp->decayLength() > anaCuts::decayLength[tmpIndex] &&
          fabs(kp->lorentzVector().rapidity()) < anaCuts::RapidityCut &&
          ((kp->decayLength()) * sin(kp->pointingAngle())) < anaCuts::dcaV0ToPv[tmpIndex];
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const* const trk, float beta, StThreeVectorF const& vtx) const
{
   bool tofKaon = false;

   if (beta > 0)
   {
      double ptot = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
      float beta_k = ptot / sqrt(ptot * ptot + M_KAON_PLUS * M_KAON_PLUS);
      tofKaon = fabs(1 / beta - 1 / beta_k) < anaCuts::kTofBetaDiff ? true : false;
   }

   return tofKaon;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofPion(StPicoTrack const* const trk, float beta, StThreeVectorF const& vtx) const
{
   bool tofPion = false;

   if (beta > 0)
   {
      double ptot = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
      float beta_pi = ptot / sqrt(ptot * ptot + M_PION_PLUS * M_PION_PLUS);
      tofPion = fabs(1 / beta - 1 / beta_pi) < anaCuts::pTofBetaDiff ? true : false;
   }

   return tofPion;
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
   int index2tof = trk->bTofPidTraitsIndex();

   float beta = std::numeric_limits<float>::quiet_NaN();

   if (index2tof >= 0)
   {
      StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

      if (tofPid)
      {
         beta = tofPid->btofBeta();

         if (beta < 1e-4)
         {
            StThreeVectorF const btofHitPos = tofPid->btofHitPos();

            StPhysicalHelixD helix = trk->helix();
            float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
         }
      }
   }

   return beta;
}
//-----------------------------------------------------------------------------
int StPicoD0AnaMaker::trkHalf(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
   StThreeVectorF mom = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField());
   if (mom.phi() > anaCuts:: rightHalfLowEdge && mom.phi() < anaCuts:: rightHalfHighEdge) return +1; //right side
   else return -1;//lest side

}
