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
#include "StRoot/StRefMultCorr/StRefMultCorr.h"

#include "StMemStat.h"

ClassImp(StPicoD0AnaMaker)

StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name, TString const inputFilesList,
                                   TString const outFileBaseName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil):
   StMaker(name), mPicoDstMaker(picoDstMaker), mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil),
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
   mGRefMultCorrUtil = new StRefMultCorr("grefmult");

   return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
   delete mGRefMultCorrUtil;
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

   if (mPicoD0Event->runId() != picoDst->event()->runId() ||
         mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
      LOG_ERROR << " StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!" << endm;
      LOG_ERROR << " StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync." << endm;
      exit(1);
   }

   // -------------- USER ANALYSIS -------------------------

   mHists->addEventBeforeCut(picoDst->event());
   if (isGoodEvent(picoDst->event()))
   {
      TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();
      if (aKaonPion->GetEntries()) mHists->addEvent(picoDst->event());

      StThreeVectorF const pVtx = picoDst->event()->primaryVertex();
      StThreeVectorF const pVtxErr = picoDst->event()->primaryVertexError();


      if (!mGRefMultCorrUtil)
      {
         LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
         return kStWarn;
      }
      mGRefMultCorrUtil->init(picoDst->event()->runId());
      mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(), pVtx.z(), picoDst->event()->ZDCx()) ;

      int centrality  = mGRefMultCorrUtil->getCentralityBin9();
      const double reweight = mGRefMultCorrUtil->getWeight();
      const double refmultCor = mGRefMultCorrUtil->getRefMultCorr();
      mHists->addCent(refmultCor, centrality, reweight, pVtx.z());

      //Basiclly add some QA plots
      UInt_t nTracks = picoDst->numberOfTracks();
      //  cout<<"mem: begin: nTracks="<<nTracks<<" ,"<<mem.Used()<<" prog size "<<mem.ProgSize()<<endl;
      for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
      {
         StPicoTrack const* trk = picoDst->track(iTrack);
         if (!trk) continue;
         StPhysicalHelixD helix = trk->helix();
         float dca = float(helix.geometricSignedDistance(pVtx));
         StThreeVectorF momentum = trk->gMom(pVtx, picoDst->event()->bField());

         bool tofMatch = getTofBeta(trk, &pVtx) > 0;

         int tofMatchFlag =  tofMatch ? 1 : 0 ;
         int hftMatchFlag =  trk->isHFTTrack() ? 1 : 0 ;

         if (!isGoodQaTrack(trk, momentum, dca)) continue;

         StThreeVectorF dcaPoint = helix.at(helix.pathLength(pVtx.x(), pVtx.y()));
         float dcaZ = dcaPoint.z() - pVtx.z();
         StThreeVectorF dcaP = helix.momentumAt(pVtx.x(), pVtx.y());
         float dcaXy = ((dcaPoint - pVtx).x() * dcaP.y() - (dcaPoint - pVtx).y() * dcaP.x()) / dcaP.perp();

         // mHists->addQaNtuple(picoDst->event()->runId(), dca, pVtx.z(), momentum.perp(), momentum.pseudoRapidity(), momentum.phi(), centrality, refmultCor, picoDst->event()->ZDCx(), tofMatchFlag, hftMatchFlag);

         bool isPion = kFALSE;
         bool isKaon = kFALSE;
         if (fabs(trk->nSigmaPion()) < anaCuts::nSigmaPion)  isPion = kTRUE;
         if (fabs(trk->nSigmaKaon()) < anaCuts::nSigmaKaon)  isKaon = kTRUE;
         if (trk && tofMatch && fabs(dca) < 1.0 && trk->isHFTTrack() && (isPion || isKaon))
         {
            mHists->addDcaPtCent(dca, dcaXy, dcaZ, isPion, isKaon, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //add Dca distribution
         }
         if (trk && tofMatch && fabs(dca) < 1.5 && (isPion || isKaon))
         {
            mHists->addTpcDenom1(isPion, isKaon, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //Dca cut on 1.5cm, add Tpc Denominator
         }
         if (trk && tofMatch && fabs(dca) < 1.5 && trk->isHFTTrack() && (isPion || isKaon))
         {
            mHists->addHFTNumer1(isPion, isKaon, momentum.perp(), centrality,  momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //Dca cut on 1.5cm, add HFT Numerator
         }
      } // .. end tracks loop

      for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
      {
         StKaonPion const* kp = (StKaonPion*)aKaonPion->UncheckedAt(idx);

         if (!isGoodPair(kp)) continue;

         StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
         StPicoTrack const* pion = picoDst->track(kp->pionIdx());

         if (!isGoodTrack(kaon, kp) || !isGoodTrack(pion, kp)) continue;

         bool tpcPion = isTpcPion(pion);
         bool tpcKaon = isTpcKaon(kaon);
         float pBeta = getTofBeta(pion, &pVtx);
         float kBeta = getTofBeta(kaon, &pVtx);
         bool pTofAvailable = pBeta > 0;
         bool kTofAvailable = kBeta > 0;
         bool tofPion = isTofPion(pion, pBeta);
         bool tofKaon = isTofKaon(kaon, kBeta);

         bool goodPion = (pTofAvailable && tofPion) || (!pTofAvailable && tpcPion);
         bool goodKaon = (kTofAvailable && tofKaon) || (!kTofAvailable && tpcKaon);
         bool tof = goodPion && goodKaon;
         bool tpc = tpcPion && tpcKaon;

         if (tpc || tof)
         {
            bool unlike = kaon->charge() * pion->charge() < 0 ? true : false;
            mHists->addKaonPion(kp, unlike, tpc, tof, centrality, reweight);
         }

      } // end of kaonPion loop
   } // end of isGoodEvent
   //      cout<<"Lomnitz: end******* event:"<<mem.Used()<<" prog size "<<mem.ProgSize()<<endl;

   return kStOK;
}
//-----------------------------------------------------------------------------
int StPicoD0AnaMaker::getD0PtIndex(StKaonPion const * kp) const
{
   for (int i = 0; i < anaCuts::nPtBins; i++)
   {
      if ((kp->pt() >= anaCuts::PtBinsEdge[i]) && (kp->pt() < anaCuts::PtBinsEdge[i + 1]))
         return i;
   }
   return anaCuts::nPtBins - 1;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodEvent(StPicoEvent const * const picoEvent) const
{
   return (picoEvent->triggerWord() & anaCuts::triggerWord) &&
          fabs(picoEvent->primaryVertex().z()) < anaCuts::vz &&
          fabs(picoEvent->primaryVertex().z() - picoEvent->vzVpd()) < anaCuts::vzVpdVz &&
          !(fabs(picoEvent->primaryVertex().x()) < anaCuts::Verror && fabs(picoEvent->primaryVertex().y()) < anaCuts::Verror && fabs(picoEvent->primaryVertex().z()) < anaCuts::Verror) &&
          sqrt(TMath::Power(picoEvent->primaryVertex().x(), 2) + TMath::Power(picoEvent->primaryVertex().y(), 2)) <=  anaCuts::Vrcut;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodQaTrack(StPicoTrack const * const trk, StThreeVectorF const momentum, const double dca) const
{
   return trk->gPt() > anaCuts::qaGPt && trk->nHitsFit() >= anaCuts::qaNHitsFit;
   // trk->nHitsDedx() >= anaCuts::qaNHitsDedx && fabs(dca) < anaCuts::qaDca &&
   // momentum.pseudoRapidity() < anaCuts::qaEta;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk, StKaonPion const * kp) const
{
   StThreeVectorF mom = trk->gMom(mPicoDstMaker->picoDst()->event()->primaryVertex(), mPicoDstMaker->picoDst()->event()->bField());
   return trk->gPt() > anaCuts::minPt &&
          trk->nHitsFit() >= anaCuts::nHitsFit &&
          fabs(mom.pseudoRapidity()) <= anaCuts::Eta;
   //    1;
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
   int tmpIndex = getD0PtIndex(kp);
   return cos(kp->pointingAngle()) > anaCuts::cosTheta[tmpIndex] &&
          kp->pionDca() > anaCuts::pDca[tmpIndex] && kp->kaonDca() > anaCuts::kDca[tmpIndex] &&
          kp->dcaDaughters() < anaCuts::dcaDaughters[tmpIndex] &&
          kp->decayLength() > anaCuts::decayLength[tmpIndex] &&
          fabs(kp->lorentzVector().rapidity()) < anaCuts::RapidityCut &&
          ((kp->decayLength()) * sin(kp->pointingAngle())) < anaCuts::dcaV0ToPv[tmpIndex];
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
   bool tofKaon = false;

   if (beta > 0)
   {
      double ptot = trk->dcaGeometry().momentum().mag();
      float beta_k = ptot / sqrt(ptot * ptot + M_KAON_PLUS * M_KAON_PLUS);
      tofKaon = fabs(1 / beta - 1 / beta_k) < anaCuts::kTofBetaDiff ? true : false;
   }

   return tofKaon;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofPion(StPicoTrack const * const trk, float beta) const
{
   bool tofPion = false;

   if (beta > 0)
   {
      double ptot = trk->dcaGeometry().momentum().mag();
      float beta_pi = ptot / sqrt(ptot * ptot + M_PION_PLUS * M_PION_PLUS);
      tofPion = fabs(1 / beta - 1 / beta_pi) < anaCuts::pTofBetaDiff ? true : false;
   }

   return tofPion;
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
   int index2tof = trk->bTofPidTraitsIndex();

   float beta = std::numeric_limits<float>::quiet_NaN();

   if (index2tof >= 0)
   {
      StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

      if (tofPid)
      {
         beta = tofPid->btofBeta();

         if (beta < 1e-4)
         {
            StThreeVectorF const btofHitPos = tofPid->btofHitPos();

            StPhysicalHelixD helix = trk->helix();
            float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
         }
      }
   }

   return beta;
}
