/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoKPiEvent
 *  simultaneously and do analysis.
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoCharmContainers/StPicoKPiXEvent.h"
#include "StPicoCharmContainers/StPicoKPiX.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"

#include "StKPiXAnaCuts.h"
#include "StPicoKPiXAnaMaker.h"

ClassImp(StPicoKPiXAnaMaker)

StPicoKPiXAnaMaker::StPicoKPiXAnaMaker(char const * name, TString const inputFilesList,
                                       std::string outFileBaseName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil):
                                       StMaker(name), mPicoDstMaker(picoDstMaker), mPicoKPiXEvent(nullptr),
                                       mGRefMultCorrUtil(grefmultCorrUtil), mChain(nullptr),
                                       mInputFilesList(inputFilesList), mOutFileBaseName(outFileBaseName),
                                       mEventCounter(0)
{}

Int_t StPicoKPiXAnaMaker::Init()
{
   mPicoKPiXEvent = new StPicoKPiXEvent();

   mChain = new TChain("KPiXTree");
   std::ifstream listOfFiles(mInputFilesList);
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoKPiXAnaMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoKPiXAnaMaker - Could not open list of files. ABORT!" << endm;
      return kStErr;
   }

   mChain->GetBranch("kPiXEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("kPiXEvent", &mPicoKPiXEvent);

   // mHists = new StPicoD0AnaHists(mOutFileBaseName,mFillQaHists);

   return kStOK;
}

StPicoKPiXAnaMaker::~StPicoKPiXAnaMaker()
{
}

Int_t StPicoKPiXAnaMaker::Finish()
{
   // mHists->closeFile();
   return kStOK;
}

Int_t StPicoKPiXAnaMaker::Make()
{
   readNextEvent();

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoKPiXAnaMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();

   if (!picoDst)
   {
      LOG_WARN << "StPicoKPiXAnaMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if (mPicoKPiXEvent->runId() != picoDst->event()->runId() ||
         mPicoKPiXEvent->eventId() != picoDst->event()->eventId())
   {
      LOG_ERROR << " StPicoKPiXAnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!" << "\n";
      LOG_ERROR << " StPicoKPiXAnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoKPiXEvent are not in sync." << endm;
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


   if(isGoodTrigger(picoDst->event()))
   {
     // mHists->addEventBeforeCut(picoDst->event());

     StThreeVectorF pVtx = picoDst->event()->primaryVertex();
     if (isGoodEvent(picoDst->event(),pVtx))
     {
       TClonesArray const* aKaonPionXaon = mPicoKPiXEvent->kaonPionXaonArray();
       // mHists->addEvent(picoDst->event());

       mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(), pVtx.z(), picoDst->event()->ZDCx()) ;
       int centrality  = mGRefMultCorrUtil->getCentralityBin9();
       const double reweight = mGRefMultCorrUtil->getWeight();
       const double refmultCor = mGRefMultCorrUtil->getRefMultCorr();

       // mHists->addCent(refmultCor, centrality, reweight, pVtx.z());

       for (int idx = 0; idx < aKaonPionXaon->GetEntries(); ++idx)
       {
         StPicoKPiX const* const kpx = (StPicoKPiX*)aKaonPionXaon->UncheckedAt(idx);
         if(!kpx) continue;

         if (!kPiXAnaCuts::isGoodKPiX(*kpx, kPiXAnaCuts::DpmCuts)) continue;

         StPicoTrack const* const kaon = picoDst->track(kpx->kaonIdx());
         StPicoTrack const* const pion = picoDst->track(kpx->pionIdx());
         StPicoTrack const* const xaon = picoDst->track(kpx->xaonIdx());

         if (!isGoodTrack(kaon, pVtx) ||
             !isGoodTrack(pion, pVtx) ||
             !isGoodTrack(xaon, pVtx)) continue;

         // PID
         if(!isTpcPion(pion) || !isTpcKaon(kaon)) continue;
         float pBeta = getTofBeta(pion, pVtx);
         float kBeta = getTofBeta(kaon, pVtx);
         bool pTofAvailable = !isnan(pBeta) && pBeta > 0;
         bool kTofAvailable = !isnan(kBeta) && kBeta > 0;
         bool tofPion = pTofAvailable ? isTofPion(pion, pBeta, pVtx) : true;//this is bybrid pid, not always require tof
         bool tofKaon = kTofAvailable ? isTofKaon(kaon, kBeta, pVtx) : true;//this is bybrid pid, not always require tof
         bool tof = tofPion && tofKaon;

         bool unlike = kaon->charge() * pion->charge() < 0 ? true : false;

         // mHists->addKaonPion(kp, unlike, true, tof, centrality, reweight);

       } // kaonPionXaon loop
     } // isGoodEvent
   } // isGoodTrigger

   return kStOK;
}

bool StPicoKPiXAnaMaker::isGoodEvent(StPicoEvent const* const picoEvent, StThreeVectorF const& pVtx) const
{
   return fabs(pVtx.z()) < kPiXAnaCuts::vz &&
          fabs(pVtx.z() - picoEvent->vzVpd()) < kPiXAnaCuts::vzVpdVz &&
          !(fabs(pVtx.x()) < kPiXAnaCuts::Verror && fabs(pVtx.y()) < kPiXAnaCuts::Verror && fabs(pVtx.z()) < kPiXAnaCuts::Verror) &&
          sqrt(pow(pVtx.x(), 2) + pow(pVtx.y(), 2)) <=  kPiXAnaCuts::Vrcut;
}

bool StPicoKPiXAnaMaker::isGoodTrigger(StPicoEvent const* const picoEvent) const
{
  for(auto trg: kPiXAnaCuts::triggers)
  {
    if(picoEvent->isTrigger(trg)) return true;
  }

  return false;
}

bool StPicoKPiXAnaMaker::isGoodTrack(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
   StThreeVectorF mom = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField());


   return mom.perp() > kPiXAnaCuts::minPt &&
          trk->nHitsFit() >= kPiXAnaCuts::nHitsFit &&
          fabs(mom.pseudoRapidity()) <= kPiXAnaCuts::Eta;
}

bool StPicoKPiXAnaMaker::isTpcPion(StPicoTrack const* const trk) const
{
   return fabs(trk->nSigmaPion()) < kPiXAnaCuts::nSigmaPion;
}

bool StPicoKPiXAnaMaker::isTpcKaon(StPicoTrack const* const trk) const
{
   return fabs(trk->nSigmaKaon()) < kPiXAnaCuts::nSigmaKaon;
}

bool StPicoKPiXAnaMaker::isTofKaon(StPicoTrack const* const trk, float beta, StThreeVectorF const& vtx) const
{
   bool tofKaon = false;

   if (beta > 0)
   {
      double ptot = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
      float beta_k = ptot / sqrt(ptot * ptot + M_KAON_PLUS * M_KAON_PLUS);
      tofKaon = fabs(1 / beta - 1 / beta_k) < kPiXAnaCuts::kTofBetaDiff ? true : false;
   }

   return tofKaon;
}

bool StPicoKPiXAnaMaker::isTofPion(StPicoTrack const* const trk, float beta, StThreeVectorF const& vtx) const
{
   bool tofPion = false;

   if (beta > 0)
   {
      double ptot = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
      float beta_pi = ptot / sqrt(ptot * ptot + M_PION_PLUS * M_PION_PLUS);
      tofPion = fabs(1 / beta - 1 / beta_pi) < kPiXAnaCuts::pTofBetaDiff ? true : false;
   }

   return tofPion;
}

float StPicoKPiXAnaMaker::getTofBeta(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
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

int StPicoKPiXAnaMaker::getIndex(float const value, std::vector<float> const& edges) const
{
  for (size_t i = 0; i < edges.size(); ++i)
  {
    if (value >= edges[i] && value < edges[i + 1])
      return i;
  }

  return edges.size() - 1;
}

bool StPicoKPiXAnaMaker::isGoodKPiX(StPicoKPiX const* const kpx, kPiXAnaCuts::TopologicalCuts const& cuts) const
{
  StLorentzVectorF const fMom = kpx->fourMom(cuts.xMassHypothesis);

  int const tmpIndex = getIndex(fMom.perp(), cuts.ptBinsEdge);

  return cos(kpx->pointingAngle()) > cuts.cosTheta[tmpIndex] &&
         kpx->pionDca() > cuts.pDca[tmpIndex] &&
         kpx->kaonDca() > cuts.kDca[tmpIndex] &&
         kpx->xaonDca() > cuts.pDca[tmpIndex] &&
         kpx->dcaDaughters() < cuts.dcaDaughters[tmpIndex] &&
         kpx->decayLength() > cuts.decayLength[tmpIndex] &&
         fabs(fMom.rapidity()) < cuts.rapidityCut &&
         kpx->perpDcaToVtx() < cuts.dcaV0ToPv[tmpIndex];
}
