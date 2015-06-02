#include <cmath>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TString.h"
#include "../StPicoDstMaker/StPicoEvent.h"
#include "../StPicoPrescales/StPicoPrescales.h"
#include "../StPicoD0EventMaker/StKaonPion.h"
#include "StAnaCuts.h"

#include "StPicoD0AnaHists.h"

ClassImp(StPicoD0AnaHists)

//-----------------------------------------------------------------------
StPicoD0AnaHists::StPicoD0AnaHists(TString fileBaseName) : mPrescales(NULL), mOutFile(NULL), 
  mh2InvariantMassVsPt(NULL), mh2InvariantMassVsPtLike(NULL), mh2InvariantMassVsPtTof(NULL), mh2InvariantMassVsPtTofLike(NULL),
  mh1Cent(NULL), mh1CentWg(NULL), mh1gRefmultCor(NULL), mh1gRefmultCorWg(NULL), mh3InvariantMassVsPtVsCent(NULL), mh3InvariantMassVsPtVsCentLike(NULL), mh3InvariantMassVsPtVsCentTof(NULL), mh3InvariantMassVsPtVsCentTofLike(NULL),
   mh2Tpc1PtCent(NULL), mh2Tpc2PtCent(NULL), mh2HFT1PtCent(NULL), mh2HFT2PtCent(NULL)
{
  mPrescales = new StPicoPrescales(anaCuts::prescalesFilesDirectoryName);

  mOutFile = new TFile(Form("%s.hists.root",fileBaseName.Data()),"RECREATE");

  int nRuns = mPrescales->numberOfRuns();
  TH1::SetDefaultSumw2();
  mh1TotalEventsInRun         = new TH1F("mh1TotalEventsInRun","totalEventsInRun;runIndex;totalEventsInRun",nRuns+1,0,nRuns+1);
  mh1TotalEventsInRunBeforeCut= new TH1F("mh1TotalEventsInRunBeforeCut","totalEventsInRun;runIndex;totalEventsInRun",nRuns+1,0,nRuns+1);
  mh2InvariantMassVsPt        = new TH2F("mh2InvariantMassVsPt","invariantMassVsPt;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})",120,0,12,50,1.6,2.1);
  mh2InvariantMassVsPtLike    = new TH2F("mh2InvariantMassVsPtLike","invariantMassVsPtLike;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})",120,0,12,50,1.6,2.1);
  mh2InvariantMassVsPtTof     = new TH2F("mh2InvariantMassVsPtTof","invariantMassVsPtTof;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})",120,0,12,50,1.6,2.1);
  mh2InvariantMassVsPtTofLike = new TH2F("mh2InvariantMassVsPtTofLike","invariantMassVsPtTofLike;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})",120,0,12,50,1.6,2.1);
  //add centrality
  mh1Cent         = new TH1F("mh1Cent","EventsVsCentrality;cent;Counts",10,-1.5,8.5);
  mh1CentWg         = new TH1F("mh1CentWg","EventsVsCentrality;cent;Counts",10,-1.5,8.5);
  mh1gRefmultCor  = new TH1F("mh1gRefmultCor","gRefmultCor;gRefmult;Counts",700,0,700);
  mh1gRefmultCorWg  = new TH1F("mh1gRefmultCorWg","gRefmultCorWg;gRefmultCorWg;Counts",700,0,700);
  mh3InvariantMassVsPtVsCent        = new TH3F("mh3InvariantMassVsPtVsCent","invariantMassVsPtVsCent;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})",120,0,12,10,-1.5,8.5,50,1.6,2.1);
  mh3InvariantMassVsPtVsCentLike    = new TH3F("mh3InvariantMassVsPtVsCentLike","invariantMassVsPtVsCentLike;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})",120,0,12,10,-1.5,8.5,50,1.6,2.1);
  mh3InvariantMassVsPtVsCentTof     = new TH3F("mh3InvariantMassVsPtVsCentTof","invariantMassVsPtVsCentTof;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})",120,0,12,10,-1.5,8.5,50,1.6,2.1);
  mh3InvariantMassVsPtVsCentTofLike = new TH3F("mh3InvariantMassVsPtVsCentTofLike","invariantMassVsPtVsCentTofLike;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})",120,0,12,10,-1.5,8.5,50,1.6,2.1);
  //Add some HFT ratio plots
  mh2Tpc1PtCent  = new TH2F("mh2Tpc1PtCent","Tpc tacks;p_{T}(GeV/c);cent",120,0,12,10,-1.5,8.5);//Dca 1.5cm
  mh2Tpc2PtCent  = new TH2F("mh2Tpc2PtCent","Tpc tacks;p_{T}(GeV/c);cent",120,0,12,10,-1.5,8.5);//Dca 0.1cm
  mh2HFT1PtCent  = new TH2F("mh2HFT1PtCent","HFT tacks;p_{T}(GeV/c);cent",120,0,12,10,-1.5,8.5);//Dca 1.5cm
  mh2HFT2PtCent  = new TH2F("mh2HFT2PtCent","HFT tacks;p_{T}(GeV/c);cent",120,0,12,10,-1.5,8.5);//Dca 1.5cm
}
StPicoD0AnaHists::~StPicoD0AnaHists()
{
  delete mPrescales;
  // note that histograms are owned by mOutFile. They will be destructed 
  // when the file is closed.
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addEvent(StPicoEvent const* const picoEvent)
{
  int runIndex = mPrescales->runIndex(picoEvent->runId());
  mh1TotalEventsInRun->Fill(runIndex);
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addEventBeforeCut(StPicoEvent const* const picoEvent)
{
  int runIndex = mPrescales->runIndex(picoEvent->runId());
  mh1TotalEventsInRunBeforeCut->Fill(runIndex);
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addCent(const double refmultCor,int centrality, const double reweight)
{
  mh1gRefmultCor->Fill(refmultCor);
  mh1gRefmultCorWg->Fill(refmultCor,reweight);
  mh1Cent->Fill(centrality);
  mh1CentWg->Fill(centrality,reweight);
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addTpcDenom1(double pt, int centrality)
{
  mh2Tpc1PtCent->Fill(pt,centrality);
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addTpcDenom2(double pt, int centrality)
{
  mh2Tpc2PtCent->Fill(pt,centrality);
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addHFTNumer1(double pt, int centrality)
{
  mh2HFT1PtCent->Fill(pt,centrality);
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addHFTNumer2(double pt, int centrality)
{
  mh2HFT2PtCent->Fill(pt,centrality);
}
//---------------------------------------------------------------------
void StPicoD0AnaHists::addKaonPion(StKaonPion const* const kp, bool unlike, bool tpc, bool tof, int centrality, const double reweight)
{
  if(unlike)
  {
    if(tpc) mh2InvariantMassVsPt->Fill(kp->pt(),kp->m());
    if(tof) mh2InvariantMassVsPtTof->Fill(kp->pt(),kp->m());
    if(tpc) mh3InvariantMassVsPtVsCent->Fill(kp->pt(),centrality,kp->m(),reweight);
    if(tof) mh3InvariantMassVsPtVsCentTof->Fill(kp->pt(),centrality,kp->m(),reweight);
  }
  else
  {
    if(tpc) mh2InvariantMassVsPtLike->Fill(kp->pt(),kp->m());
    if(tof) mh2InvariantMassVsPtTofLike->Fill(kp->pt(),kp->m());
    if(tpc) mh3InvariantMassVsPtVsCentLike->Fill(kp->pt(),centrality,kp->m(),reweight);
    if(tof) mh3InvariantMassVsPtVsCentTofLike->Fill(kp->pt(),centrality,kp->m(),reweight);
  }
}
//---------------------------------------------------------------------
void StPicoD0AnaHists::closeFile()
{
  mOutFile->cd();
  mh1TotalEventsInRun->Write();
  mh1TotalEventsInRunBeforeCut->Write();
  mh2InvariantMassVsPt->Write();
  mh2InvariantMassVsPtLike->Write();
  mh2InvariantMassVsPtTof->Write();
  mh2InvariantMassVsPtTofLike->Write();
  //centrality
  mh1Cent->Write();
  mh1CentWg->Write();
  mh1gRefmultCor->Write();
  mh1gRefmultCorWg->Write();
  mh3InvariantMassVsPtVsCent->Write();
  mh3InvariantMassVsPtVsCentLike->Write();
  mh3InvariantMassVsPtVsCentTof->Write();
  mh3InvariantMassVsPtVsCentTofLike->Write();
  //HFT ratio QA
  mh2Tpc1PtCent->Write();
  mh2Tpc2PtCent->Write();
  mh2HFT1PtCent->Write();
  mh2HFT2PtCent->Write();

  mOutFile->Close();
}
