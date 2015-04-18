#include <cmath>
#include "TH1F.h"
#include "TH2F.h"
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
  mh2InvariantMassVsPt(NULL), mh2InvariantMassVsPtLike(NULL), mh2InvariantMassVsPtTof(NULL), mh2InvariantMassVsPtTofLike(NULL)
{
  mPrescales = new StPicoPrescales(anaCuts::prescalesFilesDirectoryName);

  mOutFile = new TFile(Form("%s.hists.root",fileBaseName.Data()),"RECREATE");

  int nRuns = mPrescales->numberOfRuns();
  TH1::SetDefaultSumw2();
  mh1TotalEventsInRun         = new TH1F("mh1TotalEventsInRun","totalEventsInRun;runIndex;totalEventsInRun",nRuns+1,0,nRuns+1);
  mh2InvariantMassVsPt        = new TH2F("mh2InvariantMassVsPt","invariantMassVsPt;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})",120,0,12,50,1.6,2.1);
  mh2InvariantMassVsPtLike    = new TH2F("mh2InvariantMassVsPtLike","invariantMassVsPtLike;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})",120,0,12,50,1.6,2.1);
  mh2InvariantMassVsPtTof     = new TH2F("mh2InvariantMassVsPtTof","invariantMassVsPtTof;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})",120,0,12,50,1.6,2.1);
  mh2InvariantMassVsPtTofLike = new TH2F("mh2InvariantMassVsPtTofLike","invariantMassVsPtTofLike;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})",120,0,12,50,1.6,2.1);
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
  int runIndex = mPrescales->runIndex(picoEvent.runId());
  mh1TotalEventsInRun->Fill(runIndex);
}
//---------------------------------------------------------------------
void StPicoD0AnaHists::addKaonPion(StKaonPion const* const kp, bool unlike, bool tpc, bool tof)
{
  if(unlike)
  {
    if(tpc) mh2InvariantMassVsPt->Fill(kp->pt(),kp->m());
    if(tof) mh2InvariantMassVsPtTof->Fill(kp->pt(),kp->m());
  }
  else
  {
    if(tpc) mh2InvariantMassVsPtLike->Fill(kp->pt(),kp->m());
    if(tof) mh2InvariantMassVsPtTofLike->Fill(kp->pt(),kp->m());
  }
}
//---------------------------------------------------------------------
void StPicoD0AnaHists::closeFile()
{
  mOutFile->cd();
  mh1TotalEventsInRun->Write();
  mh2InvariantMassVsPt->Write();
  mh2InvariantMassVsPtLike->Write();
  mh2InvariantMassVsPtTof->Write();
  mh2InvariantMassVsPtTofLike->Write();
  mOutFile->Close();
}
