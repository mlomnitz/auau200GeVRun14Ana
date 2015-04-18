#include <cmath>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TString.h"
#include "../StPicoDstMaker/StPicoEvent.h"
#include "../StPicoPrescales/StPicoPrescales.h"
#include "StKaonPion.h"
#include "StAnaCuts.h"

#include "StPicoD0AnaHists.h"

ClassImp(StPicoD0AnaHists)

//-----------------------------------------------------------------------
StPicoD0AnaHists::StPicoD0AnaHists(TString fileBaseName) : mPrescales(NULL), mOutFile, 
  mh2InvariantMassVsPt(NULL), mh2InvariantMassVsPtLike(NULL), mh2InvariantMassVsPtTof(NULL), mh2InvariantMassVsPtTofLike(NULL)
{
  mPrescales = new StPicoPrescales(cuts::prescalesFilesDirectoryName);

  mOutFile = new TFile(Form("%s.picoD0.hists.root",fileBaseName.Data()),"RECREATE");

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
void StPicoD0AnaHists::addEvent(StPicoEvent const& picoEvent)
{
  int runIndex = mPrescales->runIndex(picoD0Event.runId());
  mh1TotalEventsInRun->Fill(runIndex);
}
//---------------------------------------------------------------------
void StPicoD0AnaHists::addKaonPion(StKaonPion const* const kp, bool unlike, bool tof)
{
  if(unlike)
  {
    mh2InvariantMassVsPt->Fill(kp->pt(),kp->m());
    if(tof) mh2InvariantMassVsPtTof->Fill(kp->pt(),kp->m());
  }
  else
  {
    mh2InvariantMassVsPtLike->Fill(kp->pt(),kp->m());
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
