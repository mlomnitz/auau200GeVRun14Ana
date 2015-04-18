#ifndef StPicoD0AnaHists__h
#define StPicoD0AnaHists__h

/* **************************************************
 *  A class to create and save my D0 analysis histograms.
 *
 *  Authors: Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

#include "TObject.h"

class TH1F;
class TH2F;
class TFile;
class TString;
class StPicoPrescales;
class StPicoEvent;
class StKaonPion;


class StPicoD0AnaHists: public TObject
{
  public:
   StPicoD0AnaHists(TString fileBaseName);
   virtual ~StPicoD0AnaHists();
   void addEvent(StPicoEvent const &);
   void addKaonPion(StKaonPion const*, bool unlike, bool tof);
   void closeFile();

  private:
   StPicoD0AnaHists(){}

   StPicoPrescales* mPrescales;
   TFile* mOutFile;
   TH1F* mh1TotalEventsInRun;
   TH2F* mh2InvariantMassVsPt;
   TH2F* mh2InvariantMassVsPtLike;
   TH2F* mh2InvariantMassVsPtTof;
   TH2F* mh2InvariantMassVsPtTofLike;

   ClassDef(StPicoD0AnaHists, 1)
};
#endif
