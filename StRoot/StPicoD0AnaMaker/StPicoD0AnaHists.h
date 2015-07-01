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
#include "StAnaCuts.h"
#include "THnSparse.h" 
#include "THn.h" 

class TH1F;
class TH2F;
class TH3F;
class TFile;
class TString;
class StPicoPrescales;
class StPicoEvent;
class StKaonPion;
class TNtuple;


class StPicoD0AnaHists: public TObject
{
  public:
   StPicoD0AnaHists(TString fileBaseName);
   virtual ~StPicoD0AnaHists();
   void addEvent(StPicoEvent const *);
   void addEventBeforeCut(StPicoEvent const *);
   void addCent(const double refmultCor,int centrality, const double reweight, const float vz);
   void addKaonPion(StKaonPion const*, bool unlike, bool tpc, bool tof, int centrality, const double reweight);
   void addTpcDenom1(bool IsPion, bool IsKaon, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX);
   void addHFTNumer1(bool IsPion, bool IsKaon, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX);
   void addQaNtuple(int, float, float, float, float, float, int, const double, float, int, int);
   void addDcaPtCent(float dca, float dcaXy, float  dcaZ, bool IsPion, bool IsKaon, float pt,  int centrality, float Eta, float Phi, float Vz, float ZdcX);
   int getEtaIndex(float Eta) ;
   int getPhiIndex(float Phi) ;
   int getVzIndex(float Vz) ;
   int getZdcxIndex(float ZdcX) ;
   void closeFile();

  private:
   StPicoD0AnaHists(){}

   StPicoPrescales* mPrescales;
   TFile* mOutFile;
   TH1F* mh1TotalEventsInRun;
   TH1F* mh1TotalEventsInRunBeforeCut;
   TH2F* mh2InvariantMassVsPt;
   TH2F* mh2InvariantMassVsPtLike;
   TH2F* mh2InvariantMassVsPtTof;
   TH2F* mh2InvariantMassVsPtTofLike;
   //centrality
   TH1F* mh1Cent;
   TH1F* mh1CentWg;
   TH1F* mh1gRefmultCor;
   TH1F* mh1gRefmultCorWg;
   TH2F* mh2CentVz;
   TH2F* mh2CentVzWg;
   TH3F* mh3InvariantMassVsPtVsCent;
   TH3F* mh3InvariantMassVsPtVsCentLike;
   TH3F* mh3InvariantMassVsPtVsCentTof;
   TH3F* mh3InvariantMassVsPtVsCentTofLike;
   //HFT ratio QA
   TH2F* mh2Tpc1PtCent;
   TH2F* mh2Tpc1PhiVz;
   TH2F* mh2HFT1PhiVz;
   TH2F* mh2HFT1PtCent;
   TH2F* mh2Tpc1PtCentPartEtaVz[anaCuts::nParticles][anaCuts::nEtas][anaCuts::nVzs];
   TH2F* mh2Tpc1PtCentPartPhi[anaCuts::nParticles][anaCuts::nPhis];
   TH2F* mh2Tpc1PtCentPartZdcx[anaCuts::nParticles][anaCuts::nZdcxs];
   TH2F* mh2HFT1PtCentPartEtaVz[anaCuts::nParticles][anaCuts::nEtas][anaCuts::nVzs];
   TH2F* mh2HFT1PtCentPartPhi[anaCuts::nParticles][anaCuts::nPhis];
   TH2F* mh2HFT1PtCentPartZdcx[anaCuts::nParticles][anaCuts::nZdcxs];

   //HFT Dca 
   TH3F* mh3DcaXyPtCentPartEtaVz[anaCuts::nParticles][anaCuts::nEtas][anaCuts::nVzs];
   TH3F* mh3DcaXyPtCentPartPhi[anaCuts::nParticles][anaCuts::nPhis];
   TH3F* mh3DcaXyPtCentPartZdcx[anaCuts::nParticles][anaCuts::nZdcxs];
   TH3F* mh3DcaZPtCentPartEtaVz[anaCuts::nParticles][anaCuts::nEtas][anaCuts::nVzs];
   TH3F* mh3DcaZPtCentPartPhi[anaCuts::nParticles][anaCuts::nPhis];
   TH3F* mh3DcaZPtCentPartZdcx[anaCuts::nParticles][anaCuts::nZdcxs];
   TH3F* mh3DcaXyPtCent;
   TH3F* mh3DcaZPtCent;

   //TNtuple* nt;

   ClassDef(StPicoD0AnaHists, 1)
};
#endif
