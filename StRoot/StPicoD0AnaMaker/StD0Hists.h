#ifndef StD0Hists__h
#define StD0Hists__h

/* **************************************************
 *  A class to create and save my D0 analysis histograms.
 *
 *  Authors: Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */
#include "TopologyCuts.h"
#include <string>
#include "TObject.h"

class TH1D;
class TH1F;
class TH2F;
class TH3F;
class THn;
class TProfile;
class TFile;
class TString;
class StKaonPion;

class StD0Hists
{
public:
  StD0Hists(std::string fileBaseName, int harmonic = 2);
   virtual ~StD0Hists();
   void closeFile();
   //void addKaonPion(StKaonPion const*, bool unlike, bool tpc, bool tof, int centrality, const double reweight);
   bool isGoodPair(StKaonPion const* pair, topoCuts::TopologicalCuts const& cuts) const;
   int getD0PtIndex(StKaonPion const* pair, std::vector<float> const& edges) const;
   //
   TH3F* hCentVzPsi;
   TH3F* hCentVzPsiSameEventNoWeight;
   TH3F* hCentVzPsiSameEvent;
   // Half eta's event plane 
   THn* hD0EtaSubCentPtMDphi[topoCuts::nCutsSets];
   THn* hD0EtaSubCentPtMDphiLikeSign[topoCuts::nCutsSets];
   // Eta Gap Event Plane
   THn* hD0CentPtMDphiEtaGap[topoCuts::nCutsSets];
   THn* hD0CentPtMDphiEtaGapLikeSign[topoCuts::nCutsSets];
};

#endif
