#ifndef StDpmHists__h
#define StDpmHists__h

/* **************************************************
 *  A class to create and save my Dpm analysis histograms.
 *
 *  Authors: Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */
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
class StPicoKPiX;
class StEventPlane;

class StDpmHists
{
public:
  StDpmHists(std::string fileBaseName, int harmonic = 2);
   virtual ~StDpmHists();
   void closeFile();
   void addEvent(int const, float const, float const, float const);
   void addKaPiX(bool unlike, double const, StPicoKPiX const*, float const, float const, float const, StEventPlane*, const int, const double);
   int getDpmPtIndex(StPicoKPiX const* triplet, std::vector<float> const& edges) const;
   //
 private:
   TH3F* hCentVzPsi;
   TH3F* hCentVzPsiNoWeight;
   // Half eta's event plane 
   THn* hDpmEtaSubCentPtMDphi;
   THn* hDpmEtaSubCentPtMDphiLikeSign;
   // Eta Gap Event Plane
   THn* hDpmCentPtMDphiEtaGap;
   THn* hDpmCentPtMDphiEtaGapLikeSign;
   int mHarmonic;
   TFile *mOutfile;
};

#endif
