#ifndef StAnaCuts_H
#define StAnaCuts_H

/* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "Rtypes.h"
#include <string>
#include <array>
#include <vector>

// #include "phys_constants.h"

namespace kPiXAnaCuts
{
  // path to lists of triggers prescales
  // lists are obtained from http://www.star.bnl.gov/protected/common/common2014/trigger2014/plots_au200gev/
  std::string const prescalesFilesDirectoryName = "./run14AuAu200GeVPrescales";

  //trigger
  std::array<unsigned int, 5> const triggers = { 450050,    // vpdmb-5-p-nobsmd-hlt
                                                 450060,    // vpdmb-5-p-nobsmd-hlt
                                                 450005,    // vpdmb-5-p-nobsmd
                                                 450015,    // vpdmb-5-p-nobsmd
                                                 450025};    // vpdmb-5-p-nobsmd

  // event
  float const vz = 6.0;// cm.
  float const vzVpdVz = 3.0; // 3 cm.
  float const Verror = 1.0e-5; //
  float const Vrcut = 2.0; //

  // all tracks
  float const minPt = 0.5;
  int   const nHitsFit = 20;
  float const Eta = 1.0;

  //pions
  float const nSigmaPion = 3.0;
  float const pTofBetaDiff = 0.03;

  //kaons
  float const nSigmaKaon = 2.0;
  float const kTofBetaDiff = 0.03;

  struct TopologicalCuts
  {
    double xMassHypothesis;
    float rapidityCut;
    int   nPtBins;
    std::vector<float> ptBinsEdge;
    std::vector<float> dcaV0ToPv;
    std::vector<float> decayLength;
    std::vector<float> cosTheta;
    std::vector<float> dcaDaughters;
    std::vector<float> kDca;
    std::vector<float> pDca;

    // TopologicalCuts() = delete;
    TopologicalCuts(int NumberOfPtBins): nPtBins(NumberOfPtBins),
                                         ptBinsEdge(nPtBins+1),
                                         dcaV0ToPv(nPtBins),
                                         decayLength(nPtBins),
                                         cosTheta(nPtBins),
                                         dcaDaughters(nPtBins),
                                         kDca(nPtBins),
                                         pDca(nPtBins) {}

    static TopologicalCuts makeCuts(int NumberOfPtBins)
    {
      return TopologicalCuts(NumberOfPtBins);
    }
  };

  TopologicalCuts DpmCuts = TopologicalCuts::makeCuts(5);
  // DpmCuts.xMassHypothesis = M_PION_PLUS;
  // DpmCuts.rapidityCut = 1.0;
  // DpmCuts.ptBinsEdge   = {0., 1., 2., 3., 5., 15.};
  // DpmCuts.dcaV0ToPv    = {0.0061, 0.0049, 0.0038, 0.0038, 0.0040};
  // DpmCuts.decayLength  = {0.0145, 0.0181, 0.0212, 0.0247, 0.0259};
  // DpmCuts.cosTheta     = {0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
  // DpmCuts.dcaDaughters = {0.0084, 0.0066, 0.0057, 0.0050, 0.0060};
  // DpmCuts.kDca         = {0.0103, 0.0091, 0.0095, 0.0079, 0.0058};
  // DpmCuts.pDca         = {0.0110, 0.0111, 0.0086, 0.0081, 0.0062};

  // utility methods

  int getIndex(float const value, std::vector<float> const& edges)
  {
    for (size_t i = 0; i < edges.size(); ++i)
    {
      if (value >= edges[i] && value < edges[i + 1])
        return i;
    }

    return edges.size() - 1;
  }

  bool isGoodKPiX(StPicoKPiX const& kpx, TopologicalCuts const& cuts)
  {
    StLorentzVectorF const fMom = kpx.fourMom(cuts.xMassHypothesis);

    int const tmpIndex = getIndex(fMom.perp(), cuts.ptBinsEdge);

    return cos(kpx.pointingAngle()) > cuts.cosTheta[tmpIndex] &&
      kpx.pionDca() > cuts.pDca[tmpIndex] &&
      kpx.kaonDca() > cuts.kDca[tmpIndex] &&
      kpx.xaonDca() > cuts.pDca[tmpIndex] &&
      kpx.dcaDaughters() < cuts.dcaDaughters[tmpIndex] &&
      kpx.decayLength() > cuts.decayLength[tmpIndex] &&
      fabs(fMom.rapidity()) < cuts.rapidityCut &&
      kpx.perpDcaToVtx() < cuts.dcaV0ToPv[tmpIndex];
  }


}
#endif
