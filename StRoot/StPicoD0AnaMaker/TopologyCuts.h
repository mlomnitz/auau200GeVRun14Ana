#ifndef TopologyCuts__H
#define TopologyCuts__H

/* **************************************************
 *  Topology Cut Sets
 *
 *  Authors:  Michael Lomnitz (mrlomnitz@lbl.gov)
 *
 * **************************************************
 */
#include <string>

#include "TopologyCuts.h"
namespace topoCuts
{
  float const massMin = 0.;
  float const massMax = 2.5;
  int const nCutsSets = 3;
  std::string const cutsSetName[nCutsSets] = {"standard","50efficiency","150efficiency"};

  struct TopologicalCuts
  {
    float  RapidityCut;
    int    nPtBins;
    std::vector<float> PtEdge;
    std::vector<float> dcaV0ToPv;
    std::vector<float> decayLength;
    std::vector<float> cosTheta;
    std::vector<float> dcaDaughters;
    std::vector<float> kDca;
    std::vector<float> pDca;
    
  TopologicalCuts(float const RapidityCut,
		  int const nPtBins,
		  std::vector<float> const& PtEdge,
		  std::vector<float> const& dcaV0ToPv,
		  std::vector<float> const& decayLength,
		  std::vector<float> const& cosTheta,
		  std::vector<float> const& dcaDaughters,
		  std::vector<float> const& kDca,
		  std::vector<float> const& pDca):
    RapidityCut(RapidityCut),
      nPtBins(nPtBins),
      PtEdge(PtEdge),
      dcaV0ToPv(dcaV0ToPv),
      decayLength(decayLength),
      cosTheta(cosTheta),
      dcaDaughters(dcaDaughters),
      kDca(kDca),
      pDca(pDca)
    {;};
  };
  TopologicalCuts const D0Cuts_50eff(1., //RapidityCut
				     5, //ptBins
				     {0., 1., 2., 3., 5., 15.}, //ptEdge
				     {0.0044, 0.0036, 0.0031, 0.0026, 0.0032}, //dcaV0ToPv
				     {0.0144, 0.0204, 0.0242, 0.0245, 0.0300}, //decayLength
				     {0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, //cosTheta
				     {0.0069, 0.0048, 0.0044, 0.0049, 0.0047}, //dcaDaughters
				     {0.0119, 0.0110, 0.0109, 0.0106, 0.0080}, //kDca
				     {0.0120, 0.0102, 0.0118, 0.0109, 0.0096}); //pDca
  //
  TopologicalCuts const D0Cuts_150eff(1., //RapidityCut
				      5, //ptBins 
				      {0., 1., 2., 3., 5., 15.}, //ptEdge
				      {0.0072, 0.0053, 0.0047, 0.0042, 0.0062}, //dcaV0ToPv
				      {0.0110, 0.0168, 0.0187, 0.0199, 0.0180}, //decayLength
				      {0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, //cosTheta
				      {0.0077, 0.0078, 0.0074, 0.0068, 0.0066}, //dcaDaughters
				      {0.0105, 0.0068, 0.0080, 0.0066, 0.0041}, //kDca
				      {0.0092, 0.0078, 0.0086, 0.0065, 0.0047}); //pDca
  //
  TopologicalCuts const D0Cuts(1., //RapidityCut
			       5, //ptBins 
			       {0., 1., 2., 3., 5., 15.}, //ptEdge
			       {0.0061, 0.0049, 0.0038, 0.0038, 0.0040}, //dcaV0ToPv
			       {0.0145, 0.0181, 0.0212, 0.0247, 0.0259}, //decayLength
			       {0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, //cosTheta
			       {0.0084, 0.0066, 0.0057, 0.0050, 0.0060}, //dcaDaughters
			       {0.0103, 0.0091, 0.0095, 0.0079, 0.0058}, //kDca
			       {0.0110, 0.0111, 0.0086, 0.0081, 0.0062}); //pDca 
}
#endif
