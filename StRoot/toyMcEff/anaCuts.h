#ifndef anaCuts_H
#define anaCuts_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

#include "TString.h"

namespace anaCuts
{
   int   const nPtBins = 5;
   float const PtEdge[nPtBins+1] = {0., 1., 2., 3., 5., 10.};

   float const pt = 0.6;
   float const eta = 1.0;

   float const rapidity = 1.0;

   float const massMin = 1.828;
   float const massMax = 1.892;

   float const dcaV0ToPv[nPtBins] = {61, 49, 38, 38, 40};
   float const decayLength[nPtBins] = {145, 181, 212, 247, 259};
   float const cosTheta[nPtBins] = {0.0000, 0.0000, 0.0000, 0.0000, 0.0000};//0.995
   float const dcaDaughters[nPtBins] = {84, 66, 57, 50, 60}; //0.0050;
   float const kDca[nPtBins] = {103, 91, 95, 79, 58};//0.008, // minimum
   float const pDca[nPtBins] = {110, 111, 86, 81, 62};//0.008

   // all variables with a phys prefix are for final physics plots
   int const physNCentralities = 4;
   int const physCentralityEdges[physNCentralities+1] = {0,3,5,6,9}; // 40-80, 20-40, 10-20, 0-10
   TString const physCentralityName[physNCentralities] = {"40-80%","20-40%","10-20%","0-10%"};

   int   const physNPtBins = 10;
   double const physPtEdge[physNPtBins+1] = {0., 0.750, 1.250, 1.750, 2.250, 2.750, 3.750, 4.750 , 5.750, 7., 10.};
}
#endif
