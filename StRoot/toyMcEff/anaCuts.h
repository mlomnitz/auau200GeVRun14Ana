#ifndef anaCuts_H
#define anaCuts_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

namespace anaCuts
{
   int   const nPtBins = 5;
   float const PtEdge[nPtBins+1] = {0., 1., 2., 3., 5., 10.};

   float const pt = 0.6;
   float const eta = 1.0;

   float const rapidity = 1.0;

   float const dcaV0ToPv[nPtBins] = {61, 49, 38, 38, 40};
   float const decayLength[nPtBins] = {145, 181, 212, 247, 259};
   float const cosTheta[nPtBins] = {0.0000, 0.0000, 0.0000, 0.0000, 0.0000};//0.995
   float const dcaDaughters[nPtBins] = {84, 66, 57, 50, 60}; //0.0050;
   float const kDca[nPtBins] = {103, 91, 95, 79, 58};//0.008, // minimum
   float const pDca[nPtBins] = {110, 111, 86, 81, 62};//0.008
}
#endif
