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

   float const dcaV0ToPv[nPtBins] = {0.0062, 0.0047, 0.0040, 0.0041, 0.0042};
   float const decayLength[nPtBins] = {0.0149, 0.0205, 0.0216, 0.0233, 0.0282};
   float const cosTheta[nPtBins] = {0.0000, 0.0000, 0.0000, 0.0000, 0.0000};//0.995
   float const dcaDaughters[nPtBins] = {0.0082, 0.0070, 0.0056, 0.0065, 0.0065}; //0.0050;
   float const kDca[nPtBins] = {0.0123, 0.0097, 0.0091, 0.0075, 0.0053};//0.008, // minimum
   float const pDca[nPtBins] = {0.0109, 0.0108, 0.0100, 0.0074, 0.0067};//0.008
}
#endif
