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

   float const dcaV0ToPv[nPtBins] = {62, 47, 40, 41, 42};
   float const decayLength[nPtBins] = {149, 205, 216, 233, 282};
   float const cosTheta[nPtBins] = {0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
   float const dcaDaughters[nPtBins] = {82, 70, 56, 65, 65};
   float const kDca[nPtBins] = {123, 97, 91, 75, 53};
   float const pDca[nPtBins] = {109,108,100, 74, 67};
}
#endif
