#ifndef CUTS_H
#define CUTS_H

#include "StCuts.h"

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */


namespace cuts
{
   // path to lists of triggers prescales
   // lists are obtained from http://www.star.bnl.gov/protected/common/common2014/trigger2014/plots_au200gev/
   std::string const prescalesFilesDirectoryName = "./run14AuAu200GeVPrescales";

   //event
   UShort_t const triggerWord = 0x1F; //first five bits see http://rnc.lbl.gov/~xdong/SoftHadron/picoDst.html
   float const vz = 6.0;// cm.
   float const vzVpdVz = 3.0; // 3 cm.

   // vertex refit track quality
   float const vtxDca = 3.0;
   size_t const vtxNumberOfFitPoints = 20;

   int const nPxlPhiEdges = 11;
   float const PxlPhiEdges[nPxlPhiEdges + 1] = {-3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159};
   int const PxlSectorNumber[nPxlPhiEdges + 1] = {3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3};

   int getPxlSectorNumber(float const phi)
   {
     for (int i = 0; i < nPxlPhiEdges; ++i)
     {
       if (phi > PxlPhiEdges[i] && phi <= PxlPhiEdges[i+1])
       {
         return PxlSectorNumber[i];
       }
     }

     return 0;
   }
}
#endif
