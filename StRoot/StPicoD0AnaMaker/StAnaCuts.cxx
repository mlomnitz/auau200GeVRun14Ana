/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

#include "StAnaCuts.h"

namespace anaCuts
{
   // path to lists of triggers prescales
   // lists are obtained from http://www.star.bnl.gov/protected/common/common2014/trigger2014/plots_au200gev/
   std::string const prescalesFilesDirectoryName = "./run14AuAu200GeVPrescales";
   //event
   UShort_t const triggerWord = 0x1F; //first five bits see http://rnc.lbl.gov/~xdong/SoftHadron/picoDst.html

   //tracking
   float const minPt = 1.2;
   int const nHitsFit = 20;

   //pions
   float const nSigmaPion = 3.0;

   //kaons
   float const nSigmaKaon = 2.0;
   float const kTofBetaDiff = 0.04;

   float const cosTheta = 0.995;
   float const dcaDaughters = 0.0050;
   float const kDca = 0.008; // minimum
   float const pDca = 0.008;
}
