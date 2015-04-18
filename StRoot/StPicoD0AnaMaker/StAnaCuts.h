#ifndef StAnaCuts_H
#define StAnaCuts_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

#include "Rtypes.h"
#include <string>

namespace anaCuts
{
   extern std::string const prescalesFilesDirectoryName;
   //event
   extern UShort_t const triggerWord;

   //tracking
   extern int const nHitsFit;

   //pions
   extern float const nSigmaPion;

   //kaons
   extern float const nSigmaKaon;
   extern float const kTofBetaDiff;

   extern float const cosTheta;
   extern float const dcaDaughters;
   extern float const kDca;
   extern float const pDca;
}
#endif
