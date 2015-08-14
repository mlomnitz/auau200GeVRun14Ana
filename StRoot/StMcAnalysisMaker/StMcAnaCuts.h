#ifndef StMcAnaCuts_H
#define StMcAnaCuts_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */
#include <vector>

#include "Rtypes.h"
#include "StEvent/StEnumerations.h"

namespace McAnaCuts
{
  std::vector<unsigned int> getAllTriggers()
  {
    std::vector<unsigned int> t;
    t.push_back(450050);
    t.push_back(450060); 
    t.push_back(450005); 
    t.push_back(450015); 
    t.push_back(450025);

    return t;
  }

  std::vector<unsigned int> const interesting_triggers = getAllTriggers();

  float const mcTrackStartVtxR = 1.0; // maximum
  int const geantId = 8;

  StDedxMethod dedxMethod = kLikelihoodFitId;

  int const maxNumberOfTriggers = 6;
}
#endif
