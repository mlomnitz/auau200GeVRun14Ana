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

namespace McAnaCuts
{
  std::vector<unsigned int> interesting_triggers;
  float const mcTrackStartVtxR = 1.0; // maximum
  int const geantId = 8;
}
#endif
