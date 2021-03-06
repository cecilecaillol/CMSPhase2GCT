#ifndef __STITCHTOWER_H__
#define __STITCHTOWER_H__ 

#include <utility>
#include "objects.h"

using namespace std;

static const size_t stitch_in_eta = 0;
static const size_t stitch_in_phi = 1;

void stitchNeighbors(size_t stitchType, Tower inA, Tower inB, Tower &outA, Tower &outB);

#endif /*!__STITCHTOWER_H__ */
