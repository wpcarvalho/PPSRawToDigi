#ifndef _RPRootDumpPatternInfo_h_
#define _RPRootDumpPatternInfo_h_

#include "RPRootDumpPattern.h"
#include <TObject.h>

/**
 * \brief Describes a pattern-recognition result.
 * This class stores pattern-recognition data. In TotemNtuple it is possible to have objects called:
 * par_patterns_rp_$number. or nonpar_patterns_rp_$number. (where $number indicates pot).
 * $number can be: 0,1,2,3,4,5 20,21,22,23,24,25 100,101,102,103,104,105 120,121,122,123,124,125
 **/
struct RPRootDumpPatternInfo: public TObject {
	RPRootDumpPatternInfo(bool _f = false);

	std::vector<RPRootDumpPattern> u, v; // arrays of recognized patterns in u and v projections
	bool fittable; // whether there is one (and only one) combined u-v pattern worth fitting

	ClassDef(RPRootDumpPatternInfo,2);
};

#endif
