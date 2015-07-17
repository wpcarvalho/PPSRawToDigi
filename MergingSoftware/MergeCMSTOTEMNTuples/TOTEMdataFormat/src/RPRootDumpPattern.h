#ifndef _RPRootDumpPattern_h_
#define _RPRootDumpPattern_h_

#include <TObject.h>

/**
 * \brief Describes a recognized linear pattern.
 * This class is used only internally inside RPPatternInfo.
 *
 **/
class RPRootDumpPattern : public TObject{
public:
	RPRootDumpPattern(double _a=0., double _b=0., double _w=0.);

	double a; // slope in rad
	double b; // intercept (at the middle of the RP) in mm
	double w; // weight

	ClassDef(RPRootDumpPattern,2);
};

#endif
