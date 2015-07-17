#ifndef _NtupleInfo_h_
#define _NtupleInfo_h_

#include <TObject.h>
#include <vector>
#include <map>

class NtupleInfo: public TObject {
public:
	NtupleInfo();
	int minOrbit;
	int maxOrbit;
	std::map<unsigned int, std::vector<int> > orbits;

	ClassDef(NtupleInfo, 1);
};

#endif
