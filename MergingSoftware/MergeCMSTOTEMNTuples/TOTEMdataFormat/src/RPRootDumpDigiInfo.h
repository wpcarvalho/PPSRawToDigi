#ifndef _RPRootDumpDigiInfo_h_
#define _RPRootDumpDigiInfo_h_

#include <vector>
#include <TObject.h>

/**
 * RPRootDumpDigiInfo is stored in TotemNtuple separately for every pot.
 * The name convention is: digi_rp_$number.
 * Where number (pot number) is: 0,1,2,3,4,5 20,21,22,23,24,25 100,101,102,103,104,105 120,121,122,123,124,125
 */
class RPRootDumpDigiInfo: public TObject {
public:
	RPRootDumpDigiInfo();

	std::vector<int> numberOfClusters; // number of clusters in a given plane (indexed from 0 to 9)
	unsigned int numberOfPlanesOn; // number of planes with at least one cluster
	unsigned int uPlanesOn; // number of U planes with at least one cluster
	unsigned int vPlanesOn; // number of V planes with at least one cluster
	std::vector<int> planeId; // plane ID for a given cluster (array index)
	std::vector<int> clusterSize; // cluster size of a given cluster
	std::vector<int> centralStrip; // central strip of a given cluster

	ClassDef(RPRootDumpDigiInfo,2);
};

#endif
