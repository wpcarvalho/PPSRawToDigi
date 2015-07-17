#include "RPRootDumpDigiInfo.h"

RPRootDumpDigiInfo::RPRootDumpDigiInfo() {
	numberOfPlanesOn = uPlanesOn = vPlanesOn = 0;
	numberOfClusters = std::vector<int>();
	numberOfClusters.reserve(10);
	numberOfClusters.assign(10, 0);
	planeId = std::vector<int>();
	clusterSize = std::vector<int>();
	centralStrip = std::vector<int>();
}
