/**
 * This file stores some common variables and definitions for merge12 and combine programs.
 */
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <map>
#include "TriggerData.h"
#include "EventMetaData.h"
#include "T1Event.h"
#include "T2Event.h"
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpDigiInfo.h"
#include "RPRootDumpPatternInfo.h"
#include "RPRootDumpTrackInfo.h"
#include "RPRootDumpReconstructedProtonPair.h"
#include <string>
#include <TTree.h>
#include <iostream>
#include <cstdlib>

using namespace std;

static const int dTOTEMoRxBX = 263; // Valentina
static const unsigned int maxEvtPerOrbit = 10;
const bool cleanBX = false; //clean with a given bcn pattern

const bool includeDigi = false; //indicates whether digi/patterns information should be stored in output ntuple
const bool includePatterns = false;
const bool includeMultiTracks = false;

//some RP variables (important only internally, not for user)
RPRootDumpReconstructedProton *recProtLeft = 0, *recProtRight = 0;
RPRootDumpReconstructedProtonPair *recProtonPair = 0;
map<unsigned int, RPRootDumpDigiInfo *> digi_info_;
map<unsigned int, RPRootDumpPatternInfo *> par_patterns_info_,
		nonpar_patterns_info_;
map<unsigned int, RPRootDumpTrackInfo *> track_info_;
map<unsigned int, std::vector<RPRootDumpTrackInfo> *> multi_track_info_;
// some T1/T2/trigger variables (important only internally, not for user)
T2Event *t2Event = 0;
T1Event *t1Event = 0;
TriggerData *trigData = 0;
EventMetaData *totemEventData = 0;

const unsigned int NofCollidingBX = 3;
const unsigned int collidingBX[NofCollidingBX] = { // cms notation - starting with bcn == 1
		1, 101, 1886 }; // 994 - non-colliding

/**
 * This structure contains basic totem information (trigger + orbit, bunch etc). It is used
 * during searching for offset and merging to speed up the procedure.
 */
struct info {
	info() :
			i(0), evId(0), orbit(0), bcn(0) {
	}
	info(const int theI, const int theEvId,
			const unsigned int theOrbit, const int theBcn) :
			i(theI), evId(theEvId), orbit(theOrbit), bcn(
					theBcn) {
	}
	void dump() const {
		std::cout << "orbit=" << orbit << " bunch=" << bcn << std::endl;
	}
	int i;
	int evId;
	unsigned int orbit;
	int bcn;
};

/**
 * This function writes error message to the standard output and then terminates the program.
 */
void error(string errorMessage) {
	cout << errorMessage << "\n";
	exit(-1);
}

/**
 * Checks if the bunch CMSbcn is colliding bunch.
 */
bool collidingBunch(unsigned int CMSbcn) {
	bool found = false;
	unsigned int bunchPosition = 0;
	while (!found && bunchPosition < NofCollidingBX) {
		if (collidingBX[bunchPosition] == CMSbcn)
			found = true;
		bunchPosition++;
	}
	return found;
}

#endif /* UTILITIES_H_ */
