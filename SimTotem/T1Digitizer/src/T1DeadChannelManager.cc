/*
 * T1DeadChannelManager.cc
 *
 *  Created on: Aug 31, 2011
 *      Author: Marcin Borratynski (mborratynski@gmail.com)
 */

#include "SimTotem/T1Digitizer/interface/T1DeadChannelManager.h"
#include "TotemCondFormats/DAQInformation/interface/TotemSymbId.h"
#include "TotemCondFormats/DAQInformation/interface/DAQMapping.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include <map>


T1DeadChannelManager::T1DeadChannelManager() {
	analysisMaskPresent = false;
}

T1DeadChannelManager::T1DeadChannelManager(edm::ESHandle<AnalysisMask> _analysisMask) {
	analysisMask = _analysisMask;
	analysisMaskPresent = true;

	//debug
	//displayMap();
}


bool T1DeadChannelManager::isChannelDead(T1DetId detectorId, unsigned short channelNumber) {

	// TotemSymId.symbolicID in decimal form:  [Arm][Plane][CSC][ChannelType=Layer]
	unsigned int symbolicId =
			detectorId.Layer() 		+
			detectorId.CSC()	* 10 +
			detectorId.Plane()	* 100 +
			detectorId.Arm()	* 1000 ;

	TotemSymbID totemSymbolicId;
	totemSymbolicId.subSystem = TotemSymbID::T1;
	totemSymbolicId.symbolicID = symbolicId;

	if (analysisMaskPresent) {
		std::map<TotemSymbID, VFATAnalysisMask>::const_iterator mask_iterator = analysisMask->analysisMask.find(
		        totemSymbolicId);
		if (mask_iterator != analysisMask->analysisMask.end()) {
			VFATAnalysisMask channelTypeMask = mask_iterator->second;
			//if channel is dead return true
			if (channelTypeMask.fullMask || channelTypeMask.maskedChannels.find(channelNumber)
					!= channelTypeMask.maskedChannels.end()) {

				//debug
				//if(channelTypeMask.fullMask)
				//	std::cout << "---> Masked symbolicID: " << symbolicId << " fullMask" << std::endl ;
				//else
				//	std::cout << "---> Masked symbolicID: " << symbolicId << " channelNumber: " << channelNumber << std::endl;

				return true;
			}
		}
	}
	return false;
}


void T1DeadChannelManager::displayMap() {
	std::cout << "-------------------channel-mask------------------------" << std::endl;
	if (analysisMaskPresent) {
		std::map<TotemSymbID, VFATAnalysisMask>::const_iterator vfatIter;
		for (vfatIter = analysisMask->analysisMask.begin(); vfatIter != analysisMask->analysisMask.end(); vfatIter++) {
			std::cout << vfatIter->first.symbolicID << "\n";
			VFATAnalysisMask am = vfatIter->second;
			if (am.fullMask) {
				std::cout << "   full mask\n";
			} else {
				std::set<unsigned char>::iterator setIterator;
				for (setIterator = am.maskedChannels.begin(); setIterator != am.maskedChannels.end(); setIterator++) {
					std::cout << "   " << (int) (*setIterator) << "\n";
				}
			}

		}
	}
	std::cout << "-----------------------------------------------------" << std::endl;
}
