#include "SimTotem/RPDigiProducer/interface/DeadChannelsManager.h"
#include "TotemCondFormats/DAQInformation/interface/TotemSymbId.h"
#include "TotemCondFormats/DAQInformation/interface/DAQMapping.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include <map>

DeadChannelsManager::DeadChannelsManager() {
	analysisMaskPresent = false;
}

DeadChannelsManager::DeadChannelsManager(edm::ESHandle<AnalysisMask> _analysisMask) {
	analysisMask = _analysisMask;
	analysisMaskPresent = true;
}

bool DeadChannelsManager::isChannelDead(RPDetId detectorId, unsigned short stripNumber) {
	unsigned int symbolicId = TotRPDetId::RawToDecId(detectorId) * 10; //convert to symbolic ID
	unsigned int vfat = stripNumber / 128;
	symbolicId += vfat; //add vfatID to symbolic ID
	stripNumber = stripNumber - vfat * 128; //convert strip number to a number from range <0; 127>
	TotemSymbID totemSymbolicId;
	totemSymbolicId.subSystem = TotemSymbID::RP;
	totemSymbolicId.symbolicID = symbolicId;
	if (analysisMaskPresent) {
		std::map<TotemSymbID, VFATAnalysisMask>::const_iterator vfatIter = analysisMask->analysisMask.find(
		        totemSymbolicId);
		if (vfatIter != analysisMask->analysisMask.end()) {
			VFATAnalysisMask vfatMask = vfatIter->second;
			//if channel is dead return true
			if (vfatMask.fullMask || vfatMask.maskedChannels.find(stripNumber)
			        != vfatMask.maskedChannels.end()) {
				return true;
			}
		}
	}
	return false;
}

void DeadChannelsManager::displayMap() {
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
}

