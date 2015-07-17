/*
 * T1DeadChannelManager.h
 *
 *  Created on: Aug 31, 2011
 *      Author: Marcin Borratynski (mborratynski@gmail.com)
 *
 * Purpose of this class is to answer the question whether a channel (given by detectorId
 * and stripNumber) is dead or not. This class uses analysisMask which is provided
 * by DAQMappingSourceXML.
 */

#ifndef T1DEADCHANNELMANAGER_H_
#define T1DEADCHANNELMANAGER_H_

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "TotemCondFormats/DAQInformation/interface/AnalysisMask.h"

class T1DeadChannelManager {

private:
	edm::ESHandle<AnalysisMask> analysisMask;
	bool analysisMaskPresent; //this variable indicates whether analysisMask is present or not

public:
	T1DeadChannelManager();
	T1DeadChannelManager(edm::ESHandle<AnalysisMask> analysisMask);

	/**
	 * This function answers the question whether given channel is dead or not.
	 * T1DetId - detector ID given in raw form, this function has to convert raw ID to symbolic
	 * stripNumber - ID of the strip, it is a number from range <0; 225>
	 *
	 * It is assumed that in T1DetId:
	 *  Arm		is a number form range  <0; 1>
	 *  Plane	is a number form range  <0; 4>
	 *  CSC		is a number form range  <0; 5>
	 *  Layer	is a number form range  <1; 3>
	 *
	 *  TotemSymId.symbolicID in decimal form:  [Arm][Plane][CSC][ChannelType=Layer]
	 *
	 */
	bool isChannelDead(T1DetId detectorId, unsigned short stripNumber);
	void displayMap();


};

#endif /* T1DEADCHANNELMANAGER_H_ */
