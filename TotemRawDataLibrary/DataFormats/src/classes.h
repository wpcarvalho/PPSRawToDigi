/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
****************************************************************************/

#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/OldVFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawVFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/ClusterizationVFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/FramePosition.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrameCollection.h"
#include "TotemRawDataLibrary/DataFormats/interface/SimpleVFATFrameCollection.h"
#include "TotemRawDataLibrary/DataFormats/interface/OptoRxVFATFrameCollection.h"
#include "TotemRawDataLibrary/DataFormats/interface/OptoRxSupplementalData.h"
#include "TotemRawDataLibrary/DataFormats/interface/Raw2DigiStatus.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include <vector>

namespace { namespace {
	edm::Wrapper<Totem::RawEvent> dummy0;

	// the class is abstract
	//Totem::VFATFrame dummy1;
	//std::vector<Totem::VFATFrame> dummy2;

    //todo check if this makes sense (Jan is not certain)
	std::vector<Totem::OldVFATFrame> dummy20;
	edm::Wrapper< Totem::OldVFATFrame > dummy21;

	std::vector<Totem::RawVFATFrame> dummy22;
	edm::Wrapper< Totem::RawVFATFrame > dummy23;

	std::vector<Totem::ClusterizationVFATFrame> dummy24;
	edm::Wrapper< Totem::ClusterizationVFATFrame > dummy25;

	Totem::FramePosition dummy4;
	Totem::OptoRxSupplementalData dummy26;
	edm::Wrapper< Totem::OptoRxSupplementalData > dummy27;

	// the class is abstract
	//Totem::VFATFrameCollection dummy5;
	//edm::Wrapper< Totem::VFATFrameCollection > dummy6;

	Totem::SimpleVFATFrameCollection dummy7;
	edm::Wrapper< Totem::SimpleVFATFrameCollection > dummy8;

	Totem::OptoRxVFATFrameCollection dummy9;
	edm::Wrapper< Totem::OptoRxVFATFrameCollection > dummy10;
	
	Totem::OptoRxMetaData dummy17;
	std::map<unsigned int, Totem::OptoRxMetaData> dummy18;

	Totem::TriggerData dummy19;

	std::map<unsigned int, long> dummy100;
}}
