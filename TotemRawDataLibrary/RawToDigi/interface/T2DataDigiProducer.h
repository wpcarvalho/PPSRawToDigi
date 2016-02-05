/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors Mirko Berretti: 
*   This code is a modification of RP code --- Author: Jan Kaspar (jan.kaspar@gmail.com)
*    
* $$RCSfile: T2DataDigiProducer.h,v $: $
* $Revision: 2135 $
* $Date: 2010-02-16 14:26:06 +0100 (Tue, 16 Feb 2010) $
*
****************************************************************************/

#ifndef T2DataDigiProducer_h
#define T2DataDigiProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "TotemRawData/RawToDigi/interface/T2Mapping.h"
#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "DataFormats/T2DigiVfat/interface/T2DigiVfat.h"
#include "DataFormats/T2DigiVfat/interface/T2DigiVfatCollection.h"


#include <map>
//#include <boost>
#include <memory>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
namespace edm {
	class ParameterSet;
	class EventSetup;
	class Event;
}

/**
 * \brief Convert raw event data to digi information
**/
class T2DataDigiProducer : public edm::EDProducer
{
	public:
	
		T2DataDigiProducer(const edm::ParameterSet& conf); 
		virtual	~T2DataDigiProducer();
	
		virtual void beginJob(const edm::EventSetup&);	
		virtual void endJob();	
		virtual void produce(edm::Event&, const edm::EventSetup&);

		/*
		 std::auto_ptr<TH1F> Padconvertedindex;
		 std::auto_ptr<TH1F> PadCOLindex;
		 std::auto_ptr<TH1F> PadROWindex;
		 std::auto_ptr<TH1F> HistActiveChannel;
		*/
	private:
		/// flag whether the VFAT frame position index information was given (and shall be used)
		bool positions;	

		/// vector of SLINK-positions to be masked
		std::vector<unsigned int> posMask;

		/// Vector to store the rawids of the vfat
		std::vector<unsigned int> rawIds;

		// map SLINK-position index --> The cms Id of the detector where the VFAT is placed
		std::map<unsigned int, uint32_t> cmsIds; 

		//map SLINK-position index --> raw ID (ChipID)
		std::map<unsigned int, unsigned int> posToRawID;
	       		
		// map raw ID (ChipID) --> symbolic (CMSSW) ID
                std::map<unsigned int, uint32_t> rawVFATtoCMSplane;

		/// map raw ID (ChipID) --> Giuseppe Convenction ID
		std::map<unsigned int, unsigned int> rawToSymID;

		T2Mapping VFATconvertDigichannell;
		unsigned char verbosity;
};
  
#endif
