/*
 * RPCoincidenceAnalyzer.cc
 *
 *  Created on: Dec 9, 2008
 *      Author: Leszek Grzanka
 */

#include <stdio.h>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "L1TriggerTotem/CoincidenceChip/interface/RPCoincidenceAnalyzer.h"

#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/TotemRPDataTypes/interface/RPRecoHit.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidate.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"


//----------------------------------------------------------------------------------------------------

RPCoincidenceAnalyzer::RPCoincidenceAnalyzer(const edm::ParameterSet& conf)
{
	verbose_ = conf.getParameter<unsigned int> ("verbose");

	// output file, where histograms will be saved
	outputFile = conf.getParameter<std::string>("outputFile");

	modulLabelRaw    = conf.getParameter<std::string>("modulLabelRaw");
	productLabelRaw  = conf.getParameter<std::string>("productLabelRaw");
	modulLabelSimu   = conf.getParameter<std::string>("modulLabelSimu");
	productLabelSimu = conf.getParameter<std::string>("productLabelSimu");
	trackCandidateCollectionLabel = conf.getParameter<edm::InputTag>("TrackCandidateCollectionLabel");
	detTriggerLabel =  conf.getParameter<edm::InputTag>("DetTriggerLabel");
}

//----------------------------------------------------------------------------------------------------

RPCoincidenceAnalyzer::~RPCoincidenceAnalyzer()
{		
}

//----------------------------------------------------------------------------------------------------

void RPCoincidenceAnalyzer::beginJob()
{
    if( verbose_ ) edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer: beginJob";
    	
	trigSectorsRawEven=std::auto_ptr<TH1D>(new TH1D("triggeredSectorsRawVdir","Triggered sectors in real chip - even planes (0,2,..8);Sector number;",16,-0.5,15.5));
  	trigSectorsRawEven->SetDirectory(0);
	trigSectorsRawOdd=std::auto_ptr<TH1D>(new TH1D("triggeredSectorsRawUdir","Triggered sectors in real chip - odd planes (1,3,..9);Sector number;",16,-0.5,15.5));
  	trigSectorsRawOdd->SetDirectory(0);
	trigSectorsSimuEven=std::auto_ptr<TH1D>(new TH1D("triggeredSectorsSimVdir","Triggered sectors in simulated chip - even planes (0,2,..8);Sector number;",16,-0.5,15.5));
  	trigSectorsSimuEven->SetDirectory(0);
	trigSectorsSimuOdd=std::auto_ptr<TH1D>(new TH1D("triggeredSectorsSimUdir","Triggered sectors in simulated chip - odd planes (1,3,..9);Sector number;",16,-0.5,15.5));
  	trigSectorsSimuOdd->SetDirectory(0);
 
 	differenceSimuVsRawEven=std::auto_ptr<TH1I>(new TH1I("differenceSimuVsRawVdir","Difference Simulated vs Real - even planes (0,2,..8);Number of different bits;",17,-0.5,16.5));
  	differenceSimuVsRawEven->SetDirectory(0);
	differenceSimuVsRawOdd=std::auto_ptr<TH1I>(new TH1I("differenceSimuVsRawUdir","Difference Simulated vs Real - odd planes (1,3,..9);Number of different bits;",17,-0.5,16.5));
  	differenceSimuVsRawOdd->SetDirectory(0);
	differenceSimuVsRawTotal=std::auto_ptr<TH1I>(new TH1I("differenceSimuVsRawTotal","Difference Simulated vs Real - all planes;Number of different bits;",33,-0.5,32.5));
  	differenceSimuVsRawTotal->SetDirectory(0);

	differenceCategories=std::auto_ptr<TH1I>(new TH1I("differenceCategories","Difference Categories",16,-0.5,15.5));
  	differenceCategories->SetDirectory(0);

	tracksSelectionCategories = std::auto_ptr<TH1I>(new TH1I("tracksSelectionCategories","Track selection categories",4,-0.5,3.5));
  	tracksSelectionCategories->SetDirectory(0);

	goodTracksCategories=std::auto_ptr<TH1I>(new TH1I("goodTracksCategories","Good tracks categories",8,-0.5,7.5));
  	goodTracksCategories->SetDirectory(0);

	goodEventsCategories=std::auto_ptr<TH1I>(new TH1I("goodEventsCategories","Good events categories",4,-0.5,3.5));
	goodEventsCategories->SetDirectory(0);

	activeChipSim=std::auto_ptr<TH1I>(new TH1I("activeChipSim","Active Chip Simulation",4,-0.5,3.5));
  	activeChipSim->SetDirectory(0);

	activeChipRaw=std::auto_ptr<TH1I>(new TH1I("activeChipRaw","Active Chip Real",4,-0.5,3.5));
  	activeChipRaw->SetDirectory(0);

	totalEfficiencyRaw=std::auto_ptr<TH1I>(new TH1I("totalEfficiencyRaw","Coincidence efficiency in real chip",5,-0.5,4.5));
  	totalEfficiencyRaw->SetDirectory(0);
	totalEfficiencySimu=std::auto_ptr<TH1I>(new TH1I("totalEfficiencySimu","Coincidence efficiency in simulated chip",5,-0.5,4.5));
  	totalEfficiencySimu->SetDirectory(0);
}

//----------------------------------------------------------------------------------------------------


/*
 * One event contains data from: 10 detectors and 2 coincidence chips
 */
void RPCoincidenceAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  using namespace edm;
  using namespace std;

  // Read real strips : output of real CC and raw data
  Handle<vector<RPCCBits> > inputRaw;
  event.getByLabel(modulLabelRaw, productLabelRaw, inputRaw ); 
  
  // We will also need trigger bits for some printouts
  Handle<DetSetVector<RPDetTrigger> > inputTrigger;
  event.getByLabel(detTriggerLabel, inputTrigger );
  
  // Read output of CC simulation module
  Handle<vector<RPCCBits> > inputSimu;
  event.getByLabel(modulLabelSimu, productLabelSimu, inputSimu ); 

  // Read good tracks
  Handle<RPTrackCandidateCollection> selHits;
  event.getByLabel(trackCandidateCollectionLabel, selHits);

  /*
   * We will need evenRawCCBits, oddRawCCBits, evenSimuCCBits, oddSimuCCBits
   * objects later, for comparison of raw CC and simu CC
   */
  RPCCBits evenRawCCBits;
  RPCCBits oddRawCCBits;
  RPCCBits evenSimuCCBits;
  RPCCBits oddSimuCCBits;

  /*
   * Let us initialize "simu" objects with zeros
   * it might happen that inputSimu collection is empty
   * and we will not step into second for loop
   */
  std::bitset<16> nullBitset;
  nullBitset.reset();
  evenSimuCCBits.setBS(nullBitset);
  evenSimuCCBits.setId(TotRPDetId(1,2,0,0));
  oddSimuCCBits.setBS(nullBitset);
  oddSimuCCBits.setId(TotRPDetId(1,2,0,1));

  evenRawCCBits.setBS(nullBitset);
  evenRawCCBits.setId(TotRPDetId(1,2,1,0));
  oddRawCCBits.setBS(nullBitset);
  oddRawCCBits.setId(TotRPDetId(1,2,1,1));

  if( verbose_ > 2) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << ": INPUT - RPCCRawDigi objects : " << inputRaw->size();
  if( verbose_ > 2) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << ": INPUT - RPCCSimuBits objects : " << inputSimu->size();


  vector<RPCCBits>::const_iterator inputRawIterator = inputRaw->begin();
  for (; inputRawIterator != inputRaw->end(); inputRawIterator++) {
  	
  	TotRPDetId detectorId(inputRawIterator->getId());

        if( verbose_ ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer detId " << detectorId;
    /*
     * Select RPCCBits object which belongs to VFAT 1200 having simulated data from even CC (V direction)
     * There is only one such VFAT in inputSimu collection
     */
	if( (detectorId.Arm() == 1) && (detectorId.Station() == 2) && (detectorId.RomanPot() == 0) && (detectorId.Detector() == 0) ){

           if( verbose_ ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " analyzing data from real CC - even (0,2,..8)";
           if( verbose_ ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " : CC output bits (even): " << inputRawIterator->getBS();

		/*
		 * Fill histogram with triggered sectors from even simulated CC
		 */ 
		for( unsigned int i = 0 ; i < (inputRawIterator->getBS()).size() ; i++ ){
			if( (inputRawIterator->getBS()).test(i) ){
				trigSectorsRawEven->Fill(i);
			}
		}

		/*
     	 * Store simulated-CC output for later comparison with raw-CC output in evenSimuCCBits object  
		 */
		evenRawCCBits = *inputRawIterator;
		
	}

    /*
     * Select RPCCBits object which belongs to VFAT 1201 having simulated data from odd CC (U direction)
     * There is only one such VFAT in inputSimu collection
     */
	if( (detectorId.Arm() == 1) && (detectorId.Station() == 2) && (detectorId.RomanPot() == 0) && (detectorId.Detector() == 1) ){

           if( verbose_ ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " analyzing data from real CC - odd (1,3,..9)";
    	   if( verbose_ ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " : CC output bits (odd): " << inputRawIterator->getBS();

		/*
		 * Fill histogram with triggered sectors from odd simulated CC
		 */ 
		for( unsigned int i = 0 ; i < (inputRawIterator->getBS()).size() ; i++ ){
			if( (inputRawIterator->getBS()).test(i) ){
				trigSectorsRawOdd->Fill(i);
			}
		}

		/*
       	 * Store simulated-CC output for later comparison with raw-CC output in oddSimuCCBits object  
		 */
		oddRawCCBits = *inputRawIterator;
	}
	
   }



   /*                           
    * ***************************************
    * ***** READ DATA FROM SIMULATED CC *****
    * ***************************************
    * 
    * Data from real detectors (10 VFATs) were taken, RPTriggerBits were calculated and served as input to CC simulator
    * 10 Real VFATs  =>  RPTriggerBits => 2 Simulated CC 
    * 
    * This loop is also dummy, we will just select two items which corresponds to VFAT 1200 and 1201
    */
  vector<RPCCBits>::const_iterator inputSimuIterator = inputSimu->begin();
  for (; inputSimuIterator != inputSimu->end(); inputSimuIterator++) {
  	
  	TotRPDetId detectorId(inputSimuIterator->getId());

    /*
     * Select RPCCBits object which belongs to VFAT 1200 having simulated data from even CC (V direction)
     * There is only one such VFAT in inputSimu collection
     */
	if( (detectorId.Arm() == 1) && (detectorId.Station() == 2) && (detectorId.RomanPot() == 0) && (detectorId.Detector() == 0) ){

           if( verbose_ ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " analyzing data from simulated CC - even (0,2,..8)";
           if( verbose_ ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " : CC output bits (even): " << inputSimuIterator->getBS();

		/*
		 * Fill histogram with triggered sectors from even simulated CC
		 */ 
		for( unsigned int i = 0 ; i < (inputSimuIterator->getBS()).size() ; i++ ){
			if( (inputSimuIterator->getBS()).test(i) ){
				trigSectorsSimuEven->Fill(i);
			}
		}

		/*
     	 * Store simulated-CC output for later comparison with raw-CC output in evenSimuCCBits object  
		 */
		evenSimuCCBits = *inputSimuIterator;
		
	}

    /*
     * Select RPCCBits object which belongs to VFAT 1201 having simulated data from odd CC (U direction)
     * There is only one such VFAT in inputSimu collection
     */
	if( (detectorId.Arm() == 1) && (detectorId.Station() == 2) && (detectorId.RomanPot() == 0) && (detectorId.Detector() == 1) ){

           if( verbose_ ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " analyzing data from simulated CC - odd (1,3,..9)";
    	   if( verbose_ ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " : CC output bits (odd): " << inputSimuIterator->getBS();

		/*
		 * Fill histogram with triggered sectors from odd simulated CC
		 */ 
		for( unsigned int i = 0 ; i < (inputSimuIterator->getBS()).size() ; i++ ){
			if( (inputSimuIterator->getBS()).test(i) ){
				trigSectorsSimuOdd->Fill(i);
			}
		}

		/*
     	 * Store simulated-CC output for later comparison with raw-CC output in oddSimuCCBits object  
		 */
		oddSimuCCBits = *inputSimuIterator;
	}
	
   }

   /* 
    * **************************************************
    * ***  CHECK IF THERE IS GOOD TRACK IN THE EVENT ***
    * **************************************************
    **/ 
  int fittableTrack = -1;
  if( selHits->size() > 0 ){
	 LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << ":  found " << (selHits->size()) << " good track";

        // process selected hits collection
        //map< unsigned int, pair< vector<const RPRecoHit *>, vector<const RPRecoHit *> > > selHitsMap;
        for (RPTrackCandidateCollection::const_iterator dit = selHits->begin(); dit != selHits->end(); ++dit) {
                if (!dit->second.Fittable()){
	 		LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << ":  track nonfittable";
			fittableTrack = 0;
			tracksSelectionCategories->Fill(1);
		} else {
	 		LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << ":  track fittable";
			fittableTrack = 1;
			tracksSelectionCategories->Fill(0);
		}

                //unsigned int RPId = dit->first;
                //const vector<RPRecoHit> &rhs = dit->second.TrackRecoHits();
                //for (vector<RPRecoHit>::const_iterator hit = rhs.begin(); hit != rhs.end(); ++hit) {
                //        unsigned int detId = TotRPDetId::RawToDecId(hit->DetId());
                //        bool uDir = TotRPDetId::IsStripsCoordinateUDirection(detId);

                  //      if (uDir) selHitsMap[RPId].first.push_back(& (*hit));
                //        else selHitsMap[RPId].second.push_back(& (*hit));
               // }
        }

  } 

   /* 
    * *************************************
    * ***  COMPARE REAL vs SIMULATED CC ***
    * *************************************
    * 
    * Compare output bits of Raw and Simulated Coincidence Chips for U and V direct.
    * If two output patterns are the same, we say they differs on 0 bits.
    * If they are completely different, we say they differs on 16 bits.
    * In following section we will create histogram of number of bits on which raw and simu patterns differs.
    * 
    */

	// calculate XOR of raw and simu bits patterns and count how many one's we have
	int numberOfDifferentBitsOdd = (oddSimuCCBits.getBS()^oddRawCCBits.getBS()).count();
	differenceSimuVsRawOdd->Fill(numberOfDifferentBitsOdd);
    
	int numberOfDifferentBitsEven = (evenSimuCCBits.getBS()^evenRawCCBits.getBS()).count();
	differenceSimuVsRawEven->Fill(numberOfDifferentBitsEven);

	// here we calculate total number of bits on which pattern differs
	//differenceSimuVsRawTotal->Fill(numberOfDifferentBitsEven+numberOfDifferentBitsOdd);
        //if ( (fittableTrack >= 0) && ( (oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) > 0 ) && (numberOfDifferentBitsEven + numberOfDifferentBitsOdd > 0) ){
        if ( (fittableTrack >= 0) && (numberOfDifferentBitsEven + numberOfDifferentBitsOdd > 0) ){
		differenceSimuVsRawTotal->Fill(numberOfDifferentBitsEven+numberOfDifferentBitsOdd);
	}

        if( verbose_ ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " comparison real vs simulated";


	// if one of the pattern differs, we print some info
	//if( (numberOfDifferentBitsOdd + numberOfDifferentBitsEven > 0) && (verbose_) ){
	if( (verbose_) ){
	    LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " : different in output between real and simulated CC : ";		
	    LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " : odd (1,3,..9)  simulation = " << oddSimuCCBits.getBS();
	    LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " : odd (1,3,..9)    hardware = " << oddRawCCBits.getBS();
	    LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " : even (0,2,..8) simulation = " << evenSimuCCBits.getBS();
	    LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " : even (0,2,..8)   hardware = " << evenRawCCBits.getBS();

	    // print Trigger bits
	    DetSetVector<RPDetTrigger>::const_iterator inputTriggerIterator = inputTrigger->begin();
	    std::bitset<16> triggerBits;
  	    for (; inputTriggerIterator != inputTrigger->end(); inputTriggerIterator++) {
    		  TotRPDetId detectorId(inputTriggerIterator->id);
    		  triggerBits.reset();
    		  unsigned int triggeredSectorNo = 0;
    		  for (unsigned int i = 0; i < (inputTriggerIterator->data).size(); ++i) {
      			triggeredSectorNo = (unsigned int) ((inputTriggerIterator->data)[i].GetSector());
      			triggerBits.set(triggeredSectorNo);
    		   }
    		   if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " : detector " << detectorId.Arm() << detectorId.Station() << detectorId.RomanPot() << detectorId.Detector() << " , trigger bits : " << triggerBits;
	    }
	}
	
   /* 
    * *************************************
    * ***  CALCULATE EFFICIENCY  ***
    * *************************************
    */ 

  // Check if event is empty (input to CC consists of zeros)
  DetSetVector<RPDetTrigger>::const_iterator inputTriggerIterator = inputTrigger->begin();
  int trigBitsNonZeroCounter = 0;
  for (; inputTriggerIterator != inputTrigger->end(); inputTriggerIterator++) {
        if( (inputTriggerIterator->data).size() > 0 ){
		trigBitsNonZeroCounter++;
	}
  }

  if( (fittableTrack == -1) && (trigBitsNonZeroCounter == 0) ){
	tracksSelectionCategories->Fill(3);
  }
  if( (fittableTrack == -1) && (trigBitsNonZeroCounter > 0) ){
	tracksSelectionCategories->Fill(2);
  }

  if( trigBitsNonZeroCounter == 0 ){
        LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " :  EventEmpty ";
	totalEfficiencySimu->Fill(4);
	totalEfficiencyRaw->Fill(4);
  }

   if(( oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) > 0 ){
	totalEfficiencyRaw->Fill(0);
   } else if ( oddRawCCBits.getBS().count() > 0 )  {
	totalEfficiencyRaw->Fill(1);
   } else if ( evenRawCCBits.getBS().count() > 0 )  {
	totalEfficiencyRaw->Fill(2);
   } else if ( trigBitsNonZeroCounter > 0 ) {
	totalEfficiencyRaw->Fill(3);
   }

	// calculate XOR of raw and simu bits patterns and count how many one's we have
   numberOfDifferentBitsOdd = (oddSimuCCBits.getBS()^oddRawCCBits.getBS()).count();
   numberOfDifferentBitsEven = (evenSimuCCBits.getBS()^evenRawCCBits.getBS()).count();


   // both = 0
   if( (evenSimuCCBits.getBS().count() == 0 ) && (evenRawCCBits.getBS().count() == 0 ) && (oddSimuCCBits.getBS().count() == 0 ) && (oddRawCCBits.getBS().count() == 0 ) ){
	//differenceCategories->Fill(0); // 0000
        //if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 0 , 0000";
   } else if( (evenSimuCCBits.getBS().count() == 0 ) && (evenRawCCBits.getBS().count() == 0 ) && (oddSimuCCBits.getBS().count() == 0 ) && (oddRawCCBits.getBS().count() > 0 ) ){
	differenceCategories->Fill(1); // 0001
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 1 , 0001";
   } else if( (evenSimuCCBits.getBS().count() == 0 ) && (evenRawCCBits.getBS().count() == 0 ) && (oddSimuCCBits.getBS().count() > 0 ) && (oddRawCCBits.getBS().count() == 0 ) ){
	differenceCategories->Fill(2); // 0010
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 2 , 0010";
   } else if( (evenSimuCCBits.getBS().count() == 0 ) && (evenRawCCBits.getBS().count() == 0 ) && (oddSimuCCBits.getBS().count() > 0 ) && (oddRawCCBits.getBS().count() > 0 ) && (numberOfDifferentBitsOdd > 0) ){
	differenceCategories->Fill(3); // 0011 !
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 3 , 0011";
   } else if( (evenSimuCCBits.getBS().count() == 0 ) && (evenRawCCBits.getBS().count() > 0 ) && (oddSimuCCBits.getBS().count() == 0 ) && (oddRawCCBits.getBS().count() == 0 ) ){
	differenceCategories->Fill(4); // 0100
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 4 , 0100";
   } else if( (evenSimuCCBits.getBS().count() == 0 ) && (evenRawCCBits.getBS().count() > 0 ) && (oddSimuCCBits.getBS().count() == 0 ) && (oddRawCCBits.getBS().count() > 0 ) ){
	differenceCategories->Fill(5); // 0101
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 5 , 0101";
   } else if( (evenSimuCCBits.getBS().count() == 0 ) && (evenRawCCBits.getBS().count() > 0 ) && (oddSimuCCBits.getBS().count() > 0 ) && (oddRawCCBits.getBS().count() == 0 ) ){
	differenceCategories->Fill(6); // 0110
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 6 , 0110";
   } else if( (evenSimuCCBits.getBS().count() == 0 ) && (evenRawCCBits.getBS().count() > 0 ) && (oddSimuCCBits.getBS().count() > 0 ) && (oddRawCCBits.getBS().count() > 0 )){
	differenceCategories->Fill(7); // 0111 !
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 7 , 0111";
   } else if( (evenSimuCCBits.getBS().count() > 0 ) && (evenRawCCBits.getBS().count() == 0 ) && (oddSimuCCBits.getBS().count() == 0 ) && (oddRawCCBits.getBS().count() == 0 ) ){
	differenceCategories->Fill(8); // 1000
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 8 , 1000";
   } else if( (evenSimuCCBits.getBS().count() > 0 ) && (evenRawCCBits.getBS().count() == 0 ) && (oddSimuCCBits.getBS().count() == 0 ) && (oddRawCCBits.getBS().count() > 0 ) ){
	differenceCategories->Fill(9); // 1001
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 9 , 1001";
   } else if( (evenSimuCCBits.getBS().count() > 0 ) && (evenRawCCBits.getBS().count() == 0 ) && (oddSimuCCBits.getBS().count() > 0 ) && (oddRawCCBits.getBS().count() == 0 ) ){
	differenceCategories->Fill(10); // 1010
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 10 , 1010";
   } else if( (evenSimuCCBits.getBS().count() > 0 ) && (evenRawCCBits.getBS().count() == 0 ) && (oddSimuCCBits.getBS().count() > 0 ) && (oddRawCCBits.getBS().count() > 0 )){
	differenceCategories->Fill(11); // 1011 !
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 11 , 1011";
   } else if( (evenSimuCCBits.getBS().count() > 0 ) && (evenRawCCBits.getBS().count() > 0 ) && (oddSimuCCBits.getBS().count() == 0 ) && (oddRawCCBits.getBS().count() == 0 ) && (numberOfDifferentBitsEven > 0) ){
	differenceCategories->Fill(12); // 1100 !!
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 12 , 1100";
   } else if( (evenSimuCCBits.getBS().count() > 0 ) && (evenRawCCBits.getBS().count() > 0 ) && (oddSimuCCBits.getBS().count() == 0 ) && (oddRawCCBits.getBS().count() > 0 )){
	differenceCategories->Fill(13); // 1101 !
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 13 , 1101";
   } else if( (evenSimuCCBits.getBS().count() > 0 ) && (evenRawCCBits.getBS().count() > 0 ) && (oddSimuCCBits.getBS().count() > 0 ) && (oddRawCCBits.getBS().count() == 0 ) ){
	differenceCategories->Fill(14); // 1110 !
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 14 , 1110";
   } else if( (evenSimuCCBits.getBS().count() > 0 ) && (evenRawCCBits.getBS().count() > 0 ) && (oddSimuCCBits.getBS().count() > 0 ) && (oddRawCCBits.getBS().count() > 0 ) && (numberOfDifferentBitsEven + numberOfDifferentBitsOdd > 0)){
	differenceCategories->Fill(15); // 1111 !!
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " difference category 15 , 1111";
   }

   if(( oddSimuCCBits.getBS().count() * evenSimuCCBits.getBS().count()) > 0 ){
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " orientation category 0 , U and V";
	totalEfficiencySimu->Fill(0);
   } else if ( oddSimuCCBits.getBS().count() > 0 )  {
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " orientation category 0 , only U";
	totalEfficiencySimu->Fill(1);
   } else if ( evenSimuCCBits.getBS().count() > 0 )  {
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " orientation category 0 , only V";
	totalEfficiencySimu->Fill(2);
   } else if ( trigBitsNonZeroCounter > 0 ) {
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " orientation category 0 , no coincidence";
	totalEfficiencySimu->Fill(3);
   }


   /* 
    * *************************************
    * ***  GOOD TRACKS ANALYSIS  ***
    * *************************************
    */ 

   if( (fittableTrack >= 0) && (numberOfDifferentBitsEven + numberOfDifferentBitsOdd == 0) ){
	goodEventsCategories->Fill(0);
   } else if( (fittableTrack < 0) && (numberOfDifferentBitsEven + numberOfDifferentBitsOdd == 0) && (trigBitsNonZeroCounter > 0 ) ){
	goodEventsCategories->Fill(1);
   } else if( trigBitsNonZeroCounter == 0 ){
	goodEventsCategories->Fill(2);
   } else if( numberOfDifferentBitsEven + numberOfDifferentBitsOdd > 0 ){
        if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " wielka dupa !";
	goodEventsCategories->Fill(3);
   } 

   if( fittableTrack >= 0 ){

          /* 
           * *************************************
           * ******  ACTIVE CHIP ANALYSIS  *******
           * *************************************
           */ 

     if ((( oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) > 0 ) ){
	activeChipRaw->Fill(3);
     } else if ( (oddRawCCBits.getBS().count() == 0 )  && (evenRawCCBits.getBS().count() > 0 )) {
	activeChipRaw->Fill(2);
     } else if ( (oddRawCCBits.getBS().count() > 0 )  && (evenRawCCBits.getBS().count() == 0 )) {
	activeChipRaw->Fill(1);
     } else if ( (oddRawCCBits.getBS().count() == 0 )  && (evenRawCCBits.getBS().count() == 0 )) {
	activeChipRaw->Fill(0);
     }

     if ((( oddSimuCCBits.getBS().count() * evenSimuCCBits.getBS().count()) > 0 ) ){
	activeChipSim->Fill(3);
     } else if ( (oddSimuCCBits.getBS().count() == 0 )  && (evenSimuCCBits.getBS().count() > 0 )) {
	activeChipSim->Fill(2);
     } else if ( (oddSimuCCBits.getBS().count() > 0 )  && (evenSimuCCBits.getBS().count() == 0 )) {
	activeChipSim->Fill(1);
     } else if ( (oddSimuCCBits.getBS().count() == 0 )  && (evenSimuCCBits.getBS().count() == 0 )) {
	activeChipSim->Fill(0);
     }


     if ((( oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) > 0 ) &&  (numberOfDifferentBitsEven + numberOfDifferentBitsOdd == 0) ){
         if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " hardware : good tracks = yes , raw coincidence = yes = simu coinc";
         // good track and raw CC on
 	goodTracksCategories->Fill(0);
     } else if ((( oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) > 0 ) && (numberOfDifferentBitsEven + numberOfDifferentBitsOdd > 0) ) {
          if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " hardware : good tracks =  yes , raw coincidence = yes <> simu coinc";
         // no good track and raw CC on
  	goodTracksCategories->Fill(1);
     } else if ((( oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) == 0 ) && (numberOfDifferentBitsEven + numberOfDifferentBitsOdd > 0) ) {
          if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " hardware : good tracks = yes , raw coincidence = no, simu coinc = yes";
          // good track and raw CC off
  	goodTracksCategories->Fill(2);
     } else if ((( oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) == 0 ) &&  (numberOfDifferentBitsEven + numberOfDifferentBitsOdd == 0) ) {
          if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " hardware : good tracks =  yes , raw coincidence = no, simu coinc = no";
          // no good track and raw CC off
  	goodTracksCategories->Fill(3);
      }
   } else {

     if ((( oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) > 0 ) &&  (numberOfDifferentBitsEven + numberOfDifferentBitsOdd == 0) ){
         if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " hardware : good tracks = yes , coincidence = yes";
         // good track and raw CC on
 	goodTracksCategories->Fill(4);
     } else if ((( oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) > 0 ) && (numberOfDifferentBitsEven + numberOfDifferentBitsOdd > 0) ) {
          if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " hardware : good tracks =  no , coincidence = yes";
         // no good track and raw CC on
  	goodTracksCategories->Fill(5);
     } else if ((( oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) == 0 ) && (numberOfDifferentBitsEven + numberOfDifferentBitsOdd > 0) ) {
          if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " hardware : good tracks = yes , coincidence =  no";
          // good track and raw CC off
  	goodTracksCategories->Fill(6);
     } else if ((( oddRawCCBits.getBS().count() * evenRawCCBits.getBS().count()) == 0 ) &&  (numberOfDifferentBitsEven + numberOfDifferentBitsOdd == 0) ) {
          if( verbose_ > 2 ) LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer " << event.id() << " hardware : good tracks =  no , coincidence =  no";
          // no good track and raw CC off
  	goodTracksCategories->Fill(7);
      }

   }


}


//----------------------------------------------------------------------------------------------------

void RPCoincidenceAnalyzer::endJob()
{

  if( verbose_ ) edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer: endJob";	
  TFile * theFile = TFile::Open(outputFile.c_str(), "recreate");
  if(!theFile || !theFile->IsWritable())
  {
    std::cout << "Output file not opened correctly!!" << std::endl;
  }
  trigSectorsRawEven->Write();
  trigSectorsRawOdd->Write();
  trigSectorsSimuEven->Write();
  trigSectorsSimuOdd->Write();

  TCanvas * trigSectorsEvenCanvas = new TCanvas("TriggeredSectorsVdir","Triggered Sectors Raw vs Simu - V direction");
  trigSectorsEvenCanvas->Divide(2,1);
  trigSectorsEvenCanvas->cd(1);
  trigSectorsRawEven->Draw();
  trigSectorsEvenCanvas->cd(2);  
  trigSectorsSimuEven->Draw();
  trigSectorsEvenCanvas->Write();

  TCanvas * trigSectorsOddCanvas = new TCanvas("TriggeredSectorsUdir","Triggered Sectors Raw vs Simu - U direction");
  trigSectorsOddCanvas->Divide(2,1);
  trigSectorsOddCanvas->cd(1);
  trigSectorsRawOdd->Draw();
  trigSectorsOddCanvas->cd(2);  
  trigSectorsSimuOdd->Draw();
  trigSectorsOddCanvas->Write();

  differenceSimuVsRawOdd->Write();
  differenceSimuVsRawEven->Write();
  differenceSimuVsRawTotal->Write();

  TCanvas * totalEfficiencyCanvas = new TCanvas("Efficiency","Efficiency of simulated and real CC");
  totalEfficiencyCanvas->Divide(2,1);
  totalEfficiencyCanvas->cd(1);
  totalEfficiencyRaw->GetXaxis()->SetBinLabel(1,"U and V dir.");
  totalEfficiencyRaw->GetXaxis()->SetBinLabel(2,"U dir. only");
  totalEfficiencyRaw->GetXaxis()->SetBinLabel(3,"V dir. only");
  totalEfficiencyRaw->GetXaxis()->SetBinLabel(4,"no coincid.");
  totalEfficiencyRaw->GetXaxis()->SetBinLabel(5,"empty");
//  totalEfficiencyRaw->GetYaxis()->SetRangeUser(0.,1.);
  totalEfficiencyRaw->Draw();
  totalEfficiencyCanvas->cd(2);
  totalEfficiencySimu->GetXaxis()->SetBinLabel(1,"U and V dir.");
  totalEfficiencySimu->GetXaxis()->SetBinLabel(2,"U dir. only");
  totalEfficiencySimu->GetXaxis()->SetBinLabel(3,"V dir. only");
  totalEfficiencySimu->GetXaxis()->SetBinLabel(4,"no coincid.");
  totalEfficiencySimu->GetXaxis()->SetBinLabel(5,"empty");
//  totalEfficiencySimu->GetYaxis()->SetRangeUser(0.,1.);
  totalEfficiencySimu->Draw();
  totalEfficiencyCanvas->Write();


  TCanvas * differenceCategoriesCanvas = new TCanvas("DifferenceCategories","Difference Categories");
  differenceCategories->Draw();
  differenceCategoriesCanvas->Write();

  TCanvas * goodTracksCategoriesCanvas = new TCanvas("GoodTracksCategories","Good tracks Categories");
  goodTracksCategories->GetXaxis()->SetBinLabel(1,"Good track, Coinc");
  goodTracksCategories->GetXaxis()->SetBinLabel(2,"No track, Coinc");
  goodTracksCategories->GetXaxis()->SetBinLabel(3,"Good track, No coinc");
  //goodTracksCategories->GetXaxis()->SetBinLabel(4,"4 No track, No coinc");
  goodTracksCategories->Draw();
  goodTracksCategoriesCanvas->Write();

  TCanvas * tracksSelectionCategoriesCanvas = new TCanvas("TracksSelectionCatagories","Tracks Selection Catagories");
  tracksSelectionCategories->GetXaxis()->SetBinLabel(1,"Good track, fittable");
  tracksSelectionCategories->GetXaxis()->SetBinLabel(2,"Good track, nonfittable");
  tracksSelectionCategories->GetXaxis()->SetBinLabel(3,"No tracks, non-empty");
  tracksSelectionCategories->GetXaxis()->SetBinLabel(4,"No tracks, empty");
  tracksSelectionCategories->Draw();
  tracksSelectionCategoriesCanvas->Write();

  TCanvas * activeChipSimCanvas = new TCanvas("ActiveChipSim","Active Chip Simulation");
  activeChipSim->GetXaxis()->SetBinLabel(1,"Both CC active");
  activeChipSim->GetXaxis()->SetBinLabel(2,"Only odd CC (1,3..) active");
  activeChipSim->GetXaxis()->SetBinLabel(3,"Only even CC (0,2..) active");
  activeChipSim->GetXaxis()->SetBinLabel(4,"Both CC non-active");
  activeChipSim->Draw();
  activeChipSimCanvas->Write();

  TCanvas * activeChipRawCanvas = new TCanvas("ActiveChipRaw","Active Chip Real");
  activeChipRaw->GetXaxis()->SetBinLabel(1,"Both CC active");
  activeChipRaw->GetXaxis()->SetBinLabel(2,"Only odd CC (1,3..) active");
  activeChipRaw->GetXaxis()->SetBinLabel(3,"Only even CC (0,2..) active");
  activeChipRaw->GetXaxis()->SetBinLabel(4,"Both CC non-active");
  activeChipRaw->Draw();
  activeChipRawCanvas->Write();

  TCanvas * goodEventsCategoriesCanvas = new TCanvas("GoodEventsCategories","Good events Categories");
  goodEventsCategories->GetXaxis()->SetBinLabel(1,"Good track, non empty event, Sim CC = HW CC");
  goodEventsCategories->GetXaxis()->SetBinLabel(2,"No track, non empty event, Sim CC = HW CC");
  goodEventsCategories->GetXaxis()->SetBinLabel(3,"Empty event, Sim CC = HW CC");
  goodEventsCategories->GetXaxis()->SetBinLabel(4,"Sim CC != HW CC");
  goodEventsCategories->Draw();
  goodEventsCategoriesCanvas->Write();

  differenceSimuVsRawTotal->Print("range");

 // totalEfficiencyRaw->Print("range");
 // totalEfficiencySimu->Print("range");
//  differenceCategories->Print("range");
//  goodTracksCategories->Print("range");

  theFile->Close();


  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer ****** SUMMARY ********";
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer 1. Total number of events : " << tracksSelectionCategories->GetEntries();
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with good track (fittable): " << tracksSelectionCategories->GetBinContent(1);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with good track (nonfittable): " << tracksSelectionCategories->GetBinContent(2);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with no tracks, non-empty: " << tracksSelectionCategories->GetBinContent(3);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with no tracks, empty: " << tracksSelectionCategories->GetBinContent(4);

  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer 2. Total number of events : " << activeChipRaw->GetEntries();
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with even and odd real CC active: " << activeChipRaw->GetBinContent(1);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with only odd (1,3..) real CC active: " << activeChipRaw->GetBinContent(2);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with only even (0,2..) real CC active: " << activeChipRaw->GetBinContent(3);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with both real CC non-active: " << activeChipRaw->GetBinContent(4);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer      --------------------------------------------------------- ";
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with even and odd simulated CC active: " << activeChipSim->GetBinContent(1);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with only odd (1,3..) simulated CC active: " << activeChipSim->GetBinContent(2);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with only even (0,2..) simulated CC active: " << activeChipSim->GetBinContent(3);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with both simulated CC non-active: " << activeChipSim->GetBinContent(4);

  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer 3. Total number of events : " << goodEventsCategories->GetEntries();
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with good tracks, non-empty, for which SIM CC = HW CC: " << goodEventsCategories->GetBinContent(1);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with no tracks, non-empty, for which SIM CC = HW CC: " << goodEventsCategories->GetBinContent(2);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events empty, for which SIM CC = HW CC: " << goodEventsCategories->GetBinContent(3);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM CC != HW CC: " << goodEventsCategories->GetBinContent(4);

  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer 4. Total number of events : " << differenceSimuVsRawTotal->GetEntries();
  for( int i = 1 ; i <= 16 ; i++ ){
	if( differenceSimuVsRawTotal->GetBinContent(i) > 0 )
  		edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events with " << (i-1) << " different bits between SIM CC output and HW CC output:" << differenceSimuVsRawTotal->GetBinContent(i);
  }

  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer 5. Total number of events : " << differenceCategories->GetEntries();
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC == 0, HW Even CC == 0, Sim Odd CC == 0, HW Odd CC > 0: " << differenceCategories->GetBinContent(2);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC == 0, HW Even CC == 0, Sim Odd CC > 0, HW Odd CC == 0: " << differenceCategories->GetBinContent(3);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC == 0, HW Even CC == 0, Sim Odd CC > 0, HW Odd CC > 0 (diff): " << differenceCategories->GetBinContent(4);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC == 0, HW Even CC > 0, Sim Odd CC == 0, HW Odd CC == 0: " << differenceCategories->GetBinContent(5);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC == 0, HW Even CC > 0, Sim Odd CC == 0, HW Odd CC > 0: " << differenceCategories->GetBinContent(6);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC == 0, HW Even CC > 0, Sim Odd CC > 0, HW Odd CC == 0: " << differenceCategories->GetBinContent(7);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC == 0, HW Even CC > 0, Sim Odd CC > 0, HW Odd CC > 0 (diff): " << differenceCategories->GetBinContent(8);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC > 0, HW Even CC == 0, Sim Odd CC == 0, HW Odd CC == 0: " << differenceCategories->GetBinContent(9);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC > 0, HW Even CC == 0, Sim Odd CC == 0, HW Odd CC > 0: " << differenceCategories->GetBinContent(10);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC > 0, HW Even CC == 0, Sim Odd CC > 0, HW Odd CC == 0: " << differenceCategories->GetBinContent(11);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC > 0, HW Even CC == 0, Sim Odd CC > 0, HW Odd CC > 0 (diff): " << differenceCategories->GetBinContent(12);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC > 0, HW Even CC > 0 (diff), Sim Odd CC == 0, HW Odd CC == 0: " << differenceCategories->GetBinContent(13);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC > 0, HW Even CC > 0 (diff), Sim Odd CC == 0, HW Odd CC > 0: " << differenceCategories->GetBinContent(14);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC > 0, HW Even CC > 0 (diff), Sim Odd CC > 0, HW Odd CC == 0: " << differenceCategories->GetBinContent(15);
  edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer    Number of events for which SIM Even CC > 0, HW Even CC > 0 (diff), Sim Odd CC > 0, HW Odd CC > 0 (diff): " << differenceCategories->GetBinContent(16);
// if( verbose_ > 2 ) edm::LogPrint("RPCoincidenceAnalyzer") << "RPCoincidenceAnalyzer Total number of events : " << differenceSimuVsRawTotal->GetSize();

}


DEFINE_FWK_MODULE(RPCoincidenceAnalyzer);
