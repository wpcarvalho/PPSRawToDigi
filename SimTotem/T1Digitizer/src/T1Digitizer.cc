/*
 * modified by Marcin Borratynski (mborratynski@gmail.com)
 */
#include "Utilities/Timing/interface/TimingReport.h" 
#include "SimTotem/T1Digitizer/interface/T1Digitizer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"
#include "SimTotem/T1Digitizer/interface/T1WireHitSim.h"
#include "SimTotem/T1Digitizer/interface/T1StripHitSim.h"
#include "SimTotem/T1Digitizer/interface/T1DriftSim.h"

#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>

// #define _DEBUG_
// #define _LOGINFO_

T1Digitizer::T1Digitizer(const edm::ParameterSet & p) {

  //  std::cout << " Creating the T1Digitizer --------------------- "<<std::endl; 
  theDriftSim = new T1DriftSim();
  theWireHitSim          = new T1WireHitSim(theDriftSim);
  theStripHitSim         = new T1StripHitSim();
  std::string thegeometryfile = "Geometry/TotemGeometry/data/T1_data_geometry.dat";
  theT1Geometry = new T1Geometry(thegeometryfile);
  theWireElectronicsSim  = new T1WireElectronicsSim(p);
  if( p.getParameter<std::string>("Electronics") == "Santiard"){
    theStripElectronicsSimSantiard = new T1StripElectronicsSimSantiard(p);
    theStripElectronicsSimVfat = 0;
  }
  else if( p.getParameter<std::string>("Electronics") == "VFAT" ){
    theStripElectronicsSimVfat = new T1StripElectronicsSimVfat(p);
    theStripElectronicsSimSantiard = 0;
  } else{
    throw cms::Exception("")<<"Electronics for cathodes not set or wrong!";
  }
}


T1Digitizer::~T1Digitizer() {
  delete theT1Geometry;
  if(theStripElectronicsSimSantiard)
    delete theStripElectronicsSimSantiard;
  if( theStripElectronicsSimVfat)
    delete theStripElectronicsSimVfat;
  delete theWireElectronicsSim;
  delete theStripHitSim;
  delete theWireHitSim;
  delete theDriftSim;
}



void T1Digitizer::doAction(
			   int event,//const std::vector<PSimHit> * simHits, 
			   MixCollection<PSimHit> & simHits, 
			   T1DigiWireCollection & wireDigis, 
			   T1DigiSantiardCollection & stripDigis,
			   T1DeadChannelManager & channelManager)
{

  // arrange the hits by layer
  std::map<int, edm::PSimHitContainer> hitMap;
  //  for(std::vector<PSimHit>::const_iterator hitItr = simHits->begin();
  //     hitItr != simHits->end(); ++hitItr) 
  for(MixCollection<PSimHit>::MixItr hitItr = simHits.begin();
      hitItr != simHits.end(); ++hitItr) 
    {
      hitMap[hitItr->detUnitId()].push_back(*hitItr);
    }

  // now loop over layers and run the simulation for each one
  for(std::map<int, edm::PSimHitContainer>::const_iterator hitMapItr = hitMap.begin();
      hitMapItr != hitMap.end(); ++hitMapItr)
    {
      uint32_t IDD=hitMapItr->first;

      const edm::PSimHitContainer & mySimHits = hitMapItr->second;

      std::vector<T1DetectorHit> newWireHits, newStripHits;
 
#ifdef _LOGINFO_
      edm::LogInfo("T1Digitizer") << "T1Digitizer: found " << mySimHits.size() <<" hit(s) in layer";
#endif
      // turn the edm::PSimHits into WireHits, using the WireHitSim
      {
	TimeMe t("CSCWireHitSim");
      
	newWireHits.swap(theWireHitSim->simulate(IDD, mySimHits,theT1Geometry));
      
	std::vector<T1DetectorHit>::iterator it_DH;
      
	float hitCharge=0;
	for(it_DH = newWireHits.begin(); it_DH != newWireHits.end(); it_DH++){
	  hitCharge=hitCharge+(*it_DH).getCharge();
	}
#ifdef _DEBUG_
	std::cout << " HIT CHARGE " << hitCharge <<"   ++++++++++++++++++++++++++++++++++++++++++++++++++++++ " <<std::endl;
#endif
      }
      if(!newWireHits.empty()) {
	TimeMe t("CSCStripHitSim");
	newStripHits.swap(theStripHitSim->simulate(IDD, newWireHits, theT1Geometry));
      }
      // turn the hits into wire digis, using the electronicsSim

      {
	TimeMe t("T1WireElectronicsSim");
	theWireElectronicsSim->fillDigis(event,IDD,newWireHits,wireDigis, channelManager);
      }  

      {
	TimeMe t("CSCStripElectronicsSim");
	theStripElectronicsSimSantiard->fillDigis(event,IDD,newStripHits, stripDigis, channelManager);
      }
    }



}




void T1Digitizer::doAction(
			   int event,//const std::vector<PSimHit> * simHits, 
			   MixCollection<PSimHit> & simHits, 
			   T1DigiWireCollection & wireDigis, 
			   T1DigiVfatCollection & stripDigis,
			   T1DeadChannelManager & channelManager)
{

  // arrange the hits by layer
  std::map<int, edm::PSimHitContainer> hitMap;
  //  for(std::vector<PSimHit>::const_iterator hitItr = simHits->begin();
  //     hitItr != simHits->end(); ++hitItr) 
  for(MixCollection<PSimHit>::MixItr hitItr = simHits.begin();
      hitItr != simHits.end(); ++hitItr) 
    {
      hitMap[hitItr->detUnitId()].push_back(*hitItr);
    }

  // now loop over layers and run the simulation for each one
  for(std::map<int, edm::PSimHitContainer>::const_iterator hitMapItr = hitMap.begin();
      hitMapItr != hitMap.end(); ++hitMapItr)
    {
      uint32_t IDD=hitMapItr->first;

      const edm::PSimHitContainer & mySimHits = hitMapItr->second;

      std::vector<T1DetectorHit> newWireHits, newStripHits;

#ifdef _LOGINFO_
      edm::LogInfo("T1Digitizer") << "T1Digitizer: found " << mySimHits.size() <<" hit(s) in layer";
#endif
      // turn the edm::PSimHits into WireHits, using the WireHitSim
      {
	TimeMe t("CSCWireHitSim");
      
	newWireHits.swap(theWireHitSim->simulate(IDD, mySimHits,theT1Geometry));
      
	std::vector<T1DetectorHit>::iterator it_DH;
      
	float hitCharge=0;
	for(it_DH = newWireHits.begin(); it_DH != newWireHits.end(); it_DH++){
	  hitCharge=hitCharge+(*it_DH).getCharge();
	}
#ifdef _DEBUG_
	std::cout << " HIT CHARGE " << hitCharge <<"   ++++++++++++++++++++++++++++++++++++++++++++++++++++++ " <<std::endl;
#endif
      }
      if(!newWireHits.empty()) {
	TimeMe t("CSCStripHitSim");
	newStripHits.swap(theStripHitSim->simulate(IDD, newWireHits, theT1Geometry));
      }

      // turn the hits into wire digis, using the electronicsSim
      {
	TimeMe t("T1WireElectronicsSim");
	theWireElectronicsSim->fillDigis(event,IDD,newWireHits,wireDigis,channelManager);
      }  

      {
	theStripElectronicsSimVfat->fillDigis(event,IDD,newStripHits, stripDigis,channelManager);
      }
    }

}


void T1Digitizer::setMagneticField(const MagneticField * field) {
  theDriftSim->setMagneticField(field);
}
