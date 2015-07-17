/*
 * modified by Marcin Borratynski (mborratynski@gmail.com)
 */
#ifndef T1___DIGITIZER
#define T1___DIGITIZER

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"
#include "DataFormats/T1DigiSantiard/interface/T1DigiSantiardCollection.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfatCollection.h"

//#include "DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
//#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimTotem/T1Digitizer/interface/T1DriftSim.h"
#include "SimTotem/T1Digitizer/interface/T1WireHitSim.h"
#include "SimTotem/T1Digitizer/interface/T1StripHitSim.h"
#include "SimTotem/T1Digitizer/interface/T1WireElectronicsSim.h"
#include "SimTotem/T1Digitizer/interface/T1StripElectronicsSimSantiard.h"
#include "SimTotem/T1Digitizer/interface/T1StripElectronicsSimVfat.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimTotem/T1Digitizer/interface/T1DeadChannelManager.h"
#include "TotemCondFormats/DAQInformation/interface/AnalysisMask.h"
#include "TotemCondFormats/DataRecord/interface/TotemDAQMappingRecord.h"
//class T1DriftSim;
//class T1WireHitSim;
//class T1StripHitSim;
//class T1WireElectronicsSim;
//class T1StripElectronicsSim;
//class T1Geometry;
//class T1NeutronFactory;


class T1Digitizer 
{
 public:
  /// configurable parameters
  explicit T1Digitizer(const edm::ParameterSet & p);
  
  ~T1Digitizer();

  /**  digitize
   */
  void doAction(int evt, MixCollection<PSimHit> & simHits,
                T1DigiWireCollection & wireDigis,
                T1DigiSantiardCollection & stripDigis,
                T1DeadChannelManager & channelManager
		);
  void doAction(int evt,
		MixCollection<PSimHit> & simHits,
                T1DigiWireCollection & wireDigis,
                T1DigiVfatCollection & stripDigis,
                T1DeadChannelManager & channelManager
	   
		);

  /// sets the magnetic field
  void setMagneticField(const MagneticField * field);

 private:
  T1DriftSim            * theDriftSim;
  T1WireHitSim          * theWireHitSim;
  T1StripHitSim         * theStripHitSim;


  T1WireElectronicsSim  * theWireElectronicsSim;
  T1StripElectronicsSimSantiard * theStripElectronicsSimSantiard;
  T1StripElectronicsSimVfat * theStripElectronicsSimVfat;
  T1Geometry      * theT1Geometry;
};

#endif

