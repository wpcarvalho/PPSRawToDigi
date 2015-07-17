#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimTotem/T1Digitizer/interface/T1WireElectronicsSim.h"
#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"
//#include "SimTotem/T1Digitizer/src/CSCAnalogSignal.h"
#include "DataFormats/T1DigiWire/interface/T1DigiWire.h"
//#include "Geometry/CSCGeometry/interface/CSCLayer.h"
//#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "Geometry/TotemGeometry/interface/T1ChamberSpecs.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandGaussQ.h"
#include <iostream>



T1WireElectronicsSim::T1WireElectronicsSim(const edm::ParameterSet & p) 
  :  _threshold(p.getParameter<double>("WIRETHR")),
     _rumore(p.getParameter<double>("WIRENOISE"))
{

  theWireThreshold = _threshold;
  theWireNoise = _rumore;
  int verbo = p.getParameter<int>("Verbosity");
  if(verbo>=1){
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@ WIRE THR="<<_threshold << std::endl; 
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@ WIRE NOISE="<< _rumore << std::endl; 
  }
}
/*
  T1WireElectronicsSim::T1WireElectronicsSim() 
 
  {
  theWireThreshold = 50.;
  theWireNoise = 0;
  }
*/


/*
  CSCWireElectronicsSim::CSCWireElectronicsSim(const edm::ParameterSet & p) 
  : CSCBaseElectronicsSim()
  {
  theSignalStartTime = p.getParameter<double>("wireSignalStartTime");
  theSignalStopTime = p.getParameter<double>("wireSignalStopTime");
  theSamplingTime = p.getParameter<double>("wireSamplingTime");
  theTimingCalibrationError = p.getParameter<double>("wireTimingError");
  init();
  }
*/

/*
  CSCWireElectronicsSim::CSCWireElectronicsSim()
  : CSCBaseElectronicsSim()
  {
  theSignalStartTime = -250.;
  theSignalStopTime = 500.;
  theSamplingTime = 2.0;
  theTimingCalibrationError = 0.;
  init();
  }


  void T1WireElectronicsSim::init() {
  theFraction = 0.5;
  theShapingTime = 30;
  theAmpGainVariance = 0.;
  thePeakTimeVariance = 0.;
  theTailShaping = RADICAL;

  theBunchTimingOffsets.resize(11);
  theBunchTimingOffsets[1] = 33.8;
  theBunchTimingOffsets[2] = 34.8;
  theBunchTimingOffsets[3] = 42.0;
  theBunchTimingOffsets[4] = 42.0;
  theBunchTimingOffsets[5] = 43.6;
  theBunchTimingOffsets[6] = 42.0;
  theBunchTimingOffsets[7] = 42.6;
  theBunchTimingOffsets[8] = 41.5;
  theBunchTimingOffsets[9] = 43.6;
  theBunchTimingOffsets[10] = 41.5;
  theNumberOfSamples = (int)((theSignalStopTime-theSignalStartTime)/theSamplingTime);
  fillAmpResponse();
  }


  void T1WireElectronicsSim::initParameters() {
  theLayerGeometry = theLayer->geometry();
  nElements = theLayerGeometry->numberOfWireGroups();
  theWireNoise = theSpecs->wireNoise(theShapingTime)
  * e_SI * pow(10.0,15);
  theWireThreshold = theWireNoise * 8;

  }


  int T1WireElectronicsSim::readoutElement(int element) const {
  return theLayerGeometry->wireGroup(element);
  }

*/

void T1WireElectronicsSim::fillDigis(int event, int detid,
		std::vector<T1DetectorHit> wireHits, T1DigiWireCollection & digis, T1DeadChannelManager & channelManager) {


  T1DetId myDetid(detid);
  

  int Arm = myDetid.Arm();
  int Plane = myDetid.Plane();
  int CSC = myDetid.CSC();
 
  myDetid.setLayer(Arm,Plane,CSC,3); // set layer = 3 => wires


  float wireCharge[226];
  for(int i=0; i<226; i++)wireCharge[i]=0;

  float threshold = theWireThreshold + CLHEP::RandGaussQ::shoot() * theWireNoise;




  std::vector<T1DetectorHit>::iterator it_WH;



  for(it_WH = wireHits.begin(); it_WH != wireHits.end(); it_WH++){

    wireCharge[(*it_WH).getElement()] = wireCharge[(*it_WH).getElement()] + (*it_WH).getCharge();
     

  }

  float const igain = 53.; // mV/fC

  for(int wireNumber=1; wireNumber<226; wireNumber++){
	if (wireCharge[wireNumber] * igain > threshold) {
		/**
		 * check if channel is not in masked on list of dead/noisy channels
		 */
		if(!channelManager.isChannelDead(myDetid, wireNumber)){

			//       std::cout << " Wire " << wireNumber << "   Charge " << wireCharge[wireNumber]*igain << std::endl;
			T1DigiWire wdigi(wireNumber, event, wireCharge[wireNumber]);
			digis.insertDigi(myDetid, wdigi);
		}
	}
  }
}




