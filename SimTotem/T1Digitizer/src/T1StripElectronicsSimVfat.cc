/*
 * modified by Marcin Borratynski (mborratynski@gmail.com)
 */
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimTotem/T1Digitizer/interface/T1StripElectronicsSimVfat.h"
#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"
//#include "SimTotem/T1Digitizer/src/CSCAnalogSignal.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfat.h"
//#include "Geometry/CSCGeometry/interface/CSCLayer.h"
//#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "Geometry/TotemGeometry/interface/T1ChamberSpecs.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandGaussQ.h"
#include <iostream>

//#define THR_1 500//30
//#define THR_2 1000//50

T1StripElectronicsSimVfat::T1StripElectronicsSimVfat(const edm::ParameterSet & p):  _threshold1(p.getParameter<double>("THR1")),
										    _threshold2(p.getParameter<double>("THR2")),
										    _rumore(p.getParameter<double>("NOISE"))

{

  theStripThreshold1 = _threshold1;
  theStripThreshold2 = _threshold2;

  theStripNoise = _rumore;
  int verbo = p.getParameter<int>("Verbosity");
  if(verbo>=1){
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@ THR1="<<theStripThreshold1 << std::endl;                 
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@ THR2="<<theStripThreshold2 << std::endl;  
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@ NOISE="<<theStripNoise << std::endl; }
}


 
void T1StripElectronicsSimVfat::fillDigis(int event, int detid,
		std::vector<T1DetectorHit> stripHits, T1DigiVfatCollection & digis, T1DeadChannelManager & channelManager) {


  //  ATTENZIONE AL NUMERO DI STRIP NEGATIVO O POSITIVO A SECONDA DEL PIANO A O B


  T1DetId myDetid(detid);

  int Arm = myDetid.Arm();
  int Plane = myDetid.Plane();
  int CSC = myDetid.CSC();
  

  float stripChargeA[226];
  float stripChargeB[226];
  for(int i=0; i<226; i++){ stripChargeA[i]=0; stripChargeB[i]=0; }

  float threshold1 = std::max(theStripThreshold1 + CLHEP::RandGaussQ::shoot() * theStripNoise,0.);
  float threshold2 = std::max(theStripThreshold2 + CLHEP::RandGaussQ::shoot() * theStripNoise,0.);




  std::vector<T1DetectorHit>::iterator it_WH;

  // ATTENZIONE: PER LE STRISCE NEGATIVE OCCORRE CONSIDERARE UN ID DIFFERENTE

  for(it_WH = stripHits.begin(); it_WH != stripHits.end(); it_WH++){

    if( (*it_WH).getElement() >0 )
      stripChargeA[(*it_WH).getElement()] = stripChargeA[(*it_WH).getElement()] + (*it_WH).getCharge();
    if ( (*it_WH).getElement() <0 )
      stripChargeB[-(*it_WH).getElement()] = stripChargeB[-(*it_WH).getElement()] + (*it_WH).getCharge();

  }

  float const igain = 53.; // mV/fC

  for(int stripNumber=1; stripNumber<225; stripNumber++){


	// piano A
	if (stripChargeA[stripNumber] * igain > threshold1	|| stripChargeA[stripNumber] * igain > threshold2) {

		myDetid.setLayer(Arm, Plane, CSC, 1); // strips plane A

		if(!channelManager.isChannelDead(myDetid, stripNumber)){
			T1DigiVfat sdigi(stripNumber, 0, event, stripChargeA[stripNumber]);

			if (stripChargeA[stripNumber] * igain > threshold1	&& stripChargeA[stripNumber] * igain < threshold2)
				sdigi.setThreshold(1);
			if (stripChargeA[stripNumber] * igain >= threshold2)
				sdigi.setThreshold(2);

			digis.insertDigi(myDetid, sdigi);
		}
	}

	// piano B
	if (stripChargeB[stripNumber] * igain > threshold1	|| stripChargeB[stripNumber] * igain > threshold2) {

		myDetid.setLayer(Arm, Plane, CSC, 2); // strips plane B
		if(!channelManager.isChannelDead(myDetid, stripNumber)){
			T1DigiVfat sdigi(stripNumber, 0, event, stripChargeB[stripNumber]);
			if (stripChargeB[stripNumber] * igain > threshold1 	&& stripChargeB[stripNumber] * igain < threshold2)
				sdigi.setThreshold(1);
			if (stripChargeB[stripNumber] * igain >= threshold2)
				sdigi.setThreshold(2);

			digis.insertDigi(myDetid, sdigi);
		}
	}
  }
}




