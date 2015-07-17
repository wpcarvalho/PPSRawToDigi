/*
 * modified by Marcin Borratynski (mborratynski@gmail.com)
 */
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimTotem/T1Digitizer/interface/T1StripElectronicsSimSantiard.h"
#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"
#include "DataFormats/T1DigiSantiard/interface/T1DigiSantiard.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "Geometry/TotemGeometry/interface/T1ChamberSpecs.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandGaussQ.h"
#include <iostream>

//#define _DEBUG_



T1StripElectronicsSimSantiard::T1StripElectronicsSimSantiard(const edm::ParameterSet & p) :  _threshold(p.getParameter<double>("THR")),
											     _rumore(p.getParameter<double>("NOISE"))
{

  theStripThreshold = _threshold;
  theStripNoise = _rumore;
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@ THR = "<<theStripThreshold << std::endl;                 
 
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@ NOISE="<<theStripNoise << std::endl;  
}




 
void T1StripElectronicsSimSantiard::fillDigis(int event, int detid,
		std::vector<T1DetectorHit> stripHits, T1DigiSantiardCollection & digis, T1DeadChannelManager & channelManager) {


  //  ATTENZIONE AL NUMERO DI STRIP NEGATIVO O POSITIVO A SECONDA DEL PIANO A O B


  T1DetId myDetid(detid);

  int Arm = myDetid.Arm();
  int Plane = myDetid.Plane();
  int CSC = myDetid.CSC();
  

  float stripChargeA[226];
  float stripChargeB[226];
  for(int i=0; i<226; i++){ stripChargeA[i]=0; stripChargeB[i]=0; }

  float threshold = theStripThreshold + CLHEP::RandGaussQ::shoot() * theStripNoise;




  std::vector<T1DetectorHit>::iterator it_WH;

  // ATTENZIONE: PER LE STRISCE NEGATIVE OCCORRE CONSIDERARE UN ID DIFFERENTE

  for(it_WH = stripHits.begin(); it_WH != stripHits.end(); it_WH++){

    if( (*it_WH).getElement() >0 )
      stripChargeA[(*it_WH).getElement()] = stripChargeA[(*it_WH).getElement()] + (*it_WH).getCharge();
    if ( (*it_WH).getElement() <0 )
      stripChargeB[-(*it_WH).getElement()] = stripChargeB[-(*it_WH).getElement()] + (*it_WH).getCharge();

  }

  float const igain = 1.; // mV/fC

  int balance = -1;

  for(int stripNumber=1; stripNumber<225; stripNumber++){


    // comparazione a tre a tre piano A
#ifdef _DEBUG_
    if(stripChargeA[stripNumber]>0) std::cout << " carica (" <<Plane<<","<<CSC<<") "<< stripNumber << " " <<stripChargeA[stripNumber]<<std::endl;
#endif

    if (stripChargeA[stripNumber] * igain > threshold
			&& stripChargeA[stripNumber] > stripChargeA[stripNumber - 1]
			&& stripChargeA[stripNumber] > stripChargeA[stripNumber + 1]) {
		//peak found
		if (stripChargeA[stripNumber - 1] >= stripChargeA[stripNumber + 1]) {
			balance = 0;
		} else {
			balance = 1;
		}

		myDetid.setLayer(Arm, Plane, CSC, 1); // strips plane A

		//check if channel is masked or not
		if(!channelManager.isChannelDead(myDetid, stripNumber)){

			// se ho un picco creo il digi
			T1DigiSantiard sdigi(stripNumber, balance, event, stripChargeA[stripNumber]);
			digis.insertDigi(myDetid, sdigi);
		}
	}

	// comparazione a tre a tre piano B

	if (stripChargeB[stripNumber] * igain > threshold
			&& stripChargeB[stripNumber] > stripChargeB[stripNumber - 1]
			&& stripChargeB[stripNumber] > stripChargeB[stripNumber + 1]) {
		//peak found
		if (stripChargeB[stripNumber - 1] >= stripChargeB[stripNumber + 1]) {
			balance = 0;
		} else {
			balance = 1;
		}

		myDetid.setLayer(Arm, Plane, CSC, 2); // strips plane B

		//check if channel is masked or not
		if(!channelManager.isChannelDead(myDetid, stripNumber)){

			// se ho un picco creo il digi
			T1DigiSantiard sdigi(stripNumber, balance, event, stripChargeB[stripNumber]);
			digis.insertDigi(myDetid, sdigi);
		}
	}
  }

}




