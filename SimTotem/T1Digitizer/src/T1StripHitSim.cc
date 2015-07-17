#include "SimTotem/T1Digitizer/interface/T1StripHitSim.h"
#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include <cmath>
#include <iostream>
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#define COS30 0.866025


//#define _DEBUG_

std::vector<T1DetectorHit> &
T1StripHitSim::simulate(uint32_t myDetId_, 
			const std::vector<T1DetectorHit> & wireHits, T1Geometry *layer)
{
  // make sure the gatti function is initialized
  //  T1Geometry * layer = new T1Geometry();
  const T1ChamberSpecs chamberSpecs;


  //%//  std::cout << " CALLING GATTI FUNCTION IN T1StripHitSim::simulate   ----------- "<<std::endl;

  initChamberSpecs();
  const int nNodes = chamberSpecs.nNodes();
  
  // for every wire hit, induce a Gatti-distributed charge on the
  // cathode strips
  newStripHits.clear();
  newStripHits.reserve((2*nNodes+1)*wireHits.size());
  std::vector<T1DetectorHit>::const_iterator wireHitI;
#ifdef _DEBUG_
  int flagg = 0;
#endif

  for(wireHitI = wireHits.begin(); wireHitI != wireHits.end(); ++wireHitI){
    int   wire        = (*wireHitI).getElement();
    float wireCharge  = (*wireHitI).getCharge();
    float wireHitTime = (*wireHitI).getTime();


    // The wire hit position is _along the wire_, measured from where
    // the wire intersects local y axis, so convert to local x...
    //    float hitX   = (*wireHitI).getPosition() * cos(geom->wireAngle());

  
#ifdef _DEBUG_
    if(flagg==0){
      T1DetId tempId(myDetId_);
      std::cout << tempId << std::endl;
      std::cout << " hitX   = (*wireHitI).getPosition() + layer->xOfWire(myDetId_) " << (*wireHitI).getPosition() << " + " << layer->xOfWire(myDetId_)<<std::endl;
    
    }
#endif
    //    float hitX   = (*wireHitI).getPosition() - layer->xOfWire(myDetId_);
    float hitX   = (*wireHitI).getPosition() + layer->xOfWire(myDetId_);
    float hitY   = layer->yOfWire(myDetId_,wire);

    // wireHitPos: posizione nel sistema di riferimento di geant4

    //   const LocalPoint wireHitPos(hitX, hitY);
    /*
      float discrepancyA = 0.083;
      float discrepancyB = 0.11;
    */
    //    std::cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    //  std::cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;

    int centerStripA = (int)layer->tellStripA_TB(myDetId_,hitX,hitY);
    int firstStripA = std::max(centerStripA - nNodes, 1);
    int lastStripA  = std::min(centerStripA + nNodes, layer->numberOfStrips(myDetId_));
    int centerStripB = (int)layer->tellStripB_TB(myDetId_,hitX,hitY);
    int firstStripB = std::max(centerStripB - nNodes, 1);
    int lastStripB  = std::min(centerStripB + nNodes, layer->numberOfStrips(myDetId_));

#ifdef _DEBUG_
    if(flagg==0){
      std::cout <<"Center strip A "<< centerStripA << "   Center strip B " << centerStripB << "  Nodes " << nNodes<< " of "<<layer->numberOfStrips(myDetId_) << " Det ID: "<<myDetId_ << std::endl;

      //    int lastStrip  = std::min(centerStrip + nNodes, geom->numberOfStrips());
 
      flagg=1;
    }
    //    std::cout << " -------------------------------------------------------------------------  STRIPS A "<<std::endl;
#endif

    for(int istrip = firstStripA; istrip <= lastStripA; istrip++) {

      float XofstripA = layer->xOfStripA(myDetId_,istrip, hitY);
      //      float offset = (hitX - XofstripA)*COS30 - discrepancyA;
      float offset = (hitX - XofstripA)*COS30;
      float stripWidth = chamberSpecs.Pitch();
      //      float binValue = theGattiFunction.binValue(offset, stripWidth);
      float binValue_ = binValue(offset, stripWidth);

#ifdef _DEBUG_
      //      std::cout << " ISTRIP E OFFSET : " << istrip << " " << offset << " "<< hitX <<" "<< XofstripA<<std::endl;
#endif
      // we divide by 2 because charge goes on either cathode.
      // if you're following the TDR, we already multiplied the
      // charge by 0.82 in the WireHitSim (well, DriftSim), so that explains 
      // their f_ind=0.41.
   
      // this seems to be folded in the Amp response, which peaks
      // around 0.14.  The difference may be because the amp response
      // convolutes in different drift times.
      //float collectionFraction = 0.19;
      //abbiamo messo f_ind = 0.43 quindi probabilmente 
      //non occorre dividere la carica per 2 
      //   const float igain = 1./0.9; // mv/fC        ?????????????????????????????

      const float igain = 1.; 
      static const double f_ind = 0.43;

      float stripCharge = wireCharge * binValue_ * igain * f_ind;
#ifdef _DEBUG_
      //    std::cout << " stripCharge = wireCharge * binValue_ * igain "  << stripCharge << " "<< wireCharge << " "<< binValue_ << " "<< igain <<std::endl;
#endif
      float stripTime = wireHitTime;
      float position = 0.0;

      //      float position = hitY / sin(geom->stripAngle(istrip));
      T1DetectorHit newStripHit(istrip, stripCharge, position, stripTime, 
				(*wireHitI).getSimHit(),myDetId_);
      newStripHits.push_back(newStripHit);
    }


#ifdef _DEBUG_
    //   std::cout << " -------------------------------------------------------------------------  STRIPS B "<<std::endl;
#endif

    for(int istrip = firstStripB; istrip <= lastStripB; istrip++) {

      float XofstripB = layer->xOfStripB(myDetId_,istrip, hitY);
      // ?????????????????????????  strips side A or B

      //      float offset = (hitX - XofstripB)*COS30 + discrepancyB;
      float offset = (hitX - XofstripB)*COS30;
    
      float stripWidth = chamberSpecs.Pitch();
      //      float binValue = theGattiFunction.binValue(offset, stripWidth);
      float binValue_ = binValue(offset, stripWidth);
#ifdef _DEBUG_
      //    std::cout << " ISTRIP E OFFSET : " << istrip << " " << offset << " "<< hitX <<" "<< XofstripB<<std::endl;
#endif

      // we divide by 2 because charge goes on either cathode.
      // if you're following the TDR, we already multiplied the
      // charge by 0.82 in the WireHitSim (well, DriftSim), so that explains 
      // their f_ind=0.41.
   
      // this seems to be folded in the Amp response, which peaks
      // around 0.14.  The difference may be because the amp response
      // convolutes in different drift times.
      //float collectionFraction = 0.19;
      //abbiamo messo f_ind = 0.43 quindi probabilmente 
      //non occorre dividere la carica per 2 
      //   const float igain = 1./0.9; // mv/fC        ?????????????????????????????

      const float igain = 1.; // mv/fC  potrebbe essere 53 per il VFAT - lo mettiamo in electronics sim
      static const double f_ind = 0.43;

      float stripCharge = wireCharge * binValue_ * igain * f_ind;
#ifdef _DEBUG_
      //      std::cout << " stripCharge = wireCharge * binValue_ * igain "  << stripCharge << " "<< wireCharge << " "<< binValue_ << " "<< igain <<std::endl;
#endif

      float stripTime = wireHitTime;
      float position = 0.0;

      //      float position = hitY / sin(geom->stripAngle(istrip));

      // !!!  strip side B with negative number, in order to put them in the same collection as strip A and distinguish them  !!!

      T1DetectorHit newStripHit(-istrip, stripCharge, position, stripTime, 
				(*wireHitI).getSimHit(),myDetId_);
      newStripHits.push_back(newStripHit);
    }



  }  // loop over wire hits



  //  std::cout << " Fine di T1StripHitSim::simulate -------------- "<<std::endl;  

  return newStripHits;
}

void T1StripHitSim::initChamberSpecs(){
  const T1ChamberSpecs chamberSpecs;

  //  std::cout << "T1GattiFunction::initChamberSpecs setting new values." <<std::endl;
  h = chamberSpecs.anodeCathodeSpacing();
  double s = chamberSpecs.wireSpacing();
  double ra = chamberSpecs.wireRadius();
  static const double parm[5] = {.1989337e-02, -.6901542e-04,  .8665786, 
				 154.6177, -.6801630e-03 };
  k3 = (parm[0]*s/h + parm[1]) 
    * (parm[2]*s/ra + parm[3] + parm[4]*s*s/ra/ra);
  sqrtk3 = sqrt(k3);
  norm = 0.5 / std::atan( sqrtk3 );
  k2 = M_PI_2 * (1. - sqrtk3/2.);
  k1 = 0.25 * k2 * sqrtk3 / std::atan(sqrtk3);
    

  //    std::cout  << "Gatti function constants k1=" <<  k1 << ", k2=" << k2 << ", k3=" << k3 <<  ", h=" << h << ", norm=" << norm<<std::endl;

  return;
}


double T1StripHitSim::binValue( double x, double stripWidth) const {

  double tanh1 = tanh(k2 * (x+stripWidth*0.5)/h );
  double tanh2 = tanh(k2 * (x-stripWidth*0.5)/h );
  return norm * ( std::atan(sqrtk3*tanh1) - std::atan(sqrtk3*tanh2) );
}


