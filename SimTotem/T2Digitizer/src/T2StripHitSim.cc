/** 
 * Class simulates hit on a strip in the T2 detector. Polynom that calculates
 * the charge distribution has been developed by Eraldo Olivieri
 *
 * Author: Erik Brücken / University of Helsinki
 * Email:  brucken@cc.helsinki.fi
 * Date:   2008-02-22
 */

#include "SimTotem/T2Digitizer/interface/T2StripHitSim.h"
#include "SimTotem/T2Digitizer/interface/T2DetectorHit.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/TotemGeometry/interface/T2Geometry.h"
#include "DataFormats/T2Cluster/interface/T2ROGeometry.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include <iostream>
#include <cmath>

#define NUMBER_OF_STRIPS 256

int min_strip, max_strip;

/**
 *
 */

T2StripHitSim::T2StripHitSim(const edm::ParameterSet & paraSet) 
  : theNewStripHits(),
    diffCoeff_(paraSet.getParameter<std::vector<double> >("diffCoeff")),
    NUMBER_OF_SIM_STEPS_(paraSet.getParameter<int>("NUMBER_OF_SIM_STEPS")),
    eIonE_(paraSet.getParameter<double>("eIonE")),
    gain_(paraSet.getParameter<std::vector<double> >("gain")),
    z_max_(paraSet.getParameter<double>("z_max")),
    z_min_(paraSet.getParameter<double>("z_min")),
    StripWidth_(paraSet.getParameter<std::vector<double> >("StripWidth")){
  
  theT2ROGeometry = new  T2ROGeometry();

  if(gain_.size()!=40)
    {
       std::cout<<"T2Digi Strip Warning: Found "<<gain_.size()<<" detector gains insted of 40. Missing will be put to default values"<<std::endl;
      if(gain_.size()<40)
	{
	  unsigned int miss=40-gain_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	      gain_.push_back(10000.);
	      // std::cout<<"Det "<<gain_.size()<<"th at 10000."<<std::endl;
	    }
	}
    }
 
 
   if(StripWidth_.size()!=40) 
     {
       std::cout<<"T2Digi Strip Warning: Found "<<StripWidth_.size()<<" detector StripWidth_ insted of 40. Missing will be put to default values"<<std::endl;
      if(StripWidth_.size()<40)
	{
	  unsigned int miss=40-StripWidth_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	     StripWidth_.push_back(0.17);
	     //  std::cout<<"Det "<<StripWidth_.size()<<"th at 0.17"<<std::endl;
	    }
	}
    }

   if(diffCoeff_.size()!=40) 
     {
       std::cout<<"T2Digi Strip Warning: Found "<<diffCoeff_.size()<<" detector diffCoeff_ insted of 40. Missing will be put to default values"<<std::endl;
      if(diffCoeff_.size()<40)
	{
	  unsigned int miss=40-diffCoeff_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	      diffCoeff_.push_back(0.43);
	      //  std::cout<<"Det "<<diffCoeff_.size()<<"th at 0.43"<<std::endl;
	    }
	}
    }
}



T2StripHitSim::~T2StripHitSim() {

  delete theT2ROGeometry;
}




unsigned int T2StripHitSim::RawtoSymb(uint32_t thedet)
{
  T2DetId converter;
  unsigned int pl=converter.plane(thedet);
  unsigned int pls=converter.planeSide(thedet);
  unsigned int ht=converter.halfTelescope(thedet);
  unsigned int arm=converter.arm(thedet);
  unsigned int symbolic=pl*2+pls+ht*10+20*arm;
  return symbolic;
}




std::vector<T2DetectorHit> &
T2StripHitSim::simulate(T2Geometry* geom, 
 			const edm::PSimHitContainer & simHits,
			std::vector<double> & chargeDiv,boost::shared_ptr<PSimHitMisaligner> thePSimHitMisaligner) {

  theNewStripHits.clear();
  for (edm::PSimHitContainer::const_iterator hitItr = simHits.begin();
       hitItr != simHits.end();  ++hitItr) {


    // retrieving information of planeSide
    //std::cout<<"iogan1 ps "<<std::endl;
    //std::cout<<geom->getUnitId()<<std::endl;
    T2DetId *dummyT2DetId =new T2DetId(geom->getUnitId());
    unsigned int ps=dummyT2DetId->planeSide(); // ps = 0 or 1 for phasing to IP or not. 
    unsigned int symbdet=RawtoSymb(geom->getUnitId());
    
    //std::cout<<"diogan2 ps "<<ps<<" symbdet "<<symbdet<<std::endl;

    Local3DPoint entry = (*hitItr).entryPoint();
    Local3DPoint exit = (*hitItr).exitPoint();
    Local3DPoint position = (*hitItr).localPosition();
    float energyDeposition = (*hitItr).energyLoss();

     //Simulate the hit displacements
    if(thePSimHitMisaligner->SimulationActivated()==true)
      {
	// std::cout<<"Before"<<(*hitItr).entryPoint().x()<<std::endl;
	  PSimHit ahit = thePSimHitMisaligner->DisplacePSimHit((*hitItr),geom->getUnitId());
	  //std::cout<<"After1"<<ahit.entryPoint().x()<<std::endl;	  
	  entry = ahit.entryPoint();
	  exit = ahit.exitPoint();
	  position = ahit.localPosition();
	  energyDeposition = ahit.energyLoss();
      }
    

    //get total number of generated electrons (E/W*GAIN, W=30 eV/electron, Gain = 8000)
    //  std::cout<<" Gain-det"<<gain_.at(symbdet)<<"-"<<symbdet<<std::endl;
    float totNrElectrons = energyDeposition * 1.0e9 / eIonE_ * gain_.at(symbdet);
    // std::cout<<"totPRIMNrElectrons: "<<totNrElectrons/gain_.at(symbdet)<<std::endl;
    LogDebug("T2StripHitSim") <<"TotNrElec: "<< totNrElectrons;
    // get radii of entry and exit and mean of PSimHit
    double r_in = sqrt( entry.x() * entry.x() + entry.y() * entry.y() );
    double r_out = sqrt( exit.x() * exit.x() + exit.y() * exit.y() ); 
    //    int nSRow = geom->getNearestStripRow(&position);
    int nSCol = geom->getNearestStripCol(&position);

    //get dimensions of strip detector
    double r_min = theT2ROGeometry->GetStripRMin(0,nSCol) + 0.04; // adding haf of the stripwidth 
    // to get the center of the strip
    double r_max = theT2ROGeometry->GetStripRMax(255,nSCol) - 0.04; // here substracting

    
    //    LogDebug("T2StripHitSim") << "LOCAL HIT: rIn " << r_in << ", rOut " << r_out 
    //			      << ", rMin " << r_min << ", rMax " << r_max;

    double charges[NUMBER_OF_STRIPS];
    
    for (int i=0; i<NUMBER_OF_STRIPS; i++) {
      
      charges[i]=0.0;
    }
    
    // Do the simmulation
    
    min_strip = NUMBER_OF_STRIPS;
    max_strip = -1;
    
 

    for (int simstep = 0; simstep<NUMBER_OF_SIM_STEPS_; simstep++) {
      
      double r = r_in + ((r_out-r_in) * (2 * simstep+1))/(2*NUMBER_OF_SIM_STEPS_);
      
      double z = (ps==0) ? z_max_ + ((z_min_-z_max_) * (2 * simstep+1))/(2*NUMBER_OF_SIM_STEPS_) : 
	z = z_min_ + ((z_max_-z_min_) * (2 * simstep+1))/(2*NUMBER_OF_SIM_STEPS_) ;
      
      // std::cout << "Simstep " << simstep << ": r=" << r << "mm, z=" << z << "mm" << std::endl;
      
      calculateStripCharges(charges, r, z, r_min, r_max,geom);
    }
    
    LogDebug("T2StripHitSim") << "min_strip: " << min_strip << ", max_strip: " << max_strip;
    
    double totalCharge = 0.0;

    for (int strip = min_strip; strip<=max_strip; strip++) {
      
      LogDebug("T2StripHitSim") << "strip " << strip << ": q=" << charges[strip];
      float nrStripElectrons = totNrElectrons / 100.0 * charges[strip];  

      T2DetectorHit newStripHit(0, nrStripElectrons, strip, nSCol, 0, &(*hitItr) );
      theNewStripHits.push_back(newStripHit);
      totalCharge += charges[strip];
      
    }

    chargeDiv.push_back(totalCharge/100);

    delete dummyT2DetId;
  }
  
  return theNewStripHits;
  
} // simulate

/**
 *
 */

double T2StripHitSim::calculateSigma(double diffCoeff, double zPosition) {

  return diffCoeff * sqrt(zPosition);

} // calculateSigma

/**
 *
 */
double T2StripHitSim::calculateM3(double sigma) {
  
  return 0.7043 * pow(sigma,-1.004);
  
} // calculateM3




double T2StripHitSim::calculateStripChargesCollPart(double m3, double DMin, double DMax) {  
	double StripCharge = 0.0;
	double erf1 = 0.0;
	double erf2 = 0.0;

	erf1 = erf( m3 * (fabs(DMin)) );
	erf2 = erf( m3 * (fabs(DMax)) );

	if ((DMax/DMin) >=0){
	StripCharge = fabs( 50*(1 - erf1) - 50*( 1 - erf2) );
	} else StripCharge = fabs( 100 - 50*(1 - erf1) - 50*( 1 - erf2) );
  
  return StripCharge;

} // calculateStripChargesPart


/**
 *  charge in percent of total number of electrons
 */


/** NEW VERSION **/
void T2StripHitSim::calculateStripCharges(double* charges, double r, double z, 
					  double r_min, double r_max,T2Geometry* geom) {

  //std::cout<<"Entro T2StripHitSim::calculateStripCharges"<<std::endl;


// We need to define the effective strip width. It represents the effective width (bigger than the real) of the strip and 
// it takes into account the focusing effects of the electric filed lines on strips.
// StripWidth will represent therefore the effective width (bigger than the real) of the strip. 
// We will put here 0.210mm. It could be better probably to insert this number in tha configuration file. For the moment we will put it here.

  //double stripWidth = 0.210;
  unsigned int symbdet=RawtoSymb(geom->getUnitId());
  double stripWidth = StripWidth_.at(symbdet) ;
  double sigma = calculateSigma(diffCoeff_.at(symbdet), z);
  //  std::cout<<"stripWidth diffCoeff_ "<<stripWidth<<diffCoeff_.at(symbdet)<<" Detector "<<symbdet<<std::endl;
  double m3 = calculateM3(sigma);

// We need to define the number of strips for which we will calculate the charge collected.
// We will obtain this number according to the sigma used (i.e. n times the sigma).
// We will put here 6. It will be better to insert this number in tha configuration file. For the moment we will put it here.

  int n_sigma = 6; 
  double delta_r = n_sigma * sigma;
 
  double r0 = r-delta_r;
  double r1 = r+delta_r;

  if (r0<r_min) { 
   
    r0 = r_min;
  }
  else if (r0>r_max) {
    
    r0 = r_max;
  }
  
  if (r1<r_min) { 
    
    r1 = r_min;
  }
  else if (r1>r_max) {

    r1 = r_max;
  }
  
  int s0 = static_cast<int>(ceil((r0-r_min)/(r_max-r_min)*(NUMBER_OF_STRIPS-1)));
  int s1 = static_cast<int>(floor((r1-r_min)/(r_max-r_min)*(NUMBER_OF_STRIPS-1)));

  if (s0 >= NUMBER_OF_STRIPS) {

    s0=NUMBER_OF_STRIPS-1;
  }

  if (s1 >= NUMBER_OF_STRIPS) {

    s1=NUMBER_OF_STRIPS-1;
  }
  
  //  std::cout << "   delta_r= " << delta_r << "mm" << std::endl;
  //  std::cout << "   s0     = " << s0 << std::endl;
  //  std::cout << "   s1     = " << s1 << std::endl;

  

  for(int strip = s0; strip<=s1; strip++) {

 //   double dist = fabs(r - (r_min + (r_max-r_min)/(NUMBER_OF_STRIPS-1)*strip));
	double dRMin = (r - (r_min + (r_max-r_min)/(NUMBER_OF_STRIPS-1)*strip)) - (stripWidth/2);
	double dRMax = (r - (r_min + (r_max-r_min)/(NUMBER_OF_STRIPS-1)*strip)) + (stripWidth/2);

	double q = calculateStripChargesCollPart( m3, dRMin, dRMax);


    charges[strip] += q / NUMBER_OF_SIM_STEPS_;

    // std::cout << "   strip " << strip << ": +"  << q / NUMBER_OF_SIM_STEPS_ << std::endl;

    if (strip<min_strip) min_strip = strip;
    if (strip>max_strip) max_strip = strip;
    
     
  }
 //std::cout<<" T2StripHitSim::calculateStripCharges laststrip esco "<<std::endl;
} // calculateStripCharges




