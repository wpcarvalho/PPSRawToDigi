/** 
 * Class simulates hit on a pad in the T2 detector
 *
 * Author: Erik Brücken / University of Helsinki
 * email:  brucken@cc.helsinki.fi
 * Date:   2007-11-26
 */

#include "SimTotem/T2Digitizer/interface/T2PadHitSim.h"
#include "SimTotem/T2Digitizer/interface/T2DetectorHit.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/TotemGeometry/interface/T2Geometry.h"
#include "DataFormats/T2Cluster/interface/T2ROGeometry.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include <iostream>
#include <cmath>

int min_pad_row, max_pad_row, min_pad_col, max_pad_col;

/**
 *
 */

T2PadHitSim::T2PadHitSim(const edm::ParameterSet & paraSet) 
  : theNewPadHits(),     
    diffCoeff_(paraSet.getParameter<std::vector<double> >("diffCoeff")),
    NUMBER_OF_SIM_STEPS_(paraSet.getParameter<int>("NUMBER_OF_SIM_STEPS")),
    eIonE_(paraSet.getParameter<double>("eIonE")),
    gain_(paraSet.getParameter<std::vector<double> >("gain")),
    z_max_(paraSet.getParameter<double>("z_max")),
    z_min_(paraSet.getParameter<double>("z_min")) {
  
  theT2ROGeometry = new  T2ROGeometry();


    if(gain_.size()!=40)
    {
      std::cout<<"T2Digi Pad Warning: Found "<<gain_.size()<<" detector gains insted of 40. Missing will be put to default values"<<std::endl;
      if(gain_.size()<40)
	{
	  unsigned int miss=40-gain_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	      gain_.push_back(10000.);
	      //std::cout<<"Det "<<gain_.size()<<"th at 10000."<<std::endl;
	    }
	}
    }
  
   

   if(diffCoeff_.size()!=40) 
     {
       std::cout<<"T2Digi Pad  Warning: Found "<<diffCoeff_.size()<<" detector diffCoeff_ insted of 40. Missing will be put to default values"<<std::endl;
      if(diffCoeff_.size()<40)
	{
	  unsigned int miss=40-diffCoeff_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	      diffCoeff_.push_back(0.43);
	      // std::cout<<"Det "<<diffCoeff_.size()<<"th at 0.43"<<std::endl;
	    }
	}
    }




} // T2PadHitSim

/**
 *
 */

T2PadHitSim::~T2PadHitSim() {

  delete theT2ROGeometry;

} // ~T2PadHitSim

/**
 *
 */

unsigned int T2PadHitSim::RawtoSymb(uint32_t thedet)
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
T2PadHitSim::simulate(T2Geometry* geom, 
		      const edm::PSimHitContainer & simHits,
		      std::vector<double> & chargeDiv,boost::shared_ptr<PSimHitMisaligner> thePSimHitMisaligner) {

  theNewPadHits.clear();


  for (edm::PSimHitContainer::const_iterator hitItr = simHits.begin();
       hitItr != simHits.end();  ++hitItr) {

    // getting the planeSide info 0/1 for phasing/not phasing the IP
    T2DetId *dummyT2DetId = new T2DetId(geom->getUnitId());
    unsigned int ps=dummyT2DetId->planeSide();
    // constructing Readout Geometry of T2-plane
    //    T2ROGeometry* geoT2 = new T2ROGeometry(geom->getUnitId());
    theT2ROGeometry->SetT2DetType(geom->getUnitId());
    
    unsigned int symbdet=RawtoSymb(geom->getUnitId());
    //std::cout<<"Pad processing for plane "<<symbdet<<std::endl;

    Local3DPoint entry = (*hitItr).entryPoint();
    Local3DPoint exit = (*hitItr).exitPoint();
    Local3DPoint position = (*hitItr).localPosition();
    float energyDeposition = (*hitItr).energyLoss();

     //Simulate the hit displacements
    if(thePSimHitMisaligner->SimulationActivated()==true)
      {
	//std::cout<<"Before"<<(*hitItr).entryPoint().x()<<std::endl;
	  PSimHit ahit = thePSimHitMisaligner->DisplacePSimHit((*hitItr),geom->getUnitId());
	  //std::cout<<"After1"<<ahit.entryPoint().x()<<std::endl;	  
	  entry = ahit.entryPoint();
	  exit = ahit.exitPoint();
	  position = ahit.localPosition();
	  energyDeposition = ahit.energyLoss();
      }



    
      

    //get total number of generated electrons (E/W*GAIN, W ~ 30 eV/electron, Gain ~ 8000)
    float totNrElectrons = energyDeposition * 1.0e9 / eIonE_ * gain_.at(symbdet);
    
    
    //Just for debugging purpose: this activate the pad only belonging to a given Geant track id. 
    //Warning: the strip digitization is unaffected by this command.
    /*
    if((*hitItr).trackId()!=67)
      totNrElectrons =0.;
    */
    
    LogDebug("T2PadHitSim") <<"TotNrElec: "<< totNrElectrons;
    // get radii of entry and exit and mean of PSimHit
    double r_in = sqrt( entry.x() * entry.x() + entry.y() * entry.y() );
    double r_out = sqrt( exit.x() * exit.x() + exit.y() * exit.y() ); 
    double x_in = entry.x(); 
    double x_out = exit.x(); 
    double y_in = entry.y(); 
    double y_out = exit.y(); 
    double phi_in = (x_in==0 && y_in==0) ? 0 : atan2( y_in, x_in );
    double phi_out = (x_out==0 && y_out==0) ? 0 : atan2( y_out, x_out );

    double padCharge[NUMBER_OF_ROWS][NUMBER_OF_COLS];
    
    for (int i=0; i<NUMBER_OF_ROWS; i++) {
      
      for (int j=0; j<NUMBER_OF_COLS; j++) {
	
	padCharge[i][j]=0.0;
      }
    }
    
    // Do the simmulation for pads
    
    min_pad_row = NUMBER_OF_ROWS;
    max_pad_row = -1;
    min_pad_col = NUMBER_OF_COLS;
    max_pad_col = -1;
    

   

    for (int simstep = 0; simstep<NUMBER_OF_SIM_STEPS_; simstep++) {

      double r = r_in + ((r_out-r_in) * (2 * simstep+1))/(2*NUMBER_OF_SIM_STEPS_);
      double phi =  phi_in + ((phi_out-phi_in) * (2 * simstep+1))/(2*NUMBER_OF_SIM_STEPS_);

      int rowNearestPad = geom->getNearestPadRow(r, phi);
      int colNearestPad = geom->getNearestPadCol(r, phi);
      
      double z = (ps==0) ? z_max_ + ((z_min_-z_max_) * (2 * simstep+1))/(2*NUMBER_OF_SIM_STEPS_) : 
	z = z_min_ + ((z_max_-z_min_) * (2 * simstep+1))/(2*NUMBER_OF_SIM_STEPS_) ;

      calculate_pad_charge(theT2ROGeometry, padCharge, rowNearestPad, colNearestPad, r, phi, z);
    }
    
    LogDebug("T2PadHitSim") << "min_pad_row: " << min_pad_row << ", max_pad_row: " << max_pad_row;
    LogDebug("T2PadHitSim") << "min_pad_col: " << min_pad_col << ", max_pad_col: " << max_pad_col;
    
    for (int padR = min_pad_row; padR<=max_pad_row; padR++) {
   
      for (int padC = min_pad_col; padC<=max_pad_col; padC++) {
	
	LogDebug("T2PadHitSim") << "pad " << padR << " "<< padC << ": q=" << padCharge[padR][padC];
	
	// charge on pad obtained by substracting electrons that had been collected by the strips
	float nrPadElectrons = totNrElectrons * (1.0 - chargeDiv[0]) * padCharge[padR][padC];  

	T2DetectorHit newPadHit(0, nrPadElectrons, padR, padC, 0, &(*hitItr) );
	theNewPadHits.push_back(newPadHit);
      }
    }

    delete dummyT2DetId;
    //  delete geoT2;
  }
  
  return theNewPadHits;
  
} // simulate

/**
 *
 */

double T2PadHitSim::calculateSigma(double diffCoeff, double zPosition) {

  return diffCoeff * sqrt(zPosition);

} // calculateSigma

/**
 *
 */

double T2PadHitSim::calculateM3(double sigma) {
  
  return 0.7043 * pow(sigma,-1.004);
  
} // calculateM3

/**
 *
 */


double T2PadHitSim::calculatePadChargeCollPart(double m3, double DMin , double DMax) {  
  double PadCharge = 0.0;
  double erf1 = 0.0;
  double erf2 = 0.0;

  erf1 = erf( m3 * (fabs(DMin)) );
  erf2 = erf( m3 * (fabs(DMax)) );

	if ((DMax/DMin) >=0){
	PadCharge = fabs( 50*(1 - erf1) - 50*( 1 - erf2) );
	} else PadCharge = fabs( 100 - 50*(1 - erf1) - 50*( 1 - erf2) );

  return PadCharge;
  

} // calculatePadChargeCollPart



/**
 *  charge in percent of total number of electrons
 */

void T2PadHitSim::calculate_pad_charge(T2ROGeometry* geometry, double (*charges)[NUMBER_OF_COLS], 
				       int row, int col, double r, double phi,  double z) {
   

  unsigned int symbdet=RawtoSymb(geometry->GetDetId());
  double sigma = calculateSigma(diffCoeff_.at(symbdet), z);
  double m3 = calculateM3(sigma);

  int pR0 = row-3;
  int pR1 = row+3;
  int pC0 = col-3;
  int pC1 = col+3;
  
  if (pR0 >= NUMBER_OF_ROWS) pR0 = NUMBER_OF_ROWS-1;
  if (pR1 >= NUMBER_OF_ROWS) pR1 = NUMBER_OF_ROWS-1;
  if (pC0 >= NUMBER_OF_COLS) pC0 = NUMBER_OF_COLS-1;
  if (pC1 >= NUMBER_OF_COLS) pC1 = NUMBER_OF_COLS-1;
  if (pR0 < 0) pR0 = 0;
  if (pR1 < 0) pR1 = 0;
  if (pC0 < 0) pC0 = 0;
  if (pC1 < 0) pC1 = 0;
  
  for (int padRow = pR0; padRow<=pR1; padRow++) {
    
    double pRadMin = geometry->GetPadRMin(padRow, 1);
    double pRadMax = geometry->GetPadRMax(padRow, 1);

// We need to define the distances of the borders of tha pads from the center 
// of the electrons cloud
	double dRMin = 0.0;
	double dRMax = 0.0;

    //    std::cout << "RowStuff: " << padRow << " " << pRadMin << " " << pRadMax << std::endl;
    
/**    if (r > ((pRadMin + pRadMax) / 2)) {
      
      dR = r - pRadMax;
      
    } else dR = pRadMin -r;
*/

	dRMin = (r - pRadMin);
	dRMax = (r - pRadMax);



    double qRow = calculatePadChargeCollPart( m3, dRMin, dRMax);

    //    std::cout << "ChargeRow " << qRow << std::endl;

    for (int padCol = pC0; padCol<=pC1; padCol++) {
      
      double pPhiMin = geometry->GetPadPhiMinLocal(padRow, padCol);
      double pPhiMax = geometry->GetPadPhiMaxLocal(padRow, padCol);

	  // We need to define the distances of the borders of tha pads from the center 
// of the electrons cloud
		double dPhiMin = 0.0;
		double dPhiMax = 0.0;

      //      std::cout << "ColStuff: " << padCol << " " << pPhiMin << " " << pPhiMax << std::endl;

/**      if (phi > ((pPhiMin + pPhiMax) / 2)) {
      
        dPhi = phi - pPhiMax;
      
      } else dPhi = pPhiMin -phi;
*/
	
		dPhiMin = (phi - pPhiMin);
		dPhiMax = (phi - pPhiMax);


      double qCol = calculatePadChargeCollPart( m3, dPhiMin*r, dPhiMax*r);

      //      std::cout << "    ChargeCol " << qCol << std::endl;//" dPhi " << dPhi << std::endl;
  
      charges[padRow][padCol] += qRow * qCol / NUMBER_OF_SIM_STEPS_ / 10000;
      
      if (padRow<min_pad_row) min_pad_row = padRow;
      if (padRow>max_pad_row) max_pad_row = padRow;
      if (padCol<min_pad_col) min_pad_col = padCol;
      if (padCol>max_pad_col) max_pad_col = padCol;
    }
  }
} // calculate_pad_charge

