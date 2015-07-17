#include "SimTotem/T1Digitizer/interface/T1WireHitSim.h"
#include "SimTotem/T1Digitizer/interface/T1DriftSim.h"
#include "SimTotem/T1Digitizer/interface/T1CrossGap.h"
#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"
#include "SimTotem/T1Digitizer/interface/T1GasCollisions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/TotemGeometry/interface/T1ChamberSpecs.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "CLHEP/Units/SystemOfUnits.h"

//#define _DEBUG_
//#define _LOGINFO_

T1WireHitSim::T1WireHitSim(T1DriftSim* driftSim) 
  : pDriftSim(driftSim),
    theGasIonizer( new T1GasCollisions() ) 
{
//  std::cout << " Creating T1WireHitSim  ----------------- "<<std::endl;
}


T1WireHitSim::~T1WireHitSim() {
  delete theGasIonizer;
}


std::vector<T1DetectorHit> &
T1WireHitSim::simulate(uint32_t mydetid, const edm::PSimHitContainer & simHits, T1Geometry *layer) 
{

  //  std::cout << " FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"<<std::endl;
  //  std::cout << mydetid <<std::endl;

  T1DetId mydetId(mydetid);
  // std::cout << mydetid << " " << mydetId << std::endl;


  newWireHits.clear();
  int iii = 0;
  for (edm::PSimHitContainer::const_iterator hitItr = simHits.begin();
       hitItr != simHits.end();  ++hitItr)
    {

      iii++;
#ifdef _DEBUG_
      std::cout << " ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  SimHit numero "<< iii <<" "<< (*hitItr).entryPoint().x() << " "<< (*hitItr).entryPoint().y() << " "<< (*hitItr).entryPoint().z() << " "<< (*hitItr).energyLoss() << std::endl;
      edm::LogInfo("T1wireHitSim") <<  " ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  SimHit numero "<< iii <<" "<< (*hitItr).entryPoint().x() << " "<< (*hitItr).entryPoint().y() << " "<< (*hitItr).entryPoint().z() << " "<<(*hitItr).energyLoss() <<"  Ptype"<<(*hitItr).particleType();
#endif
 

      std::vector<LocalPoint> ionClusters = getIonizationClusters(*hitItr, layer);
#ifdef _DEBUG_
      bool fff=true;
#endif

#ifdef _DEBUG_
      std::cout << "ionCluster size"<< ionClusters.size() << std::endl;
#endif
      for(int icl = 0; icl < int( ionClusters.size() ); ++icl) {

	// Drift the electrons in the cluster to the nearest wire...
    

	int nearestWire =layer->tellWireTB(mydetid,ionClusters[icl].x(),ionClusters[icl].y());

	//   int nearestWire =layer->tellWireTB(mydetid,0,200);

#ifdef _DEBUG_
	if(fff){
	  //%// std::cout << "Nearest wire y : " << ionClusters[icl].y() <<std::endl;
	  //	if(nearestWire > 50 && mydetId.Plane()==0)
	  std::cout << "Nearest wire : " << nearestWire <<" in plane " << mydetId.Plane()<<std::endl;
	  fff=false;
	}
#endif
	// The wire hit contains wire# and position measured _along the wire_
	// from where it intersects local y axis.
	//    std::cout << "T1WireHitSim: simulate:  nearestWire " << nearestWire <<"  **********************" << std::endl;

	T1DetectorHit temp_hit = pDriftSim->getWireHit(ionClusters[icl], layer, nearestWire,*hitItr);


	newWireHits.push_back(temp_hit);
	//     newWireHits.push_back(pDriftSim->getWireHit(ionClusters[icl], layer, nearestWire,*hitItr) );

      }
    } 
  //%// std::cout << "T1WireHitSim: simulate:  end  **********************" << std::endl;
  return newWireHits;
}

std::vector<LocalPoint> 
T1WireHitSim::getIonizationClusters(const PSimHit & simHit, T1Geometry * layer) {

  // std::cout << "getIonizationClusters: simulate:  start  **********************" << std::endl;
  //  const LocalPoint & entryPoint = simHit.entryPoint();
  //  const LocalPoint & exitPoint  = simHit.exitPoint();

#ifdef _LOGINFO_
  edm::LogInfo("T1WireHitSim") << "T1WireHitSim:" 
			       << " type=" << simHit.particleType() 
			       << " mom=" << simHit.pabs()
			       << "\n Local entry " << simHit.entryPoint << " exit " << exitPoint;
#endif

  std::vector<LocalPoint> positions;
  std::vector<int> electrons;
  theGasIonizer->simulate( simHit, positions, electrons );

  //%// std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@ " << positions.size() << std::endl;

  //  std::vector<LocalPoint> results( positions ); //copy
  std::vector<LocalPoint> results; // start empty

  int j = 0;
  for( std::vector<LocalPoint>::const_iterator i = positions.begin(); 
       i != positions.end(); ++i ) {
    ++j;
    LocalPoint newPoint( (*i).x()*10.,(*i).y()*10.,(*i).z()*10. );
    // some verification
    // if(layer->geometry()->inside(newPoint) ) {



    if(layer->inside(simHit.detUnitId(),newPoint) ) {
      //%// std::cout << "getIonizationCluster: if inside is true **********************" << std::endl;
      // push the point for each electron at this point
      
      for( int ie = 1;  ie <= electrons[j-1]; ++ie ) {
        results.push_back(newPoint);
      }
    }
  }
  //%// std::cout << "getIonizationCluster: number of cluster "<<  results.size() <<" **********************" << std::endl;
#ifdef _LOGINFO_
  edm::LogInfo("T1WireHitSim") << "T1WHS: there are " << results.size()
			       << " clusters identified with each electron.";
#endif
  return results;
}


