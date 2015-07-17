#include "RecoTotemT1T2/T1RecHit/interface/T1RecHitProducer.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"
#include "DataFormats/T1Cluster/interface/T1Cluster.h"
#include "DataFormats/T1Cluster/interface/T1ClusterCollection.h"

#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
//#include "RecoTotemT1T2/T1RecHit/interface/T1Clusterizer.h"
//#include "RecoTotemT1T2/T1RecHit/interface/T1ClusterContainer.h"
//#include "RecoTotemT1T2/T1RecHit/interface/T1Cluster.h"

#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include <string>


using namespace edm;
using namespace std;

//#define _DEBUG_


T1RecHitProducer::T1RecHitProducer(const ParameterSet& config){
  // Set verbose output

  t1ClusterCollectionLabel = config.getParameter<InputTag>("T1ClusterCollectionLabel");
  t1DigiWireCollectionLabel = config.getParameter<InputTag>("T1DigiWireCollectionLabel");
  produces<T1RecHit2DCollection>();
  std::string thegeometryfile = "Geometry/TotemGeometry/data/T1_data_geometry.dat";
  theT1Geometry = new T1Geometry(thegeometryfile);
  theChiSquare = config.getParameter<double>("ChiSquare");

  theVerbosity = config.getParameter<int>("Verbosity");

  if(theVerbosity >=1 )cout << " Inside T1RecHitProducer "<<endl;

}

T1RecHitProducer::~T1RecHitProducer(){
  delete theT1Geometry;
  //delete RecHitMatrix;
}



void T1RecHitProducer::produce(Event& event, const EventSetup& setup) {

  // Get the digis from the event
  Handle<T1DigiWireCollection> Wdigis;
  Handle<T1ClusterCollection> Sdigis;

  event.getByLabel(t1DigiWireCollectionLabel, Wdigis);
  event.getByLabel(t1ClusterCollectionLabel, Sdigis);
  
  std::vector< edm::OwnVector<T1RecHit2D> > RecHitMatrix;

  // Create the pointer to the collection which will store the rechits
  auto_ptr<T1RecHit2DCollection> recHitCollection(new T1RecHit2DCollection());

  // load the geometry
  //  T1Geometry * theT1Geometry = new T1Geometry;

  T1ChamberSpecs parametri;

  // Iterate through all digi collections ordered by LayerId   
  T1DigiWireCollection::DigiRangeIterator t1wdgIt;

  T1ClusterCollection::DigiRangeIterator t1sdgIt;
  T1ClusterCollection::DigiRangeIterator t1sdgIt2;

  for(int g=0; g<60; g++){
    RecHitMatrix.push_back( OwnVector<T1RecHit2D>());
  }
  
  
  for (t1wdgIt = Wdigis->begin(); t1wdgIt != Wdigis->end();
       ++t1wdgIt){
    
    // The layerId
    const T1DetId& t1Id = (*t1wdgIt).first;
	
    const T1DigiWireCollection::Range& range = (*t1wdgIt).second;
    
    for (T1DigiWireCollection::const_iterator digiIt = range.first;
	 digiIt!=range.second;++digiIt){
    // std::cout << t1Id << std::endl;
      for (t1sdgIt = Sdigis->begin(); t1sdgIt != Sdigis->end();
	   ++t1sdgIt){
      //   std::cout << t1Id << " ------ " << (*t1sdgIt).first << std::endl;
	//	std::cout << t1Id.rawId()<< std::endl;
	if( 
	   (*t1sdgIt).first.Arm()== t1Id.Arm() &&
	   (*t1sdgIt).first.Plane()== t1Id.Plane() &&
	   (*t1sdgIt).first.CSC()== t1Id.CSC() &&
	   (*t1sdgIt).first.Layer() == 1
	   ){
	   
	  const T1ClusterCollection::Range& srange = (*t1sdgIt).second;
	  
	  for (T1ClusterCollection::const_iterator sdigiIt = srange.first;
	       sdigiIt!=srange.second;++sdigiIt){

	    for (t1sdgIt2 = Sdigis->begin(); t1sdgIt2 != Sdigis->end();
		 ++t1sdgIt2){
	      //   std::cout << t1Id << " ------ " << (*t1sdgIt).first << std::endl;
	      //	std::cout << t1Id.rawId()<< std::endl;
	      if( 
		 (*t1sdgIt2).first.Arm()== t1Id.Arm() &&
		 (*t1sdgIt2).first.Plane()== t1Id.Plane() &&
		 (*t1sdgIt2).first.CSC()== t1Id.CSC() &&
		 (*t1sdgIt2).first.Layer() == 2
		 ){
		  
		const T1ClusterCollection::Range& srange2 = (*t1sdgIt2).second;
		  
		for (T1ClusterCollection::const_iterator sdigiIt2 = srange2.first;
		     sdigiIt2!=srange2.second;++sdigiIt2){
		  
		  float yW = theT1Geometry->yOfWire(t1Id.rawId(), (*digiIt).wire() );
//		  float strisciaA=0.0;
		  /*
		    if((*sdigiIt).balance()==0) strisciaA=(float)(*sdigiIt).Center()-0.25;
		    if((*sdigiIt).balance()==1) strisciaA=(float)(*sdigiIt).Center()+0.25;
		    float strisciaB=0.0;
		    if((*sdigiIt2).balance()==0) strisciaB=(float)(*sdigiIt2).Center()-0.25;
		    if((*sdigiIt2).balance()==1) strisciaB=(float)(*sdigiIt2).Center()+0.25;

		    float yS = theT1Geometry->yOfStripCrossing(t1Id.rawId(),strisciaA,strisciaB );
		    float xS = theT1Geometry->xOfStripCrossing(t1Id.rawId(),strisciaA,strisciaB );
		    float xSAW = theT1Geometry->xOfStripAWireCrossing(t1Id.rawId(),strisciaA, (*digiIt).wire());
		    float xSBW = theT1Geometry->xOfStripBWireCrossing(t1Id.rawId(),strisciaB, (*digiIt).wire());
		  */
		  //   std::cout << "###################### FLAG 2 ########################"<<std::endl;
		      
		  float yS = theT1Geometry->yOfStripCrossing(t1Id.rawId(),(*sdigiIt).Center(),(*sdigiIt2).Center() );
		  float xS = theT1Geometry->xOfStripCrossing(t1Id.rawId(),(*sdigiIt).Center(),(*sdigiIt2).Center() );
		  ///
		  //  float xSAW = theT1Geometry->xOfStripAWireCrossing(t1Id.rawId(),(*it).Center(), (*digiIt).wire());
		  //	  float xSBW = theT1Geometry->xOfStripBWireCrossing(t1Id.rawId(),(*it2).Center(), (*digiIt).wire());
		      
		  ////
		  //     std::cout << "###################### FLAG 3 ########################"<<std::endl;
		  //#//		  std::cout << "y del filo " << yW <<"  y delle strip cross " <<yS << "  x delle strip cross " <<xS << "  x strip A e wire "<< xSAW << "  x strip B e wire "<< xSBW <<  std::endl;
		      
		  float sigma_yW = (1.5);// /RADQ3;
		  //  float sigma_yS = 2./4. * parametri.Pitch();// / RADQ3;
		  //	  float sigma_xS = (2./4. * parametri.Pitch() / RADQ3); // /RADQ3;
		  float sigma_yS = (*sdigiIt).Sigma();// / RADQ3;
		  float sigma_xS = (*sdigiIt).Sigma() / RADQ3; // /RADQ3;
		  // float sigma_xSW = 1.5 + 1./4. * parametri.Pitch() / RADQ3;
		      
		  // float mean_x = (xS/sigma_xS/sigma_xS + xSAW/sigma_xSW/sigma_xSW + xSBW/sigma_xSW/sigma_xSW )/ (2./sigma_xSW/sigma_xSW + 1./sigma_xS/sigma_xS);
		  float mean_x=xS;
		  float weighted_mean_y = (yW/sigma_yW/sigma_yW + yS/sigma_yS/sigma_yS)/(1./sigma_yW/sigma_yW + 1./sigma_yS/sigma_yS);
		  float mean_y = yS;
		      
		  // float sigma_x = sqrt( 1./( 1./sigma_xS/sigma_xS + 2./sigma_xSW/sigma_xSW ) );
		  float sigma_x=sigma_xS;
		  float sigma_y =  sqrt( 1./(1./sigma_yW/sigma_yW + 1./sigma_yS/sigma_yS) );
		  //   std::cout << "###################### FLAG 4 ########################"<<std::endl;
		  //			    float chi_squared =sqrt( (mean_x - xS)*(mean_x - xS)/sigma_xS/sigma_xS + (mean_x - xSAW)*(mean_x - xSAW)/sigma_xSW/sigma_xSW + (mean_x-xSBW)*(mean_x-xSBW)/sigma_xSW/sigma_xSW + (mean_y - yW)*(mean_y - yW)/sigma_yW/sigma_yW + (mean_y - yS)*(mean_y - yS)/sigma_yS/sigma_yS );
		      
		  float lamda = fabs(yW-yS)/sqrt((sigma_yW)*(sigma_yW)+(sigma_yS)*(sigma_yS));
		      
		  //	    std::cout << "Mean X " << mean_x << "   Mean Y " << mean_y << "   chi "<<chi_squared<< std::endl;
if(theVerbosity >=2){		 
		  std::cout << "XS " <<xS << "  YS " << yS << "  YW "<<yW << std::endl;
		  std::cout << "Mean X " << mean_x << "   Mean Y " << mean_y <<" Lamda  " << lamda << std::endl;
}
		 
		  // inserire un controllo sulla geometria !!!
		  // gli hits ricostruiti stanno dentro alla camera ??? 
		      
		      
		  //con lamda verifichiamo la compatibilita delle misure
		  if( lamda <= 3 ) {

		    LocalPoint lp(mean_x,mean_y,0.0);
		    LocalPoint lpWeight(mean_x,weighted_mean_y,0.0);
		    LocalPoint lpWire(mean_x,yW,0.0);
		    if(theVerbosity >= 2)
		    std::cout << "T1RecHit sigma_x " << sigma_x << std::endl;
		    LocalError le(sigma_x*sigma_x,0.0,sigma_y*sigma_y); // da verificare la covarianza
			
		    T1RecHit2D* newHit_lp = new T1RecHit2D(t1Id, lp, le, lamda);
		    T1RecHit2D* newHit_lpWeight = new T1RecHit2D(t1Id, lpWeight, le, lamda);
		    T1RecHit2D* newHit_lpWire = new T1RecHit2D(t1Id, lpWire, le, lamda);

		    if((theT1Geometry->insideRECO(t1Id,lp))) {
		      //-------------------------------------------------------------	
//		      T1RecHit2D* newHit = new T1RecHit2D(t1Id, lp, le, lamda);
		      index= t1Id.Arm()*30+ t1Id.Plane()*6+t1Id.CSC();
if(theVerbosity >=2){
		      std::cout << "New rechit in "<<t1Id<<" "<<lp <<std::endl;
}
		      edm::OwnVector<T1RecHit2D>::iterator ittico;
		      bool bandiera = true;
		      for(ittico = RecHitMatrix[index].begin(); ittico !=  RecHitMatrix[index].end() && bandiera; ittico++){
			if((*ittico).localPosition().x() == newHit_lp->localPosition().x() && (*ittico).localPosition().y() == newHit_lp->localPosition().y())
			  bandiera = false;
			    
		      }
		      if(bandiera){
//			if(newHit_lp->localPosition().y()<-220)std::cout << " HIT BASSO " << t1Id << " I " << lp <<std::endl;
			RecHitMatrix[index].push_back(newHit_lp);
			delete newHit_lpWeight;
			delete newHit_lpWire;
		      }else{
		      delete newHit_lp;
		      delete newHit_lpWeight;
		      delete newHit_lpWire;
		      }
		      //-----------------------------------------			  
		    }else if((theT1Geometry->insideRECO(t1Id,lpWeight))) {
		      //-------------------------------------------------------------	
//		      T1RecHit2D* newHit = new T1RecHit2D(t1Id, lpWeight, le, lamda);
		      index= t1Id.Arm()*30+ t1Id.Plane()*6+t1Id.CSC();
if(theVerbosity >=2){
		      std::cout << "New rechit in "<<t1Id<<" "<<lpWeight <<std::endl;
}
		      edm::OwnVector<T1RecHit2D>::iterator ittico;
		      bool bandiera = true;
		      for(ittico = RecHitMatrix[index].begin(); ittico !=  RecHitMatrix[index].end() && bandiera; ittico++){
			if((*ittico).localPosition().x() == newHit_lpWeight->localPosition().x() && (*ittico).localPosition().y() == newHit_lpWeight->localPosition().y())
			  bandiera = false;
			    
		      }
		      if(bandiera){
//	if(newHit_lpWeight->localPosition().y()<-220)std::cout << " HIT BASSO " << t1Id << " II " << lpWeight <<std::endl;
			RecHitMatrix[index].push_back(newHit_lpWeight);
			delete newHit_lp;
			delete newHit_lpWire;
		      }else{
		      delete newHit_lp;
		      delete newHit_lpWeight;
		      delete newHit_lpWire;
		      }
			
		      //-----------------------------------------		
		    }else if((theT1Geometry->insideRECO(t1Id,lpWire))) {
		    //	      T1RecHit2D* newHit = new T1RecHit2D(t1Id, lpWire, le, lamda);
		      index= t1Id.Arm()*30+ t1Id.Plane()*6+t1Id.CSC();
if(theVerbosity >=2){
		      std::cout << "New rechit in "<<t1Id<<" "<<lpWire <<std::endl;
}
		      edm::OwnVector<T1RecHit2D>::iterator ittico;
		      bool bandiera = true;
		      for(ittico = RecHitMatrix[index].begin(); ittico !=  RecHitMatrix[index].end() && bandiera; ittico++){
			if((*ittico).localPosition().x() == newHit_lpWire->localPosition().x() && (*ittico).localPosition().y() == newHit_lpWire->localPosition().y())
			  bandiera = false;
			    
		      }
		      if(bandiera)
			{
//	if(newHit_lpWire->localPosition().y()<-220)std::cout << " HIT BASSO " << t1Id << " III  " << lpWire <<std::endl;
			RecHitMatrix[index].push_back(newHit_lpWire);
			delete newHit_lpWeight;
			delete newHit_lp;
			}			else{
		      delete newHit_lp;
		      delete newHit_lpWeight;
		      delete newHit_lpWire;
			}

		      //-----------------------------------------		
		    }else{
		      if(theVerbosity >= 1){
		      std::cout << " ********************** RECONSTRUCTED HIT OUTSIDE CHAMBER ************************* "<< std::endl;
		      std::cout << "XS " <<xS << "  YS " << yS << "  YW "<<yW << std::endl; }		
		      delete newHit_lp;
		      delete newHit_lpWeight;
		      delete newHit_lpWire;
		      //		       std::cout << (*digiIt).wire() <<std::endl;
		    }//else
			
		  }// if( lamda <= 2 ) 
		}
		  
	      }
	    }
	      
	  }	
	}
      }
    
    }
  }// loop DigiWire
  for(int jj=0; jj<60; jj++){
    int braccio=jj/30;
    int piano=(jj%30)/6;
    int camera=((jj%30)%6);
    //	std::cout << "###################### FLAG 7 ########################"<<std::endl;

    T1DetId temp(braccio,piano,camera);
    if(  RecHitMatrix[jj].size() > 0) 
      recHitCollection->put(temp,RecHitMatrix[jj].begin(),RecHitMatrix[jj].end());
    //			std::cout << "###################### FLAG 8 ########################"<<std::endl;
  }
  event.put(recHitCollection);
  //	std::cout << "###################### FLAG 9 ########################"<<std::endl;

  //  delete theT1Geometry;

}

