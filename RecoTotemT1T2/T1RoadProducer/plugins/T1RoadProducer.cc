/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
*/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"
#include "DataFormats/T1Cluster/interface/T1Cluster.h"
#include "DataFormats/T1Cluster/interface/T1ClusterCollection.h"

#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"

#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"

#include "RecoTotemT1T2/T1RoadProducer/interface/T1RoadProducer.h"

#include <string>


using namespace edm;
using namespace std;

//#define _DEBUG_


T1RoadProducer::T1RoadProducer(const ParameterSet& config):theVerbosity(0),theAlignment(false){
  // Set verbose output
  t1RecHit2DCollectionLabel = config.getParameter<edm::InputTag>("T1RecHit2DCollectionLabel");

  produces<T1RoadCollection>();
  std::string thegeometryfile = "Geometry/TotemGeometry/data/T1_data_geometry.dat";
  theT1Geometry = new T1Geometry(thegeometryfile);
  theDeltaEta = config.getParameter<double>("DeltaEta");
  theDeltaPhi = config.getParameter<double>("DeltaPhi");
  theMinHits = config.getParameter<int>("MinHits");
  theMaxHits = config.getParameter<int>("MaxHits");
  theVerbosity = config.getParameter<int>("Verbosity");
  theAlignment = config.getParameter<bool>("Alignment");



  if(theVerbosity >=1 )cout << " Inside T1RoadProducer "<<endl;

  if(theAlignment)cout << "ALIGNMENT: applying loaded alignment." << endl;
  if(!theAlignment)cout << "T1RoadProducer: --NO EXTERNAL ALIGNMENT--" << endl;
}

T1RoadProducer::~T1RoadProducer(){
  delete theT1Geometry;
  //delete RecHitMatrix;
//  if(file!=NULL)
//  fclose(file);
}

void T1RoadProducer::produce(Event& event, const EventSetup& setup) {
  //  cout << "Inside T1RoadProducer::produce " << endl;
  // Get the reco hits from the event
  if(theVerbosity >=2 ){
    std::cout <<"Evt "<<event.id().event() <<std::endl;
  }else{
  if(event.id().event()%1000==0)std::cout <<"Evt "<<event.id().event() <<std::endl;
  }

  edm::Handle<T1RecHit2DCollection>myRecoColl;
  event.getByLabel(t1RecHit2DCollectionLabel, myRecoColl);

  vector<pair<T1RecHitGlobal,int> > * global_copy_of_myRecoColl = new vector<pair<T1RecHitGlobal,int> >;

  T1RecHit2DCollection::const_iterator T1Reco_it;
  vector<pair<T1RecHitGlobal,int> >::iterator T1RecoGlobal_it;
  vector<pair<T1RecHitGlobal,int> >::iterator T1RecoGlobal_it2;
  vector<T1RecHitGlobal>::iterator T1RecoGlobal_it3;
  vector<pair<T1RecHitGlobal,int> >::iterator T1RecoGlobal_it4;

  // Create the pointer to the collection which will store the roads
  auto_ptr<T1RoadCollection> roadCollection(new T1RoadCollection());

  T1ChamberSpecs parametri;


  int numHits=0;

  for(T1Reco_it=myRecoColl->begin();T1Reco_it!=myRecoColl->end();T1Reco_it++){
 
    (*T1Reco_it).t1DetId().Arm();
    (*T1Reco_it).t1DetId().Plane();
    (*T1Reco_it).t1DetId().CSC();

    //store rec hits in temporary collection
    if(theVerbosity >=2){
      std::cout << "                                                             " << ++numHits<<std::endl;
    }

    float xxxR=0;
    float yyyR=0;
    float zzzR=0;

    theT1Geometry->Local2Global((*T1Reco_it).t1DetId().rawId(),(*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y(),xxxR,yyyR,zzzR);

    float xxxxR=xxxR;
    float yyyyR=yyyR;
    float zzzzR=zzzR;
    if(theAlignment)
    theT1Geometry->Align((*T1Reco_it).t1DetId().rawId(),xxxR,yyyR,xxxxR,yyyyR,zzzzR);

    float errorX =0;
    float errorY =0;
    float errorZ = 25.;

    theT1Geometry->RotationLocalError2Global((*T1Reco_it).t1DetId().rawId(),(*T1Reco_it).localPositionError().xx(),(*T1Reco_it).localPositionError().yy(),errorX,errorY);
    GlobalError theError (errorX,0.,errorY,0.,0.,errorZ);


    GlobalPoint thePos( xxxxR, yyyyR, zzzzR );
	  
    //	  T1RecHitGlobal *myGlobalHit = new T1RecHitGlobal(thePos, theError);
    T1RecHitGlobal myGlobalHit((*T1Reco_it).t1DetId(),thePos, theError);
    //	  pair<T1RecHitGlobal, int> *myGlobalHitPair=new pair<T1RecHitGlobal, int>(myGlobalHit,0);
    //	  global_copy_of_myRecoColl->push_back(*myGlobalHitPair);
   
    global_copy_of_myRecoColl->push_back(pair<T1RecHitGlobal, int>(myGlobalHit,0) );

  }

if(theVerbosity >=2){
  for(T1RecoGlobal_it=global_copy_of_myRecoColl->begin();T1RecoGlobal_it!=global_copy_of_myRecoColl->end();T1RecoGlobal_it++){
    std::cout << (*T1RecoGlobal_it).first <<std::endl;
  }
}

  for(T1RecoGlobal_it=global_copy_of_myRecoColl->begin();T1RecoGlobal_it!=global_copy_of_myRecoColl->end();T1RecoGlobal_it++){
 
    //    T1DetId t1Id((*T1Reco_it).t1DetId());
    if((*T1RecoGlobal_it).second == 0){
      T1Road *road = new T1Road();

      //    T1Road road;

      pair<T1RecHitGlobal,int> tempRecHitPair( (*T1RecoGlobal_it) );

      float tempX = (*T1RecoGlobal_it).first.GlobalPosition().x();
      float tempY = (*T1RecoGlobal_it).first.GlobalPosition().y();
      float tempZ = (*T1RecoGlobal_it).first.GlobalPosition().z();

      for(T1RecoGlobal_it2=global_copy_of_myRecoColl->begin();T1RecoGlobal_it2!=global_copy_of_myRecoColl->end();T1RecoGlobal_it2++){
	if((*T1RecoGlobal_it2).second == 0){

if(theVerbosity >=2){
	  cout <<"Reco hit: "<< tempX << " " << tempY << " " << tempZ <<endl;
	  //    if(event.id().event()==10)
	  //	cout << (*T1RecoGlobal_it2).first.GlobalPosition().x() << " " << (*T1RecoGlobal_it2).first.GlobalPosition().y() << " " << (*T1RecoGlobal_it2).first.GlobalPosition().z() <<endl;
}
	  float DE = Eta( tempX,tempY,tempZ ) - Eta( (*T1RecoGlobal_it2).first.GlobalPosition().x(),(*T1RecoGlobal_it2).first.GlobalPosition().y(),(*T1RecoGlobal_it2).first.GlobalPosition().z() );
	  float DF = fabs( Phi( tempX,tempY ) - Phi( (*T1RecoGlobal_it2).first.GlobalPosition().x(),(*T1RecoGlobal_it2).first.GlobalPosition().y() ) );
	  if(DF>3.14159)DF=2*3.14159-DF;
if(theVerbosity >=2){
	  cout << "DE= " << DE << "     DF= "<< DF<<endl;
}

	  if(fabs(DE) < theDeltaEta && fabs(DF)< theDeltaPhi ){
	    //	  T1RecHitGlobal *tempRecHit2 = new T1RecHitGlobal( (*T1RecoGlobal_it2).first );
	    //       	  road->push_back(*tempRecHit2);


	    road->push_back(T1RecHitGlobal( (*T1RecoGlobal_it2).first ));


	    if(DE != 0 && DF != 0){
	      //	  T1RecoGlobal_it2 = global_copy_of_myRecoColl->erase(T1RecoGlobal_it2);
	  
	      //  T1RecoGlobal_it2--;
	    }
	  }
	}
      }
      if(road->size()>=(unsigned)theMinHits && road->size()<=(unsigned)theMaxHits){

	roadCollection->push_back(*road);
      
	for(T1RecoGlobal_it3=road->begin();T1RecoGlobal_it3!=road->end();T1RecoGlobal_it3++){
	  for(T1RecoGlobal_it4=global_copy_of_myRecoColl->begin();T1RecoGlobal_it4!=global_copy_of_myRecoColl->end();T1RecoGlobal_it4++){

	    //cout << (*T1RecoGlobal_it3)<<endl;
	    //    cout << (*T1RecoGlobal_it4).first<<endl;
	    //
	    //      cout << endl;
	    if((*T1RecoGlobal_it3).GlobalPosition().x()==(*T1RecoGlobal_it4).first.GlobalPosition().x() && (*T1RecoGlobal_it3).GlobalPosition().y()==(*T1RecoGlobal_it4).first.GlobalPosition().y() && (*T1RecoGlobal_it3).GlobalPosition().z()==(*T1RecoGlobal_it4).first.GlobalPosition().z() ){
	      (*T1RecoGlobal_it4).second = 1 ;
	      //	      cout << "                                                  Messo a 1 " <<endl;
	    }
	  }
	}
	delete road;
      }
      else{

	delete road;
      
      }
    }
  }
if(theVerbosity >=2){
  std::cout << "Road Collection size (taglia): " << roadCollection->size()<< std::endl;

  T1RoadCollection::iterator RC_it;
  for(RC_it=roadCollection->begin(); RC_it!=roadCollection->end(); RC_it++){
    std::cout << "Road size: " << (*RC_it).size() << std::endl;
    for(unsigned int iy = 0; iy < (*RC_it).size(); iy++){
      cout << (*RC_it)[iy] <<endl;
    }
  }
}
  event.put(roadCollection);

  delete global_copy_of_myRecoColl;
}

float T1RoadProducer::Eta(float x,float y,float z){
  float xyt;
  float c=0;
  float eta2;
  xyt = sqrt(x*x + y*y);
  //theta
  if(z>0) c = atan(xyt/z);
  if(z<0) c = atan(xyt/z)+3.14159;
  if(z==0) {c = 3.14159;}
  //pseudorapidity
  eta2 = -log(tan(c/2.));
  return eta2;
}
float T1RoadProducer::Phi(float x,float y){
  float c=0;
  if(x>0 && y>0) c = atan(y/x);
  if(x<0) c = atan(y/x)+3.14159;
  if(x>0 && y<0) c = atan(y/x)+6.28318;
  return c;
}
