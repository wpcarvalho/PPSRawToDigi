/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
*/
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TotemAlignment/T1Alignment/interface/T1Sex.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TVector3.h"
#include "TVector2.h"
#include <fstream>
#include <sstream>


#ifndef PI
#define PI 3.14159
#endif
#define _30DEG 0.5236
#define _25DEG 0.4363
#define _60DEG 1.0472
#define _90DEG 1.5708
#define _120DEG 2.0944
#define _150DEG 2.618
#define _210DEG 3.665
#define _240DEG 4.1888
#define _270DEG 4.7124
#define _300DEG 5.2360
#define _330DEG 5.7596
#define _3DEG 0.05239878

  using namespace edm;
  using namespace std;

//#define _DEBUG_

static const char* ROAD_COLLECTION_LABEL = "roadCollection_Label";


T1Sex::T1Sex(const edm::ParameterSet& iConfig):_Verbosity(0)
{
 
  _Verbosity =  iConfig.getParameter<int>("Verbosity");
  if(!iConfig.exists(ROAD_COLLECTION_LABEL)){
	  cout<<"T1Sex expecting input label: "<<ROAD_COLLECTION_LABEL
			  <<" for data type: T1RoadCollection"<<endl;
  }
  inputLabelRoadCollection_ = iConfig.getParameter<edm::InputTag>(ROAD_COLLECTION_LABEL);

}


T1Sex::~T1Sex()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
T1Sex::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

 


 edm::Handle<T1RoadCollection> roadCollection;
//  iEvent.getByType(roadCollection);
  iEvent.getByLabel(inputLabelRoadCollection_, roadCollection);
if(_Verbosity >= 2){
  std::cout << " Road collection size: " << roadCollection->size() <<std::endl;
}
  T1RoadCollection::const_iterator RC_it;

 
  for(RC_it=roadCollection->begin(); RC_it!=roadCollection->end(); RC_it++){
    


    unsigned int taglia = (*RC_it).size();
    if(taglia > 5 && taglia < 12)
    for(unsigned int i =0; i<taglia-1; i++){
      for(unsigned int j =i+1; j<taglia; j++){
      
	float DDDx = (*RC_it)[i].GlobalPosition().x()-(*RC_it)[j].GlobalPosition().x();
	float DDDy = (*RC_it)[i].GlobalPosition().y()-(*RC_it)[j].GlobalPosition().y();
	float DDDz = (*RC_it)[i].GlobalPosition().z()-(*RC_it)[j].GlobalPosition().z();
      
      
      

	  if ( fabs(DDDz) < 100 && fabs(DDDz) > 0 ){
	    if( sqrt(DDDx*DDDx + DDDy*DDDy) < 50){

	      if( (*RC_it)[i].GlobalPosition().z() > 0){
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_90DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_90DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
		  cout << (*RC_it)[i].GlobalPosition().x() << " " << (*RC_it)[i].GlobalPosition().y() << " " << (*RC_it)[i].GlobalPosition().z() << " "<<endl;
		  cout << (*RC_it)[j].GlobalPosition().x() << " " << (*RC_it)[j].GlobalPosition().y() << " " << (*RC_it)[j].GlobalPosition().z() << " "<<endl;
		  cout << endl;
 
		  hDDDx_0_0_01->Fill(DDDx);
		  hDDDy_0_0_01->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_150DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_150DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_0_0_12->Fill(DDDx);
		  hDDDy_0_0_12->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_210DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_210DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_0_0_23->Fill(DDDx);
		  hDDDy_0_0_23->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_270DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_270DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_0_0_34->Fill(DDDx);
		  hDDDy_0_0_34->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_330DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_330DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_0_0_45->Fill(DDDx);
		  hDDDy_0_0_45->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_30DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_30DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_0_0_50->Fill(DDDx);
		  hDDDy_0_0_50->Fill(DDDy);
		}





		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_90DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_90DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_0_1_01->Fill(DDDx);
		  hDDDy_0_1_01->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_150DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_150DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_0_1_12->Fill(DDDx);
		  hDDDy_0_1_12->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_210DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_210DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_0_1_23->Fill(DDDx);
		  hDDDy_0_1_23->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_270DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_270DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_0_1_34->Fill(DDDx);
		  hDDDy_0_1_34->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_330DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_330DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_0_1_45->Fill(DDDx);
		  hDDDy_0_1_45->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_30DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_30DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_0_1_50->Fill(DDDx);
		  hDDDy_0_1_50->Fill(DDDy);
		}









		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_90DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_90DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_0_2_01->Fill(DDDx);
		  hDDDy_0_2_01->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_150DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_150DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_0_2_12->Fill(DDDx);
		  hDDDy_0_2_12->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_210DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_210DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_0_2_23->Fill(DDDx);
		  hDDDy_0_2_23->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_270DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_270DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_0_2_34->Fill(DDDx);
		  hDDDy_0_2_34->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_330DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_330DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_0_2_45->Fill(DDDx);
		  hDDDy_0_2_45->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_30DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_30DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_0_2_50->Fill(DDDx);
		  hDDDy_0_2_50->Fill(DDDy);
		}









		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_90DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_90DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_0_3_01->Fill(DDDx);
		  hDDDy_0_3_01->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_150DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_150DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_0_3_12->Fill(DDDx);
		  hDDDy_0_3_12->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_210DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_210DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_0_3_23->Fill(DDDx);
		  hDDDy_0_3_23->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_270DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_270DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_0_3_34->Fill(DDDx);
		  hDDDy_0_3_34->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_330DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_330DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_0_3_45->Fill(DDDx);
		  hDDDy_0_3_45->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_30DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_30DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_0_3_50->Fill(DDDx);
		  hDDDy_0_3_50->Fill(DDDy);
		}





		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_90DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_90DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_0_4_01->Fill(DDDx);
		  hDDDy_0_4_01->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_150DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_150DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_0_4_12->Fill(DDDx);
		  hDDDy_0_4_12->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_210DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_210DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_0_4_23->Fill(DDDx);
		  hDDDy_0_4_23->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_270DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_270DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_0_4_34->Fill(DDDx);
		  hDDDy_0_4_34->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_330DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_330DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_0_4_45->Fill(DDDx);
		  hDDDy_0_4_45->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_30DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_30DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_0_4_50->Fill(DDDx);
		  hDDDy_0_4_50->Fill(DDDy);
		}







	      }



	      if( (*RC_it)[i].GlobalPosition().z() < 0){
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_90DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_90DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_1_0_01->Fill(DDDx);
		  hDDDy_1_0_01->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_150DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_150DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_1_0_12->Fill(DDDx);
		  hDDDy_1_0_12->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_210DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_210DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_1_0_23->Fill(DDDx);
		  hDDDy_1_0_23->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_270DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_270DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_1_0_34->Fill(DDDx);
		  hDDDy_1_0_34->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_330DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_330DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_1_0_45->Fill(DDDx);
		  hDDDy_1_0_45->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_30DEG-_3DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_30DEG-_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 0 ){
	      
		  hDDDx_1_0_50->Fill(DDDx);
		  hDDDy_1_0_50->Fill(DDDy);
		}





		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_90DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_90DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_1_1_01->Fill(DDDx);
		  hDDDy_1_1_01->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_150DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_150DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_1_1_12->Fill(DDDx);
		  hDDDy_1_1_12->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_210DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_210DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_1_1_23->Fill(DDDx);
		  hDDDy_1_1_23->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_270DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_270DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_1_1_34->Fill(DDDx);
		  hDDDy_1_1_34->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_330DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_330DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_1_1_45->Fill(DDDx);
		  hDDDy_1_1_45->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_30DEG-_3DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_30DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 1 ){
	      
		  hDDDx_1_1_50->Fill(DDDx);
		  hDDDy_1_1_50->Fill(DDDy);
		}









		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_90DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_90DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_1_2_01->Fill(DDDx);
		  hDDDy_1_2_01->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_150DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_150DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_1_2_12->Fill(DDDx);
		  hDDDy_1_2_12->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_210DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_210DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_1_2_23->Fill(DDDx);
		  hDDDy_1_2_23->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_270DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_270DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_1_2_34->Fill(DDDx);
		  hDDDy_1_2_34->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_330DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_330DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_1_2_45->Fill(DDDx);
		  hDDDy_1_2_45->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_30DEG+_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_30DEG+_3DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 2 ){
	      
		  hDDDx_1_2_50->Fill(DDDx);
		  hDDDy_1_2_50->Fill(DDDy);
		}









		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_90DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_90DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_1_3_01->Fill(DDDx);
		  hDDDy_1_3_01->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_150DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_150DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_1_3_12->Fill(DDDx);
		  hDDDy_1_3_12->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_210DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_210DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_1_3_23->Fill(DDDx);
		  hDDDy_1_3_23->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_270DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_270DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_1_3_34->Fill(DDDx);
		  hDDDy_1_3_34->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_330DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_330DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_1_3_45->Fill(DDDx);
		  hDDDy_1_3_45->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_30DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_30DEG+_3DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 3 ){
	      
		  hDDDx_1_3_50->Fill(DDDx);
		  hDDDy_1_3_50->Fill(DDDy);
		}





		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_90DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_90DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_1_4_01->Fill(DDDx);
		  hDDDy_1_4_01->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_150DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_150DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_1_4_12->Fill(DDDx);
		  hDDDy_1_4_12->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_210DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_210DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_1_4_23->Fill(DDDx);
		  hDDDy_1_4_23->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_270DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_270DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_1_4_34->Fill(DDDx);
		  hDDDy_1_4_34->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_330DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_330DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_1_4_45->Fill(DDDx);
		  hDDDy_1_4_45->Fill(DDDy);
		}
		if(Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) > (_30DEG-_3DEG) && 
		   Phi((*RC_it)[i].GlobalPosition().x(),(*RC_it)[i].GlobalPosition().y()) < (_30DEG+_3DEG) &&
		   plane((*RC_it)[i].GlobalPosition().z() ) == 4 ){
	      
		  hDDDx_1_4_50->Fill(DDDx);
		  hDDDy_1_4_50->Fill(DDDy);
		}





	      }
	    }
	  }
      }
    

    
    }
  }





}


// ------------ method called once each job just before starting event loop  ------------
void 
T1Sex::beginJob()
{
  if(_Verbosity>0)
    cout << "T1Sex::beginJob" << endl;




  hDDDx_0_0_01 = new TH1D("DeltaX_0_0_01","DeltaX_0_0_01",100,-50.,50.);
  hDDDy_0_0_01 = new TH1D("DeltaY_0_0_01","DeltaY_0_0_01",100,-50.,50.);

  hDDDx_0_0_12 = new TH1D("DeltaX_0_0_12","DeltaX_0_0_12",100,-50.,50.);
  hDDDy_0_0_12 = new TH1D("DeltaY_0_0_12","DeltaY_0_0_12",100,-50.,50.);

  hDDDx_0_0_23 = new TH1D("DeltaX_0_0_23","DeltaX_0_0_23",100,-50.,50.);
  hDDDy_0_0_23 = new TH1D("DeltaY_0_0_23","DeltaY_0_0_23",100,-50.,50.);

  hDDDx_0_0_34 = new TH1D("DeltaX_0_0_34","DeltaX_0_0_34",100,-50.,50.);
  hDDDy_0_0_34 = new TH1D("DeltaY_0_0_34","DeltaY_0_0_34",100,-50.,50.);

  hDDDx_0_0_45 = new TH1D("DeltaX_0_0_45","DeltaX_0_0_45",100,-50.,50.);
  hDDDy_0_0_45 = new TH1D("DeltaY_0_0_45","DeltaY_0_0_45",100,-50.,50.);

  hDDDx_0_0_50 = new TH1D("DeltaX_0_0_50","DeltaX_0_0_50",100,-50.,50.);
  hDDDy_0_0_50 = new TH1D("DeltaY_0_0_50","DeltaY_0_0_50",100,-50.,50.);





  hDDDx_0_1_01 = new TH1D("DeltaX_0_1_01","DeltaX_0_1_01",100,-50.,50.);
  hDDDy_0_1_01 = new TH1D("DeltaY_0_1_01","DeltaY_0_1_01",100,-50.,50.);

  hDDDx_0_1_12 = new TH1D("DeltaX_0_1_12","DeltaX_0_1_12",100,-50.,50.);
  hDDDy_0_1_12 = new TH1D("DeltaY_0_1_12","DeltaY_0_1_12",100,-50.,50.);

  hDDDx_0_1_23 = new TH1D("DeltaX_0_1_23","DeltaX_0_1_23",100,-50.,50.);
  hDDDy_0_1_23 = new TH1D("DeltaY_0_1_23","DeltaY_0_1_23",100,-50.,50.);

  hDDDx_0_1_34 = new TH1D("DeltaX_0_1_34","DeltaX_0_1_34",100,-50.,50.);
  hDDDy_0_1_34 = new TH1D("DeltaY_0_1_34","DeltaY_0_1_34",100,-50.,50.);

  hDDDx_0_1_45 = new TH1D("DeltaX_0_1_45","DeltaX_0_1_45",100,-50.,50.);
  hDDDy_0_1_45 = new TH1D("DeltaY_0_1_45","DeltaY_0_1_45",100,-50.,50.);

  hDDDx_0_1_50 = new TH1D("DeltaX_0_1_50","DeltaX_0_1_50",100,-50.,50.);
  hDDDy_0_1_50 = new TH1D("DeltaY_0_1_50","DeltaY_0_1_50",100,-50.,50.);



  hDDDx_0_2_01 = new TH1D("DeltaX_0_2_01","DeltaX_0_2_01",100,-50.,50.);
  hDDDy_0_2_01 = new TH1D("DeltaY_0_2_01","DeltaY_0_2_01",100,-50.,50.);

  hDDDx_0_2_12 = new TH1D("DeltaX_0_2_12","DeltaX_0_2_12",100,-50.,50.);
  hDDDy_0_2_12 = new TH1D("DeltaY_0_2_12","DeltaY_0_2_12",100,-50.,50.);

  hDDDx_0_2_23 = new TH1D("DeltaX_0_2_23","DeltaX_0_2_23",100,-50.,50.);
  hDDDy_0_2_23 = new TH1D("DeltaY_0_2_23","DeltaY_0_2_23",100,-50.,50.);

  hDDDx_0_2_34 = new TH1D("DeltaX_0_2_34","DeltaX_0_2_34",100,-50.,50.);
  hDDDy_0_2_34 = new TH1D("DeltaY_0_2_34","DeltaY_0_2_34",100,-50.,50.);

  hDDDx_0_2_45 = new TH1D("DeltaX_0_2_45","DeltaX_0_2_45",100,-50.,50.);
  hDDDy_0_2_45 = new TH1D("DeltaY_0_2_45","DeltaY_0_2_45",100,-50.,50.);

  hDDDx_0_2_50 = new TH1D("DeltaX_0_2_50","DeltaX_0_2_50",100,-50.,50.);
  hDDDy_0_2_50 = new TH1D("DeltaY_0_2_50","DeltaY_0_2_50",100,-50.,50.);





  hDDDx_0_3_01 = new TH1D("DeltaX_0_3_01","DeltaX_0_3_01",100,-50.,50.);
  hDDDy_0_3_01 = new TH1D("DeltaY_0_3_01","DeltaY_0_3_01",100,-50.,50.);

  hDDDx_0_3_12 = new TH1D("DeltaX_0_3_12","DeltaX_0_3_12",100,-50.,50.);
  hDDDy_0_3_12 = new TH1D("DeltaY_0_3_12","DeltaY_0_3_12",100,-50.,50.);

  hDDDx_0_3_23 = new TH1D("DeltaX_0_3_23","DeltaX_0_3_23",100,-50.,50.);
  hDDDy_0_3_23 = new TH1D("DeltaY_0_3_23","DeltaY_0_3_23",100,-50.,50.);

  hDDDx_0_3_34 = new TH1D("DeltaX_0_3_34","DeltaX_0_3_34",100,-50.,50.);
  hDDDy_0_3_34 = new TH1D("DeltaY_0_3_34","DeltaY_0_3_34",100,-50.,50.);

  hDDDx_0_3_45 = new TH1D("DeltaX_0_3_45","DeltaX_0_3_45",100,-50.,50.);
  hDDDy_0_3_45 = new TH1D("DeltaY_0_3_45","DeltaY_0_3_45",100,-50.,50.);

  hDDDx_0_3_50 = new TH1D("DeltaX_0_3_50","DeltaX_0_3_50",100,-50.,50.);
  hDDDy_0_3_50 = new TH1D("DeltaY_0_3_50","DeltaY_0_3_50",100,-50.,50.);





  hDDDx_0_4_01 = new TH1D("DeltaX_0_4_01","DeltaX_0_4_01",100,-50.,50.);
  hDDDy_0_4_01 = new TH1D("DeltaY_0_4_01","DeltaY_0_4_01",100,-50.,50.);

  hDDDx_0_4_12 = new TH1D("DeltaX_0_4_12","DeltaX_0_4_12",100,-50.,50.);
  hDDDy_0_4_12 = new TH1D("DeltaY_0_4_12","DeltaY_0_4_12",100,-50.,50.);

  hDDDx_0_4_23 = new TH1D("DeltaX_0_4_23","DeltaX_0_4_23",100,-50.,50.);
  hDDDy_0_4_23 = new TH1D("DeltaY_0_4_23","DeltaY_0_4_23",100,-50.,50.);

  hDDDx_0_4_34 = new TH1D("DeltaX_0_4_34","DeltaX_0_4_34",100,-50.,50.);
  hDDDy_0_4_34 = new TH1D("DeltaY_0_4_34","DeltaY_0_4_34",100,-50.,50.);

  hDDDx_0_4_45 = new TH1D("DeltaX_0_4_45","DeltaX_0_4_45",100,-50.,50.);
  hDDDy_0_4_45 = new TH1D("DeltaY_0_4_45","DeltaY_0_4_45",100,-50.,50.);

  hDDDx_0_4_50 = new TH1D("DeltaX_0_4_50","DeltaX_0_4_50",100,-50.,50.);
  hDDDy_0_4_50 = new TH1D("DeltaY_0_4_50","DeltaY_0_4_50",100,-50.,50.);









  hDDDx_1_0_01 = new TH1D("Deltax_1_0_01","Deltax_1_0_01",100,-50.,50.);
  hDDDy_1_0_01 = new TH1D("Deltay_1_0_01","Deltay_1_0_01",100,-50.,50.);

  hDDDx_1_0_12 = new TH1D("Deltax_1_0_12","Deltax_1_0_12",100,-50.,50.);
  hDDDy_1_0_12 = new TH1D("Deltay_1_0_12","Deltay_1_0_12",100,-50.,50.);

  hDDDx_1_0_23 = new TH1D("Deltax_1_0_23","Deltax_1_0_23",100,-50.,50.);
  hDDDy_1_0_23 = new TH1D("Deltay_1_0_23","Deltay_1_0_23",100,-50.,50.);

  hDDDx_1_0_34 = new TH1D("Deltax_1_0_34","Deltax_1_0_34",100,-50.,50.);
  hDDDy_1_0_34 = new TH1D("Deltay_1_0_34","Deltay_1_0_34",100,-50.,50.);

  hDDDx_1_0_45 = new TH1D("Deltax_1_0_45","Deltax_1_0_45",100,-50.,50.);
  hDDDy_1_0_45 = new TH1D("Deltay_1_0_45","Deltay_1_0_45",100,-50.,50.);

  hDDDx_1_0_50 = new TH1D("Deltax_1_0_50","Deltax_1_0_50",100,-50.,50.);
  hDDDy_1_0_50 = new TH1D("Deltay_1_0_50","Deltay_1_0_50",100,-50.,50.);





  hDDDx_1_1_01 = new TH1D("Deltax_1_1_01","Deltax_1_1_01",100,-50.,50.);
  hDDDy_1_1_01 = new TH1D("Deltay_1_1_01","Deltay_1_1_01",100,-50.,50.);

  hDDDx_1_1_12 = new TH1D("Deltax_1_1_12","Deltax_1_1_12",100,-50.,50.);
  hDDDy_1_1_12 = new TH1D("Deltay_1_1_12","Deltay_1_1_12",100,-50.,50.);

  hDDDx_1_1_23 = new TH1D("Deltax_1_1_23","Deltax_1_1_23",100,-50.,50.);
  hDDDy_1_1_23 = new TH1D("Deltay_1_1_23","Deltay_1_1_23",100,-50.,50.);

  hDDDx_1_1_34 = new TH1D("Deltax_1_1_34","Deltax_1_1_34",100,-50.,50.);
  hDDDy_1_1_34 = new TH1D("Deltay_1_1_34","Deltay_1_1_34",100,-50.,50.);

  hDDDx_1_1_45 = new TH1D("Deltax_1_1_45","Deltax_1_1_45",100,-50.,50.);
  hDDDy_1_1_45 = new TH1D("Deltay_1_1_45","Deltay_1_1_45",100,-50.,50.);

  hDDDx_1_1_50 = new TH1D("Deltax_1_1_50","Deltax_1_1_50",100,-50.,50.);
  hDDDy_1_1_50 = new TH1D("Deltay_1_1_50","Deltay_1_1_50",100,-50.,50.);



  hDDDx_1_2_01 = new TH1D("Deltax_1_2_01","Deltax_1_2_01",100,-50.,50.);
  hDDDy_1_2_01 = new TH1D("Deltay_1_2_01","Deltay_1_2_01",100,-50.,50.);

  hDDDx_1_2_12 = new TH1D("Deltax_1_2_12","Deltax_1_2_12",100,-50.,50.);
  hDDDy_1_2_12 = new TH1D("Deltay_1_2_12","Deltay_1_2_12",100,-50.,50.);

  hDDDx_1_2_23 = new TH1D("Deltax_1_2_23","Deltax_1_2_23",100,-50.,50.);
  hDDDy_1_2_23 = new TH1D("Deltay_1_2_23","Deltay_1_2_23",100,-50.,50.);

  hDDDx_1_2_34 = new TH1D("Deltax_1_2_34","Deltax_1_2_34",100,-50.,50.);
  hDDDy_1_2_34 = new TH1D("Deltay_1_2_34","Deltay_1_2_34",100,-50.,50.);

  hDDDx_1_2_45 = new TH1D("Deltax_1_2_45","Deltax_1_2_45",100,-50.,50.);
  hDDDy_1_2_45 = new TH1D("Deltay_1_2_45","Deltay_1_2_45",100,-50.,50.);

  hDDDx_1_2_50 = new TH1D("Deltax_1_2_50","Deltax_1_2_50",100,-50.,50.);
  hDDDy_1_2_50 = new TH1D("Deltay_1_2_50","Deltay_1_2_50",100,-50.,50.);





  hDDDx_1_3_01 = new TH1D("Deltax_1_3_01","Deltax_1_3_01",100,-50.,50.);
  hDDDy_1_3_01 = new TH1D("Deltay_1_3_01","Deltay_1_3_01",100,-50.,50.);

  hDDDx_1_3_12 = new TH1D("Deltax_1_3_12","Deltax_1_3_12",100,-50.,50.);
  hDDDy_1_3_12 = new TH1D("Deltay_1_3_12","Deltay_1_3_12",100,-50.,50.);

  hDDDx_1_3_23 = new TH1D("Deltax_1_3_23","Deltax_1_3_23",100,-50.,50.);
  hDDDy_1_3_23 = new TH1D("Deltay_1_3_23","Deltay_1_3_23",100,-50.,50.);

  hDDDx_1_3_34 = new TH1D("Deltax_1_3_34","Deltax_1_3_34",100,-50.,50.);
  hDDDy_1_3_34 = new TH1D("Deltay_1_3_34","Deltay_1_3_34",100,-50.,50.);

  hDDDx_1_3_45 = new TH1D("Deltax_1_3_45","Deltax_1_3_45",100,-50.,50.);
  hDDDy_1_3_45 = new TH1D("Deltay_1_3_45","Deltay_1_3_45",100,-50.,50.);

  hDDDx_1_3_50 = new TH1D("Deltax_1_3_50","Deltax_1_3_50",100,-50.,50.);
  hDDDy_1_3_50 = new TH1D("Deltay_1_3_50","Deltay_1_3_50",100,-50.,50.);





  hDDDx_1_4_01 = new TH1D("Deltax_1_4_01","Deltax_1_4_01",100,-50.,50.);
  hDDDy_1_4_01 = new TH1D("Deltay_1_4_01","Deltay_1_4_01",100,-50.,50.);

  hDDDx_1_4_12 = new TH1D("Deltax_1_4_12","Deltax_1_4_12",100,-50.,50.);
  hDDDy_1_4_12 = new TH1D("Deltay_1_4_12","Deltay_1_4_12",100,-50.,50.);

  hDDDx_1_4_23 = new TH1D("Deltax_1_4_23","Deltax_1_4_23",100,-50.,50.);
  hDDDy_1_4_23 = new TH1D("Deltay_1_4_23","Deltay_1_4_23",100,-50.,50.);

  hDDDx_1_4_34 = new TH1D("Deltax_1_4_34","Deltax_1_4_34",100,-50.,50.);
  hDDDy_1_4_34 = new TH1D("Deltay_1_4_34","Deltay_1_4_34",100,-50.,50.);

  hDDDx_1_4_45 = new TH1D("Deltax_1_4_45","Deltax_1_4_45",100,-50.,50.);
  hDDDy_1_4_45 = new TH1D("Deltay_1_4_45","Deltay_1_4_45",100,-50.,50.);

  hDDDx_1_4_50 = new TH1D("Deltax_1_4_50","Deltax_1_4_50",100,-50.,50.);
  hDDDy_1_4_50 = new TH1D("Deltay_1_4_50","Deltay_1_4_50",100,-50.,50.);






}

// ------------ method called once each job just after ending the event loop  ------------
void T1Sex::endJob() {
  if(_Verbosity>0)
    cout << "T1Sex::endJob" << endl;


  theFile = new TFile("sexFile.root","RECREATE");





  hDDDx_0_0_01 ->Write(); // TH1D("DeltaX_0_0_01","DeltaX_0_0_01",100,-50.,50.);
  hDDDy_0_0_01 ->Write(); // TH1D("DeltaY_0_0_01","DeltaY_0_0_01",100,-50.,50.);

  hDDDx_0_0_12 ->Write(); // TH1D("DeltaX_0_0_12","DeltaX_0_0_12",100,-50.,50.);
  hDDDy_0_0_12 ->Write(); // TH1D("DeltaY_0_0_12","DeltaY_0_0_12",100,-50.,50.);

  hDDDx_0_0_23 ->Write(); // TH1D("DeltaX_0_0_23","DeltaX_0_0_23",100,-50.,50.);
  hDDDy_0_0_23 ->Write(); // TH1D("DeltaY_0_0_23","DeltaY_0_0_23",100,-50.,50.);

  hDDDx_0_0_34 ->Write(); // TH1D("DeltaX_0_0_34","DeltaX_0_0_34",100,-50.,50.);
  hDDDy_0_0_34 ->Write(); // TH1D("DeltaY_0_0_34","DeltaY_0_0_34",100,-50.,50.);

  hDDDx_0_0_45 ->Write(); // TH1D("DeltaX_0_0_45","DeltaX_0_0_45",100,-50.,50.);
  hDDDy_0_0_45 ->Write(); // TH1D("DeltaY_0_0_45","DeltaY_0_0_45",100,-50.,50.);

  hDDDx_0_0_50 ->Write(); // TH1D("DeltaX_0_0_50","DeltaX_0_0_50",100,-50.,50.);
  hDDDy_0_0_50 ->Write(); // TH1D("DeltaY_0_0_50","DeltaY_0_0_50",100,-50.,50.);





  hDDDx_0_1_01 ->Write(); // TH1D("DeltaX_0_1_01","DeltaX_0_1_01",100,-50.,50.);
  hDDDy_0_1_01 ->Write(); // TH1D("DeltaY_0_1_01","DeltaY_0_1_01",100,-50.,50.);

  hDDDx_0_1_12 ->Write(); // TH1D("DeltaX_0_1_12","DeltaX_0_1_12",100,-50.,50.);
  hDDDy_0_1_12 ->Write(); // TH1D("DeltaY_0_1_12","DeltaY_0_1_12",100,-50.,50.);

  hDDDx_0_1_23 ->Write(); // TH1D("DeltaX_0_1_23","DeltaX_0_1_23",100,-50.,50.);
  hDDDy_0_1_23 ->Write(); // TH1D("DeltaY_0_1_23","DeltaY_0_1_23",100,-50.,50.);

  hDDDx_0_1_34 ->Write(); // TH1D("DeltaX_0_1_34","DeltaX_0_1_34",100,-50.,50.);
  hDDDy_0_1_34 ->Write(); // TH1D("DeltaY_0_1_34","DeltaY_0_1_34",100,-50.,50.);

  hDDDx_0_1_45 ->Write(); // TH1D("DeltaX_0_1_45","DeltaX_0_1_45",100,-50.,50.);
  hDDDy_0_1_45 ->Write(); // TH1D("DeltaY_0_1_45","DeltaY_0_1_45",100,-50.,50.);

  hDDDx_0_1_50 ->Write(); // TH1D("DeltaX_0_1_50","DeltaX_0_1_50",100,-50.,50.);
  hDDDy_0_1_50 ->Write(); // TH1D("DeltaY_0_1_50","DeltaY_0_1_50",100,-50.,50.);



  hDDDx_0_2_01 ->Write(); // TH1D("DeltaX_0_2_01","DeltaX_0_2_01",100,-50.,50.);
  hDDDy_0_2_01 ->Write(); // TH1D("DeltaY_0_2_01","DeltaY_0_2_01",100,-50.,50.);

  hDDDx_0_2_12 ->Write(); // TH1D("DeltaX_0_2_12","DeltaX_0_2_12",100,-50.,50.);
  hDDDy_0_2_12 ->Write(); // TH1D("DeltaY_0_2_12","DeltaY_0_2_12",100,-50.,50.);

  hDDDx_0_2_23 ->Write(); // TH1D("DeltaX_0_2_23","DeltaX_0_2_23",100,-50.,50.);
  hDDDy_0_2_23 ->Write(); // TH1D("DeltaY_0_2_23","DeltaY_0_2_23",100,-50.,50.);

  hDDDx_0_2_34 ->Write(); // TH1D("DeltaX_0_2_34","DeltaX_0_2_34",100,-50.,50.);
  hDDDy_0_2_34 ->Write(); // TH1D("DeltaY_0_2_34","DeltaY_0_2_34",100,-50.,50.);

  hDDDx_0_2_45 ->Write(); // TH1D("DeltaX_0_2_45","DeltaX_0_2_45",100,-50.,50.);
  hDDDy_0_2_45 ->Write(); // TH1D("DeltaY_0_2_45","DeltaY_0_2_45",100,-50.,50.);

  hDDDx_0_2_50 ->Write(); // TH1D("DeltaX_0_2_50","DeltaX_0_2_50",100,-50.,50.);
  hDDDy_0_2_50 ->Write(); // TH1D("DeltaY_0_2_50","DeltaY_0_2_50",100,-50.,50.);





  hDDDx_0_3_01 ->Write(); // TH1D("DeltaX_0_3_01","DeltaX_0_3_01",100,-50.,50.);
  hDDDy_0_3_01 ->Write(); // TH1D("DeltaY_0_3_01","DeltaY_0_3_01",100,-50.,50.);

  hDDDx_0_3_12 ->Write(); // TH1D("DeltaX_0_3_12","DeltaX_0_3_12",100,-50.,50.);
  hDDDy_0_3_12 ->Write(); // TH1D("DeltaY_0_3_12","DeltaY_0_3_12",100,-50.,50.);

  hDDDx_0_3_23 ->Write(); // TH1D("DeltaX_0_3_23","DeltaX_0_3_23",100,-50.,50.);
  hDDDy_0_3_23 ->Write(); // TH1D("DeltaY_0_3_23","DeltaY_0_3_23",100,-50.,50.);

  hDDDx_0_3_34 ->Write(); // TH1D("DeltaX_0_3_34","DeltaX_0_3_34",100,-50.,50.);
  hDDDy_0_3_34 ->Write(); // TH1D("DeltaY_0_3_34","DeltaY_0_3_34",100,-50.,50.);

  hDDDx_0_3_45 ->Write(); // TH1D("DeltaX_0_3_45","DeltaX_0_3_45",100,-50.,50.);
  hDDDy_0_3_45 ->Write(); // TH1D("DeltaY_0_3_45","DeltaY_0_3_45",100,-50.,50.);

  hDDDx_0_3_50 ->Write(); // TH1D("DeltaX_0_3_50","DeltaX_0_3_50",100,-50.,50.);
  hDDDy_0_3_50 ->Write(); // TH1D("DeltaY_0_3_50","DeltaY_0_3_50",100,-50.,50.);





  hDDDx_0_4_01 ->Write(); // TH1D("DeltaX_0_4_01","DeltaX_0_4_01",100,-50.,50.);
  hDDDy_0_4_01 ->Write(); // TH1D("DeltaY_0_4_01","DeltaY_0_4_01",100,-50.,50.);

  hDDDx_0_4_12 ->Write(); // TH1D("DeltaX_0_4_12","DeltaX_0_4_12",100,-50.,50.);
  hDDDy_0_4_12 ->Write(); // TH1D("DeltaY_0_4_12","DeltaY_0_4_12",100,-50.,50.);

  hDDDx_0_4_23 ->Write(); // TH1D("DeltaX_0_4_23","DeltaX_0_4_23",100,-50.,50.);
  hDDDy_0_4_23 ->Write(); // TH1D("DeltaY_0_4_23","DeltaY_0_4_23",100,-50.,50.);

  hDDDx_0_4_34 ->Write(); // TH1D("DeltaX_0_4_34","DeltaX_0_4_34",100,-50.,50.);
  hDDDy_0_4_34 ->Write(); // TH1D("DeltaY_0_4_34","DeltaY_0_4_34",100,-50.,50.);

  hDDDx_0_4_45 ->Write(); // TH1D("DeltaX_0_4_45","DeltaX_0_4_45",100,-50.,50.);
  hDDDy_0_4_45 ->Write(); // TH1D("DeltaY_0_4_45","DeltaY_0_4_45",100,-50.,50.);

  hDDDx_0_4_50 ->Write(); // TH1D("DeltaX_0_4_50","DeltaX_0_4_50",100,-50.,50.);
  hDDDy_0_4_50 ->Write(); // TH1D("DeltaY_0_4_50","DeltaY_0_4_50",100,-50.,50.);









  hDDDx_1_0_01 ->Write(); // TH1D("Deltax_1_0_01","Deltax_1_0_01",100,-50.,50.);
  hDDDy_1_0_01 ->Write(); // TH1D("Deltay_1_0_01","Deltay_1_0_01",100,-50.,50.);

  hDDDx_1_0_12 ->Write(); // TH1D("Deltax_1_0_12","Deltax_1_0_12",100,-50.,50.);
  hDDDy_1_0_12 ->Write(); // TH1D("Deltay_1_0_12","Deltay_1_0_12",100,-50.,50.);

  hDDDx_1_0_23 ->Write(); // TH1D("Deltax_1_0_23","Deltax_1_0_23",100,-50.,50.);
  hDDDy_1_0_23 ->Write(); // TH1D("Deltay_1_0_23","Deltay_1_0_23",100,-50.,50.);

  hDDDx_1_0_34 ->Write(); // TH1D("Deltax_1_0_34","Deltax_1_0_34",100,-50.,50.);
  hDDDy_1_0_34 ->Write(); // TH1D("Deltay_1_0_34","Deltay_1_0_34",100,-50.,50.);

  hDDDx_1_0_45 ->Write(); // TH1D("Deltax_1_0_45","Deltax_1_0_45",100,-50.,50.);
  hDDDy_1_0_45 ->Write(); // TH1D("Deltay_1_0_45","Deltay_1_0_45",100,-50.,50.);

  hDDDx_1_0_50 ->Write(); // TH1D("Deltax_1_0_50","Deltax_1_0_50",100,-50.,50.);
  hDDDy_1_0_50 ->Write(); // TH1D("Deltay_1_0_50","Deltay_1_0_50",100,-50.,50.);





  hDDDx_1_1_01 ->Write(); // TH1D("Deltax_1_1_01","Deltax_1_1_01",100,-50.,50.);
  hDDDy_1_1_01 ->Write(); // TH1D("Deltay_1_1_01","Deltay_1_1_01",100,-50.,50.);

  hDDDx_1_1_12 ->Write(); // TH1D("Deltax_1_1_12","Deltax_1_1_12",100,-50.,50.);
  hDDDy_1_1_12 ->Write(); // TH1D("Deltay_1_1_12","Deltay_1_1_12",100,-50.,50.);

  hDDDx_1_1_23 ->Write(); // TH1D("Deltax_1_1_23","Deltax_1_1_23",100,-50.,50.);
  hDDDy_1_1_23 ->Write(); // TH1D("Deltay_1_1_23","Deltay_1_1_23",100,-50.,50.);

  hDDDx_1_1_34 ->Write(); // TH1D("Deltax_1_1_34","Deltax_1_1_34",100,-50.,50.);
  hDDDy_1_1_34 ->Write(); // TH1D("Deltay_1_1_34","Deltay_1_1_34",100,-50.,50.);

  hDDDx_1_1_45 ->Write(); // TH1D("Deltax_1_1_45","Deltax_1_1_45",100,-50.,50.);
  hDDDy_1_1_45 ->Write(); // TH1D("Deltay_1_1_45","Deltay_1_1_45",100,-50.,50.);

  hDDDx_1_1_50 ->Write(); // TH1D("Deltax_1_1_50","Deltax_1_1_50",100,-50.,50.);
  hDDDy_1_1_50 ->Write(); // TH1D("Deltay_1_1_50","Deltay_1_1_50",100,-50.,50.);



  hDDDx_1_2_01 ->Write(); // TH1D("Deltax_1_2_01","Deltax_1_2_01",100,-50.,50.);
  hDDDy_1_2_01 ->Write(); // TH1D("Deltay_1_2_01","Deltay_1_2_01",100,-50.,50.);

  hDDDx_1_2_12 ->Write(); // TH1D("Deltax_1_2_12","Deltax_1_2_12",100,-50.,50.);
  hDDDy_1_2_12 ->Write(); // TH1D("Deltay_1_2_12","Deltay_1_2_12",100,-50.,50.);

  hDDDx_1_2_23 ->Write(); // TH1D("Deltax_1_2_23","Deltax_1_2_23",100,-50.,50.);
  hDDDy_1_2_23 ->Write(); // TH1D("Deltay_1_2_23","Deltay_1_2_23",100,-50.,50.);

  hDDDx_1_2_34 ->Write(); // TH1D("Deltax_1_2_34","Deltax_1_2_34",100,-50.,50.);
  hDDDy_1_2_34 ->Write(); // TH1D("Deltay_1_2_34","Deltay_1_2_34",100,-50.,50.);

  hDDDx_1_2_45 ->Write(); // TH1D("Deltax_1_2_45","Deltax_1_2_45",100,-50.,50.);
  hDDDy_1_2_45 ->Write(); // TH1D("Deltay_1_2_45","Deltay_1_2_45",100,-50.,50.);

  hDDDx_1_2_50 ->Write(); // TH1D("Deltax_1_2_50","Deltax_1_2_50",100,-50.,50.);
  hDDDy_1_2_50 ->Write(); // TH1D("Deltay_1_2_50","Deltay_1_2_50",100,-50.,50.);





  hDDDx_1_3_01 ->Write(); // TH1D("Deltax_1_3_01","Deltax_1_3_01",100,-50.,50.);
  hDDDy_1_3_01 ->Write(); // TH1D("Deltay_1_3_01","Deltay_1_3_01",100,-50.,50.);

  hDDDx_1_3_12 ->Write(); // TH1D("Deltax_1_3_12","Deltax_1_3_12",100,-50.,50.);
  hDDDy_1_3_12 ->Write(); // TH1D("Deltay_1_3_12","Deltay_1_3_12",100,-50.,50.);

  hDDDx_1_3_23 ->Write(); // TH1D("Deltax_1_3_23","Deltax_1_3_23",100,-50.,50.);
  hDDDy_1_3_23 ->Write(); // TH1D("Deltay_1_3_23","Deltay_1_3_23",100,-50.,50.);

  hDDDx_1_3_34 ->Write(); // TH1D("Deltax_1_3_34","Deltax_1_3_34",100,-50.,50.);
  hDDDy_1_3_34 ->Write(); // TH1D("Deltay_1_3_34","Deltay_1_3_34",100,-50.,50.);

  hDDDx_1_3_45 ->Write(); // TH1D("Deltax_1_3_45","Deltax_1_3_45",100,-50.,50.);
  hDDDy_1_3_45 ->Write(); // TH1D("Deltay_1_3_45","Deltay_1_3_45",100,-50.,50.);

  hDDDx_1_3_50 ->Write(); // TH1D("Deltax_1_3_50","Deltax_1_3_50",100,-50.,50.);
  hDDDy_1_3_50 ->Write(); // TH1D("Deltay_1_3_50","Deltay_1_3_50",100,-50.,50.);





  hDDDx_1_4_01 ->Write(); // TH1D("Deltax_1_4_01","Deltax_1_4_01",100,-50.,50.);
  hDDDy_1_4_01 ->Write(); // TH1D("Deltay_1_4_01","Deltay_1_4_01",100,-50.,50.);

  hDDDx_1_4_12 ->Write(); // TH1D("Deltax_1_4_12","Deltax_1_4_12",100,-50.,50.);
  hDDDy_1_4_12 ->Write(); // TH1D("Deltay_1_4_12","Deltay_1_4_12",100,-50.,50.);

  hDDDx_1_4_23 ->Write(); // TH1D("Deltax_1_4_23","Deltax_1_4_23",100,-50.,50.);
  hDDDy_1_4_23 ->Write(); // TH1D("Deltay_1_4_23","Deltay_1_4_23",100,-50.,50.);

  hDDDx_1_4_34 ->Write(); // TH1D("Deltax_1_4_34","Deltax_1_4_34",100,-50.,50.);
  hDDDy_1_4_34 ->Write(); // TH1D("Deltay_1_4_34","Deltay_1_4_34",100,-50.,50.);

  hDDDx_1_4_45 ->Write(); // TH1D("Deltax_1_4_45","Deltax_1_4_45",100,-50.,50.);
  hDDDy_1_4_45 ->Write(); // TH1D("Deltay_1_4_45","Deltay_1_4_45",100,-50.,50.);

  hDDDx_1_4_50 ->Write(); // TH1D("Deltax_1_4_50","Deltax_1_4_50",100,-50.,50.);
  hDDDy_1_4_50 ->Write(); //     ("Deltay_1_4_50","Deltay_1_4_50",100,-50.,50.);






 

  







  theFile->Close();
}




float T1Sex::Eta(float x,float y,float z){
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
float T1Sex::Phi(float x,float y){
  float c=0;
  if(x>0 && y>0) c = atan(y/x);
  if(x<0) c = atan(y/x)+3.14159;
  if(x>0 && y<0) c = atan(y/x)+6.28318;
  return c;
}

int T1Sex::sextant(float x, float y, float z){
    int sestante = -1;
    float fi = Phi(x,y);
//  cout << " FI " << fi << endl;
    if(z>0){
      if(plane(z)==4){
	if( fi < _30DEG || fi > (2*PI - _30DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG) && fi < (_60DEG +_30DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG) && fi < (_120DEG +_30DEG) ) sestante = 1;
	if(fi > (PI -_30DEG) && fi < (PI +_30DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG) && fi < (_240DEG +_30DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG) && fi < (_300DEG +_30DEG) ) sestante = 4;
      } 
      if(plane(z)==3){
	if( fi < (_30DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG-_3DEG) && fi < (PI +_30DEG-_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG) ) sestante = 4;
      } 
      if(plane(z)==2){
	if( fi < (_30DEG-_3DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG-_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG-_3DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG-_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG-_3DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG-_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG-_3DEG-_3DEG) && fi < (PI +_30DEG-_3DEG-_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG-_3DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG-_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG-_3DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG-_3DEG) ) sestante = 4;
      } 
      if(plane(z)==1){
	if( fi < (_30DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG+_3DEG) && fi < (PI +_30DEG+_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG) ) sestante = 4;
      } 
      if(plane(z)==0){
	if( fi < (_30DEG+_3DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG+_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG+_3DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG+_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG+_3DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG+_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG+_3DEG+_3DEG) && fi < (PI +_30DEG+_3DEG+_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG+_3DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG+_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG+_3DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG+_3DEG) ) sestante = 4;
      } 
    
    }


    if(z<0){
      if(plane(z)==4){
	if( fi < _30DEG || fi > (2*PI - _30DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG) && fi < (_60DEG +_30DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG) && fi < (_120DEG +_30DEG) ) sestante = 0;
	if(fi > (PI -_30DEG) && fi < (PI +_30DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG) && fi < (_240DEG +_30DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG) && fi < (_300DEG +_30DEG) ) sestante = 3;
      } 
      if(plane(z)==3){
	if( fi < (_30DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG+_3DEG) && fi < (PI +_30DEG+_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG) ) sestante = 3;
      } 
      if(plane(z)==2){
	if( fi < (_30DEG+_3DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG+_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG+_3DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG+_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG+_3DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG+_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG+_3DEG+_3DEG) && fi < (PI +_30DEG+_3DEG+_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG+_3DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG+_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG+_3DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG+_3DEG) ) sestante = 3;
      } 
      if(plane(z)==1){
	if( fi < (_30DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG-_3DEG) && fi < (PI +_30DEG-_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG) ) sestante = 3;
      } 
      if(plane(z)==0){
	if( fi < (_30DEG-_3DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG-_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG-_3DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG-_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG-_3DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG-_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG-_3DEG-_3DEG) && fi < (PI +_30DEG-_3DEG-_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG-_3DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG-_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG-_3DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG-_3DEG) ) sestante = 3;
      } 
    
    }

  
    return sestante;
  }
  int T1Sex::plane(float z) {

    int piano = -1;
    if(fabs(z)<8000){
      piano = 0;
 
    }       
    if(fabs(z)<8400 && fabs(z)>8000){
      piano = 1;
    }
    if(fabs(z)<9000 && fabs(z)>8400){
      piano = 2;
    }
    if(fabs(z)<9500 && fabs(z)>9000){
      piano = 3;
    }
    if(fabs(z)>9500){
      piano = 4;
    }
    return piano;
 	
  }
