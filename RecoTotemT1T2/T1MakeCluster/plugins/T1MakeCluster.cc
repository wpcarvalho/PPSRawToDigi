#include "RecoTotemT1T2/T1MakeCluster/interface/T1MakeCluster.h"



#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"
#include "DataFormats/T1DigiSantiard/interface/T1DigiSantiardCollection.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfatCollection.h"
#include "DataFormats/T1Cluster/interface/T1Cluster.h"
#include "DataFormats/T1Cluster/interface/T1ClusterCollection.h"


#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "DataFormats/T1DetId/interface/T1DetId.h"
//#include "T1Clusterizer/T1MakeCluster/interface/T1Clusterizer.h"

#include <string>


using namespace edm;
using namespace std;

//#define _DEBUG_
////
/*
  #ifndef _weight2
  #define _weight2 1000
  #endif
  #ifndef _weight1
  #define _weight1 500
  #endif
*/
#ifndef STRMAX
#define STRMAX 90
#endif



////

T1MakeCluster::T1MakeCluster(const ParameterSet& config) :
  _electronics(config.getParameter<std::string>("Electronics")),

  _weight1(config.getParameter<double>("WEIGHT1")),
  _weight2(config.getParameter<double>("WEIGHT2")),
  _ActivateDeadChannels(config.getParameter<bool>("ActivateDeadChannels"))
{
  if(_electronics == "Santiard")
  t1DigiSantiardCollectionLabel = config.getParameter<edm::InputTag>("T1DigiSantiardCollectionLabel");
  else if(_electronics == "VFAT")
  t1DigiVfatCollectionLabel = config.getParameter<edm::InputTag>("T1DigiVfatCollectionLabel");
  // Set verbose output

  produces<T1ClusterCollection>();
  std::string thegeometryfile = "Geometry/TotemGeometry/data/T1_data_geometry.dat";
  theT1Geometry = new T1Geometry(thegeometryfile);
  theVerbosity = config.getParameter<int>("Verbosity");

  if(theVerbosity >=1 )
{cout << " Inside T1MakeCluster "<<endl;
  // theChiSquare = config.getParameter<double>("ChiSquare");

  
  if(_electronics == "Santiard")
    std::cout << "###################### Using SANTIARD electronics ########################"<<std::endl;


  if(_electronics == "VFAT")
    std::cout << "###################### Using VFAT electronics ########################"<<std::endl;

  }

}

T1MakeCluster::~T1MakeCluster(){
  delete theT1Geometry;
  //delete RecHitMatrix;
}



void T1MakeCluster::produce(Event& event, const EventSetup& setup) {
  // Pass the EventSetup to the algo
  //  theAlgo->setES(setup);

  if(theVerbosity >=2)
    std::cout <<"Evt "<<event.id().event() <<std::endl;

  // Create the pointer to the collection which will store the rechits
  auto_ptr<T1ClusterCollection> clusterCollection(new T1ClusterCollection());

  T1ChamberSpecs chamberspecs;

  if(_electronics == "Santiard"){
    Handle<T1DigiSantiardCollection> Sdigis;
    event.getByLabel(t1DigiSantiardCollectionLabel, Sdigis);
    T1DigiSantiardCollection::DigiRangeIterator t1sdgIt;

    for (t1sdgIt = Sdigis->begin(); t1sdgIt != Sdigis->end();
	 ++t1sdgIt)

      {
	//         std::cout << t1Id << " ------ " << (*t1sdgIt).first << std::endl;
	//	std::cout << t1Id.rawId()<< std::endl;
	
	const T1DigiSantiardCollection::Range& srange = (*t1sdgIt).second;

	for (T1DigiSantiardCollection::const_iterator digi = srange.first;
	     digi != srange.second;
	     ++digi) {
	  float striscia=0.0;
	  if((*digi).balance()==0) striscia=(float)(*digi).strip()-0.25;
	  if((*digi).balance()==1) striscia=(float)(*digi).strip()+0.25;

	  T1Cluster cl(0,0,striscia,(1./4.*chamberspecs.Pitch()),0.5,(*digi).bx());
	  clusterCollection->insertDigi((*t1sdgIt).first,cl);
	}

      }
  }

  ////'''''''''''''''''''''''''

  if(_electronics == "VFAT"){
    Handle<T1DigiVfatCollection> Sdigis;
    event.getByLabel(t1DigiVfatCollectionLabel, Sdigis);
  
    T1DigiVfatCollection::DigiRangeIterator t1sdgIt;

        
    for (t1sdgIt = Sdigis->begin(); t1sdgIt != Sdigis->end();
	 ++t1sdgIt){
    //    std::cout  << (*t1sdgIt).first << std::endl;
      //	std::cout << t1Id.rawId()<< std::endl;
      assert((*t1sdgIt).first.Layer()==1 || (*t1sdgIt).first.Layer()==2 ); // check that it is a strip ID
      const T1DigiVfatCollection::Range& srange = (*t1sdgIt).second;
      int flag1=0;  
//flag2=0, flag4=0, flag5=0, flag6=0;
    //       int j=0, l=0;
      int st[240];
      int original_st[240];
     
      float peakST[240];
      float peakW[240];
      int peakN=0;
      int firststrip=0;
      int laststrip=0;
      float centerstrip=0.0;
      float sigma=0.0;
      float width = 0.0;
//      float peso=1.4;

      int bxx=0;
      for (int g=0; g<240; g++){
	st[g]=0;
	original_st[g]=0;
   
	peakST[g]=0.0;
	peakW[g]=0.0;
      }

   
      for (T1DigiVfatCollection::const_iterator digi = srange.first;
	   digi != srange.second;
	   ++digi) {
	if(digi->threshold()>0){

	st[digi->strip()]=digi->threshold();
	original_st[digi->strip()]=digi->threshold();
	
	}else{
	
	  st[digi->strip()]=1;
	  original_st[digi->strip()]=1;
      
	}

	bxx=digi->bx();
      }

//   (*t1sdgIt).first id the T1DetId

// turn on dead channels adjacent to good channels already on 
      if(_ActivateDeadChannels)
      for( int hh=1; hh<240; hh++){
	if(_T1deadChannelManager.isChannelDead(   (*t1sdgIt).first, hh) && (original_st[hh-1]>0 || original_st[hh+1]>0 ) ){
      
	  if(original_st[hh-1]>0){
	    st[hh]=original_st[hh-1];
	  }else{
	    st[hh]=original_st[hh+1];
	  }
     
	    if(theVerbosity >=2	 
	       &&   (*t1sdgIt).first.Arm()==1
	       &&   (*t1sdgIt).first.Plane()==3
	       &&   (*t1sdgIt).first.CSC()==5)
	  cout << " Added " << (*t1sdgIt).first << " " << hh << endl;
	}
      }


      for (int jj=0; jj<240; jj++){
	if(jj>0){
	  if(st[jj] > st[jj-1]){
	
	    flag1=1;
	    firststrip=jj;
	  }
	  if((st[jj] < st[jj-1]) && flag1 ==1 ){
	    flag1=0;
	    laststrip=jj-1;
	
	    peakST[peakN]=(laststrip + firststrip)/2.;
	    peakW[peakN]=laststrip-firststrip+1;
	    
	
#ifndef _MORE_THAN_ONE_THRESHOLD_
	    centerstrip=peakST[peakN];
	    width = peakW[peakN];
	    sigma = width/sqrt(12.)*chamberspecs.Pitch();
	    if(theVerbosity >=2   
	       &&   (*t1sdgIt).first.Arm()==1
	       &&   (*t1sdgIt).first.Plane()==3
	       &&   (*t1sdgIt).first.CSC()==5
	       ){
	      std::cout << "    @@@@@@@@@@@@@@@  Cluster @@@@@@@@@@@@@@@@@@@@" << std::endl;
	      std::cout << "    @@ First Strip Cluster   "<<firststrip<<       "  @@" << std::endl;
	      std::cout << "    @@ Last Strip Cluster    "<<laststrip<<       "  @@" << std::endl;
	      std::cout << "    @@ Center Strip Cluster  "<<centerstrip<<"  @@" << std::endl;
	      std::cout << "    @@ Cluster width "<<width<<"  @@" << std::endl;
	      std::cout << "    @@ Sigma Cluster  "<<sigma<<"  @@" << std::endl;
	      
	      
	      std::cout << "    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
	      

	      
	    }
	    T1Cluster temp(firststrip,laststrip,centerstrip,sigma,width,bxx);
	    clusterCollection->insertDigi((*t1sdgIt).first,temp);
#endif
peakN++;
	  }

	}//if(jj>0)
    

      }//for( int jj=0; jj<theT1Geometry -> numberOfStrips(t1detid) + 1l; jj++)
      //
#ifdef _MORE_THAN_ONE_THRESHOLD_

      int peak=0;
      int alfa=0;
      float Smin=0.0, Smax=0.0;
      int flag3=0;
      int strip=0;
      if(peakN > 0){
	for(int ii=0; ii<peakN; ii++){

	  Smin=0.0;
	  Smax=0.0;
	  peak=(int)peakST[ii];
	  flag3=0;
	  for(int kk=0; kk<STRMAX; kk++){
	    if(flag3==1){ 
	      break;
	    }else{

	      Smax=peak+STRMAX;
	      if(Smax >  theT1Geometry -> numberOfStrips((*t1sdgIt).first)) Smax = theT1Geometry -> numberOfStrips((*t1sdgIt).first);
	      strip=(kk+peak);
	      if(st[strip]==0){
		Smax=strip-1;
		flag3=1;
	      }//st[strip]==0
	      if(strip == (int)peakST[ii+1] && (ii+1) < peakN){

		flag3=1;
		if((int)peakW[ii]%2==0 && (int)peakW[ii+1]%2==0){

		  alfa=abs((int)peakST[ii]-(int)peakST[ii+1])- (int)peakW[ii]/2 - (int)peakW[ii+1]/2;
		  Smax=(int)peakST[ii] + (int)peakW[ii]/2 - alfa/2;


		}

		if((int)peakW[ii]%2==0 && (int)peakW[ii+1]%2!=0){

		  alfa=abs((int)peakST[ii]-(int)peakST[ii+1])- (int)peakW[ii]/2 - (int)peakW[ii+1]/2 -1;
		  Smax=(int)peakST[ii] + (int)peakW[ii]/2 - alfa/2;


		}
		if((int)peakW[ii]%2!=0 && (int)peakW[ii+1]%2==0){

		  alfa=abs((int)peakST[ii]-(int)peakST[ii+1])- (int)peakW[ii]/2 - (int)peakW[ii+1]/2;
		  Smax=(int)peakST[ii] + (int)peakW[ii]/2  - alfa/2;


		}
		if((int)peakW[ii]%2!=0 && (int)peakW[ii+1]%2!=0){

		  alfa=abs((int)peakST[ii]-(int)peakST[ii+1])- (int)peakW[ii]/2 - (int)peakW[ii+1]/2 -1;
		  Smax=(int)peakST[ii] + (int)peakW[ii]/2 - alfa/2;


		}

	      }// if(strip == peakST[ii+1])

	    }
	  }//	for(int k=0; k<STRMAX && flag3==0; k++)

	  flag3=0;
	  for(int k=0; k> - STRMAX; k--){

	    if(flag3==1){ break;}else{

	      Smin=peak - STRMAX;
	      if(Smin < 0) Smin = 0;
	      strip=(k+peak);

	      if(st[strip]==0){
		Smin=strip + 1;
		flag3=1;
	      }//st[strip]==0
	      if(strip == (int)peakST[ii - 1] && (ii-1) > 0){

		flag3=1;
		if((int)peakW[ii]%2==0 && (int)peakW[ii-1]%2==0){

		  alfa=abs((int)peakST[ii]-(int)peakST[ii-1])- (int)peakW[ii]/2 - (int)peakW[ii-1]/2;
		  Smin=(int)peakST[ii] - (int)peakW[ii]/2 - alfa/2 + 1;

		}

		if((int)peakW[ii]%2==0 && (int)peakW[ii-1]%2!=0){

		  alfa=abs((int)peakST[ii]-(int)peakST[ii-1])- (int)peakW[ii]/2 - (int)peakW[ii-1]/2 -1;
		  Smin=(int)peakST[ii] - (int)peakW[ii]/2 - alfa/2 + 1;


		}
		if((int)peakW[ii]%2!=0 && (int)peakW[ii-1]%2==0){

		  alfa=abs((int)peakST[ii]-(int)peakST[ii-1])- (int)peakW[ii]/2 - (int)peakW[ii-1]/2;
		  Smin=(int)peakST[ii] - (int)peakW[ii]/2  - alfa/2;


		}
		if((int)peakW[ii]%2!=0 && (int)peakW[ii-1]%2!=0){

		  alfa=abs((int)peakST[ii]-(int)peakST[ii-1])- (int)peakW[ii]/2 - (int)peakW[ii-1]/2 -1;
		  Smin=(int)peakST[ii] - (int)peakW[ii]/2 - alfa/2 ;


		}

	      }// if(strip == peakST[ii-1])

	    }
	  }//	for(int k=0; k< - STRMAX && flag3==0; k--)
	
	  centerstrip=Center((int)Smin,(int)Smax,st);
	  sigma=Sigma((int)Smin,(int)Smax,centerstrip,st)*chamberspecs.Pitch();
	  width = (int)Smax - (int)Smin + 1;
if(theVerbosity >=2){
	  std::cout << "@@@@@@@@@@@@@@@  Cluster @@@@@@@@@@@@@@@@@@@@" << std::endl;
	  std::cout << "@@ First Strip Cluster   "<<Smin<<       "  @@" << std::endl;
	  std::cout << "@@ Last Strip Cluster    "<<Smax<<       "  @@" << std::endl;
	  std::cout << "@@ Center Strip Cluster  "<<centerstrip<<"  @@" << std::endl;
	  std::cout << "@@ Sigma Cluster  "<<sigma<<"  @@" << std::endl;


	  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;



}
	  T1Cluster temp(Smin,Smax,centerstrip,sigma,width,bxx);
	  clusterCollection->insertDigi((*t1sdgIt).first,temp);


	}// for(int ii=0; ii<peakN; ii+=)

      }//if(peakN >0)
#endif
    }
  }
  //metetre nell,evento

  event.put(clusterCollection);
}



//calcola il centro di un cluster la firststrip e laststrip devono essere globali
float T1MakeCluster::Center(int fstrip, int lstrip, int * st){

  float  media=0.0;
  float peso=0.0;



  for(int ii=fstrip; ii<=lstrip; ii++){
	
    if(st[ii] == 1) {
      media += (float)ii*_weight1;
      peso+=_weight1; 
    }
	
    if(st[ii] == 2){
      media += (float)(ii)*_weight2; 
      peso+=_weight2; 
    }
	
	
  }
  return media/peso;


}  

float T1MakeCluster::Sigma(int fstrip, int lstrip, float center, int * st){

  float sigma=0.0;
 

  if(fstrip == lstrip){
    sigma = 0.5/sqrt(3.);
  }else{
    sigma = (float)(lstrip - fstrip + 1)/sqrt(12.);
  }
  return sigma;

}


void T1MakeCluster::beginRun(edm::Run const&, edm::EventSetup const& es){
  t1DataBeginJob(es);


}
void T1MakeCluster::t1DataBeginJob(const edm::EventSetup &es) 
{

  ESHandle<DAQMapping> mapping;
  es.get<TotemDAQMappingRecord>().get(mapping);
	edm::ESHandle<AnalysisMask> analysisMask;
	es.get<TotemDAQMappingRecord> ().get(analysisMask);
	_T1deadChannelManager = T1DeadChannelManager(analysisMask);
}

