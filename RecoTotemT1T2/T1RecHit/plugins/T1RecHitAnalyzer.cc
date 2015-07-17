/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
*/
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include "RecoTotemT1T2/T1RecHit/interface/T1RecHitAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"
#include "DataFormats/T1Cluster/interface/T1Cluster.h"
#include "DataFormats/T1Cluster/interface/T1ClusterCollection.h"

#include <vector>

#include "TFile.h"
#include "TH1.h"

#include <CLHEP/Vector/LorentzVector.h>

#ifndef PI
#define PI 3.141592653589793
#endif 

//#define _DEBUG_

T1RecHitAnalyzer::T1RecHitAnalyzer(const edm::ParameterSet& iConfig)
{
  t1RecHit2DCollectionLabel = iConfig.getParameter<edm::InputTag>("T1RecHit2DCollectionLabel");
  t1DigiWireCollectionLabel = iConfig.getParameter<edm::InputTag>("T1DigiWireCollectionLabel");
  t1ClusterCollectionLabel = iConfig.getParameter<edm::InputTag>("T1ClusterCollectionLabel");

  std::string thegeometryfile = "Geometry/TotemGeometry/data/T1_data_geometry.dat";
  layer = new T1Geometry(thegeometryfile);
 outputFileName = iConfig.getParameter<std::string>("OutputFile");
}


T1RecHitAnalyzer::~T1RecHitAnalyzer()
{
 
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
T1RecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;



  if(iEvent.id().event() < 100 || (iEvent.id().event() %100 == 0) )
  std::cout <<"Evento "<<iEvent.id().event() <<std::endl;


  edm::Handle<T1RecHit2DCollection>myRecoColl;
  iEvent.getByLabel(t1RecHit2DCollectionLabel, myRecoColl);


  T1RecHit2DCollection::const_iterator T1Reco_it;


  
  T1ChamberSpecs parametri;
 

  int numHits=0;
  

  for(T1Reco_it=myRecoColl->begin();T1Reco_it!=myRecoColl->end();T1Reco_it++){
    ++numHits; 
    //store rec hits in temporary collection
#ifdef _DEBUG_
    std::cout << "                                                             " << numHits<<std::endl;
#endif

    if((*T1Reco_it).t1DetId().CSC()==0 && (*T1Reco_it).t1DetId().Plane()==0)
      hCSC000->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==4 &&(*T1Reco_it).t1DetId().Plane()==0 )
      hCSC004->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==5 && (*T1Reco_it).t1DetId().Plane()==0)
      hCSC005->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());

    if((*T1Reco_it).t1DetId().CSC()==0 && (*T1Reco_it).t1DetId().Plane()==1)
      hCSC010->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==4 &&(*T1Reco_it).t1DetId().Plane()==1 )
      hCSC014->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==5 && (*T1Reco_it).t1DetId().Plane()==1)
      hCSC015->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());

    if((*T1Reco_it).t1DetId().CSC()==0 && (*T1Reco_it).t1DetId().Plane()==2)
      hCSC020->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==4 &&(*T1Reco_it).t1DetId().Plane()==2 )
      hCSC024->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==5 && (*T1Reco_it).t1DetId().Plane()==2)
      hCSC025->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());

    if((*T1Reco_it).t1DetId().CSC()==0 && (*T1Reco_it).t1DetId().Plane()==3)
      hCSC030->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==4 &&(*T1Reco_it).t1DetId().Plane()==3 )
      hCSC034->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==5 && (*T1Reco_it).t1DetId().Plane()==3)
      hCSC035->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());

    if((*T1Reco_it).t1DetId().CSC()==0 && (*T1Reco_it).t1DetId().Plane()==4)
      hCSC040->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==4 &&(*T1Reco_it).t1DetId().Plane()==4 )
      hCSC044->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==5 && (*T1Reco_it).t1DetId().Plane()==4)
      hCSC045->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());






    if((*T1Reco_it).t1DetId().CSC()==1 && (*T1Reco_it).t1DetId().Plane()==0)
      hCSC001->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==2 &&(*T1Reco_it).t1DetId().Plane()==0 )
      hCSC002->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==3 && (*T1Reco_it).t1DetId().Plane()==0)
      hCSC003->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());

    if((*T1Reco_it).t1DetId().CSC()==1 && (*T1Reco_it).t1DetId().Plane()==1)
      hCSC011->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==2 &&(*T1Reco_it).t1DetId().Plane()==1 )
      hCSC012->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==3 && (*T1Reco_it).t1DetId().Plane()==1)
      hCSC013->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());

    if((*T1Reco_it).t1DetId().CSC()==1 && (*T1Reco_it).t1DetId().Plane()==2)
      hCSC021->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==2 &&(*T1Reco_it).t1DetId().Plane()==2 )
      hCSC022->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==3 && (*T1Reco_it).t1DetId().Plane()==2)
      hCSC023->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());

    if((*T1Reco_it).t1DetId().CSC()==1 && (*T1Reco_it).t1DetId().Plane()==3)
      hCSC031->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==2 &&(*T1Reco_it).t1DetId().Plane()==3 )
      hCSC032->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==3 && (*T1Reco_it).t1DetId().Plane()==3)
      hCSC033->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());

    if((*T1Reco_it).t1DetId().CSC()==1 && (*T1Reco_it).t1DetId().Plane()==4)
      hCSC041->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==2 &&(*T1Reco_it).t1DetId().Plane()==4 )
      hCSC042->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());
    if((*T1Reco_it).t1DetId().CSC()==3 && (*T1Reco_it).t1DetId().Plane()==4)
      hCSC043->Fill((*T1Reco_it).localPosition().x(),(*T1Reco_it).localPosition().y());




//  if((*T1Reco_it).localPosition().y()< -240.)std::cout << " READING : y = " << (*T1Reco_it).localPosition().y() << " in CSC " << (*T1Reco_it).t1DetId().CSC() << std::endl;

  }

//=================================================================================



  // Get the digis from the event
  Handle<T1DigiWireCollection> Wdigis;
  Handle<T1ClusterCollection> Sdigis;

  iEvent.getByLabel(t1DigiWireCollectionLabel, Wdigis);
  iEvent.getByLabel(t1ClusterCollectionLabel, Sdigis);
  
  std::vector< edm::OwnVector<T1RecHit2D> > RecHitMatrix;

  // Create the pointer to the collection which will store the rechits

  // load the geometry
  //  T1Geometry * layer = new T1Geometry;

  // Iterate through all digi collections ordered by LayerId   
  T1DigiWireCollection::DigiRangeIterator t1wdgIt;

  T1ClusterCollection::DigiRangeIterator t1sdgIt;
  T1ClusterCollection::DigiRangeIterator t1sdgIt2;

    
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
		  
		  float yW = layer->yOfWire(t1Id.rawId(), (*digiIt).wire() );
//		  float strisciaA=0.0;

		      
		  float yS = layer->yOfStripCrossing(t1Id.rawId(),(*sdigiIt).Center(),(*sdigiIt2).Center() );
//		  float xS = layer->xOfStripCrossing(t1Id.rawId(),(*sdigiIt).Center(),(*sdigiIt2).Center() );

		  float sigma_yW = (1.5);// /RADQ3;
		  //  float sigma_yS = 2./4. * parametri.Pitch();// / RADQ3;
		  //	  float sigma_xS = (2./4. * parametri.Pitch() / RADQ3); // /RADQ3;
		  float sigma_yS = (*sdigiIt).Sigma();// / RADQ3;
//		  float sigma_xS = (*sdigiIt).Sigma() / RADQ3; // /RADQ3;
		  // float sigma_xSW = 1.5 + 1./4. * parametri.Pitch() / RADQ3;
		      
		  // float mean_x = (xS/sigma_xS/sigma_xS + xSAW/sigma_xSW/sigma_xSW + xSBW/sigma_xSW/sigma_xSW )/ (2./sigma_xSW/sigma_xSW + 1./sigma_xS/sigma_xS);
//		  float mean_x=xS;
//		  float weighted_mean_y = (yW/sigma_yW/sigma_yW + yS/sigma_yS/sigma_yS)/(1./sigma_yW/sigma_yW + 1./sigma_yS/sigma_yS);
//		  float mean_y = yS;
		      
		  // float sigma_x = sqrt( 1./( 1./sigma_xS/sigma_xS + 2./sigma_xSW/sigma_xSW ) );
//		  float sigma_x=sigma_xS;
//		  float sigma_y =  sqrt( 1./(1./sigma_yW/sigma_yW + 1./sigma_yS/sigma_yS) );
		  //   std::cout << "###################### FLAG 4 ########################"<<std::endl;
		  //			    float chi_squared =sqrt( (mean_x - xS)*(mean_x - xS)/sigma_xS/sigma_xS + (mean_x - xSAW)*(mean_x - xSAW)/sigma_xSW/sigma_xSW + (mean_x-xSBW)*(mean_x-xSBW)/sigma_xSW/sigma_xSW + (mean_y - yW)*(mean_y - yW)/sigma_yW/sigma_yW + (mean_y - yS)*(mean_y - yS)/sigma_yS/sigma_yS );
		      
		  float lamda = fabs(yW-yS)/sqrt((sigma_yW)*(sigma_yW)+(sigma_yS)*(sigma_yS));
		      
		  //	    std::cout << "Mean X " << mean_x << "   Mean Y " << mean_y << "   chi "<<chi_squared<< std::endl;
//if(theVerbosity >=2){		 
//  std::cout << (*t1sdgIt2).first.Arm() << " " << (*t1sdgIt2).first.Plane()<< " " << (*t1sdgIt2).first.CSC()<< "XS " <<xS << "  YS " << yS << "  YW "<<yW << std::endl;
//		  std::cout << "Mean X " << mean_x << "   Mean Y " << mean_y <<" Lamda  " << lamda << std::endl;
//}
		 
		  // inserire un controllo sulla geometria !!!
		  // gli hits ricostruiti stanno dentro alla camera ??? 
		      
		      
		  //con lamda verifichiamo la compatibilita delle misure
		  if( lamda <= 5 && sigma_yS < 10 ) {


		    

    if(t1Id.CSC()==0 && t1Id.Plane()==0)
      hCSCcsc000->Fill(yW-yS);
    if(t1Id.CSC()==4 &&t1Id.Plane()==0 )
      hCSCcsc004->Fill(yW-yS);
    if(t1Id.CSC()==5 && t1Id.Plane()==0)
      hCSCcsc005->Fill(yW-yS);

    if(t1Id.CSC()==0 && t1Id.Plane()==1)
      hCSCcsc010->Fill(yW-yS);
    if(t1Id.CSC()==4 &&t1Id.Plane()==1 )
      hCSCcsc014->Fill(yW-yS);
    if(t1Id.CSC()==5 && t1Id.Plane()==1)
      hCSCcsc015->Fill(yW-yS);

    if(t1Id.CSC()==0 && t1Id.Plane()==2)
      hCSCcsc020->Fill(yW-yS);
    if(t1Id.CSC()==4 &&t1Id.Plane()==2 )
      hCSCcsc024->Fill(yW-yS);
    if(t1Id.CSC()==5 && t1Id.Plane()==2)
      hCSCcsc025->Fill(yW-yS);

    if(t1Id.CSC()==0 && t1Id.Plane()==3)
      hCSCcsc030->Fill(yW-yS);
    if(t1Id.CSC()==4 &&t1Id.Plane()==3 )
      hCSCcsc034->Fill(yW-yS);
    if(t1Id.CSC()==5 && t1Id.Plane()==3)
      hCSCcsc035->Fill(yW-yS);

    if(t1Id.CSC()==0 && t1Id.Plane()==4)
      hCSCcsc040->Fill(yW-yS);
    if(t1Id.CSC()==4 &&t1Id.Plane()==4 )
      hCSCcsc044->Fill(yW-yS);
    if(t1Id.CSC()==5 && t1Id.Plane()==4)
      hCSCcsc045->Fill(yW-yS);






    if(t1Id.CSC()==1 && t1Id.Plane()==0)
      hCSCcsc001->Fill(yW-yS);
    if(t1Id.CSC()==2 &&t1Id.Plane()==0 )
      hCSCcsc002->Fill(yW-yS);
    if(t1Id.CSC()==3 && t1Id.Plane()==0)
      hCSCcsc003->Fill(yW-yS);

    if(t1Id.CSC()==1 && t1Id.Plane()==1)
      hCSCcsc011->Fill(yW-yS);
    if(t1Id.CSC()==2 &&t1Id.Plane()==1 )
      hCSCcsc012->Fill(yW-yS);
    if(t1Id.CSC()==3 && t1Id.Plane()==1)
      hCSCcsc013->Fill(yW-yS);

    if(t1Id.CSC()==1 && t1Id.Plane()==2)
      hCSCcsc021->Fill(yW-yS);
    if(t1Id.CSC()==2 &&t1Id.Plane()==2 )
      hCSCcsc022->Fill(yW-yS);
    if(t1Id.CSC()==3 && t1Id.Plane()==2)
      hCSCcsc023->Fill(yW-yS);

    if(t1Id.CSC()==1 && t1Id.Plane()==3)
      hCSCcsc031->Fill(yW-yS);
    if(t1Id.CSC()==2 &&t1Id.Plane()==3 )
      hCSCcsc032->Fill(yW-yS);
    if(t1Id.CSC()==3 && t1Id.Plane()==3)
      hCSCcsc033->Fill(yW-yS);

    if(t1Id.CSC()==1 && t1Id.Plane()==4)
      hCSCcsc041->Fill(yW-yS);
    if(t1Id.CSC()==2 &&t1Id.Plane()==4 )
      hCSCcsc042->Fill(yW-yS);
    if(t1Id.CSC()==3 && t1Id.Plane()==4)
      hCSCcsc043->Fill(yW-yS);





		    
		  }


 		}
		  
	      }
	    }
	      
	  }	
	}
      }
    
    }
  }// loop DigiWire


//=================================================================================

}


// ------------ method called once each job just before starting event loop  ------------
void 
T1RecHitAnalyzer::beginJob()
{

//  theFile = new TFile("hitsFile.root","RECREATE");



  hCSC000 = std::auto_ptr<TH2D>( new TH2D("CSC000","",500,-700,450,500,-600,550));
  hCSC004 = std::auto_ptr<TH2D>( new TH2D("CSC004","",500,-700,450,500,-600,550));
  hCSC005 = std::auto_ptr<TH2D>( new TH2D("CSC005","",500,-700,450,500,-600,550));

  hCSC010 = std::auto_ptr<TH2D>( new TH2D("CSC010","",500,-700,450,500,-600,550));
  hCSC014 = std::auto_ptr<TH2D>( new TH2D("CSC014","",500,-700,450,500,-600,550));
  hCSC015 = std::auto_ptr<TH2D>( new TH2D("CSC015","",500,-700,450,500,-600,550));

  hCSC020 = std::auto_ptr<TH2D>( new TH2D("CSC020","",500,-700,450,500,-600,550));
  hCSC024 = std::auto_ptr<TH2D>( new TH2D("CSC024","",500,-700,450,500,-600,550));
  hCSC025 = std::auto_ptr<TH2D>( new TH2D("CSC025","",500,-700,450,500,-600,550));

  hCSC030 = std::auto_ptr<TH2D>( new TH2D("CSC030","",500,-700,450,500,-600,550));
  hCSC034 = std::auto_ptr<TH2D>( new TH2D("CSC034","",500,-700,450,500,-600,550));
  hCSC035 = std::auto_ptr<TH2D>( new TH2D("CSC035","",500,-700,450,500,-600,550));

  hCSC040 = std::auto_ptr<TH2D>( new TH2D("CSC040","",500,-700,450,500,-600,550));
  hCSC044 = std::auto_ptr<TH2D>( new TH2D("CSC044","",500,-700,450,500,-600,550));
  hCSC045 = std::auto_ptr<TH2D>( new TH2D("CSC045","",500,-700,450,500,-600,550));


  hCSC001 = std::auto_ptr<TH2D>( new TH2D("CSC001","",500,-700,450,500,-600,550));
  hCSC002 = std::auto_ptr<TH2D>( new TH2D("CSC002","",500,-700,450,500,-600,550));
  hCSC003 = std::auto_ptr<TH2D>( new TH2D("CSC003","",500,-700,450,500,-600,550));

  hCSC011 = std::auto_ptr<TH2D>( new TH2D("CSC011","",500,-700,450,500,-600,550));
  hCSC012 = std::auto_ptr<TH2D>( new TH2D("CSC012","",500,-700,450,500,-600,550));
  hCSC013 = std::auto_ptr<TH2D>( new TH2D("CSC013","",500,-700,450,500,-600,550));

  hCSC021 = std::auto_ptr<TH2D>( new TH2D("CSC021","",500,-700,450,500,-600,550));
  hCSC022 = std::auto_ptr<TH2D>( new TH2D("CSC022","",500,-700,450,500,-600,550));
  hCSC023 = std::auto_ptr<TH2D>( new TH2D("CSC023","",500,-700,450,500,-600,550));

  hCSC031 = std::auto_ptr<TH2D>( new TH2D("CSC031","",500,-700,450,500,-600,550));
  hCSC032 = std::auto_ptr<TH2D>( new TH2D("CSC032","",500,-700,450,500,-600,550));
  hCSC033 = std::auto_ptr<TH2D>( new TH2D("CSC033","",500,-700,450,500,-600,550));

  hCSC041 = std::auto_ptr<TH2D>( new TH2D("CSC041","",500,-700,450,500,-600,550));
  hCSC042 = std::auto_ptr<TH2D>( new TH2D("CSC042","",500,-700,450,500,-600,550));
  hCSC043 = std::auto_ptr<TH2D>( new TH2D("CSC043","",500,-700,450,500,-600,550));

//-----------------------

  hCSC000->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC000","",500,-700,450,500,-600,550));
  hCSC004->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC004","",500,-700,450,500,-600,550));
  hCSC005->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC005","",500,-700,450,500,-600,550));

  hCSC010->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC010","",500,-700,450,500,-600,550));
  hCSC014->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC014","",500,-700,450,500,-600,550));
  hCSC015->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC015","",500,-700,450,500,-600,550));

  hCSC020->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC020","",500,-700,450,500,-600,550));
  hCSC024->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC024","",500,-700,450,500,-600,550));
  hCSC025->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC025","",500,-700,450,500,-600,550));

  hCSC030->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC030","",500,-700,450,500,-600,550));
  hCSC034->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC034","",500,-700,450,500,-600,550));
  hCSC035->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC035","",500,-700,450,500,-600,550));

  hCSC040->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC040","",500,-700,450,500,-600,550));
  hCSC044->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC044","",500,-700,450,500,-600,550));
  hCSC045->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC045","",500,-700,450,500,-600,550));


  hCSC001->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC001","",500,-700,450,500,-600,550));
  hCSC002->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC002","",500,-700,450,500,-600,550));
  hCSC003->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC003","",500,-700,450,500,-600,550));

  hCSC011->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC011","",500,-700,450,500,-600,550));
  hCSC012->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC012","",500,-700,450,500,-600,550));
  hCSC013->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC013","",500,-700,450,500,-600,550));

  hCSC021->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC021","",500,-700,450,500,-600,550));
  hCSC022->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC022","",500,-700,450,500,-600,550));
  hCSC023->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC023","",500,-700,450,500,-600,550));

  hCSC031->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC031","",500,-700,450,500,-600,550));
  hCSC032->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC032","",500,-700,450,500,-600,550));
  hCSC033->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC033","",500,-700,450,500,-600,550));

  hCSC041->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC041","",500,-700,450,500,-600,550));
  hCSC042->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC042","",500,-700,450,500,-600,550));
  hCSC043->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSC043","",500,-700,450,500,-600,550));

//============================================================================



  hCSCcsc000 = std::auto_ptr<TH1D>( new TH1D("CSCcsc000","",500,-100,100));
  hCSCcsc004 = std::auto_ptr<TH1D>( new TH1D("CSCcsc004","",500,-100,100));
  hCSCcsc005 = std::auto_ptr<TH1D>( new TH1D("CSCcsc005","",500,-100,100));

  hCSCcsc010 = std::auto_ptr<TH1D>( new TH1D("CSCcsc010","",500,-100,100));
  hCSCcsc014 = std::auto_ptr<TH1D>( new TH1D("CSCcsc014","",500,-100,100));
  hCSCcsc015 = std::auto_ptr<TH1D>( new TH1D("CSCcsc015","",500,-100,100));

  hCSCcsc020 = std::auto_ptr<TH1D>( new TH1D("CSCcsc020","",500,-100,100));
  hCSCcsc024 = std::auto_ptr<TH1D>( new TH1D("CSCcsc024","",500,-100,100));
  hCSCcsc025 = std::auto_ptr<TH1D>( new TH1D("CSCcsc025","",500,-100,100));

  hCSCcsc030 = std::auto_ptr<TH1D>( new TH1D("CSCcsc030","",500,-100,100));
  hCSCcsc034 = std::auto_ptr<TH1D>( new TH1D("CSCcsc034","",500,-100,100));
  hCSCcsc035 = std::auto_ptr<TH1D>( new TH1D("CSCcsc035","",500,-100,100));

  hCSCcsc040 = std::auto_ptr<TH1D>( new TH1D("CSCcsc040","",500,-100,100));
  hCSCcsc044 = std::auto_ptr<TH1D>( new TH1D("CSCcsc044","",500,-100,100));
  hCSCcsc045 = std::auto_ptr<TH1D>( new TH1D("CSCcsc045","",500,-100,100));


  hCSCcsc001 = std::auto_ptr<TH1D>( new TH1D("CSCcsc001","",500,-100,100));
  hCSCcsc002 = std::auto_ptr<TH1D>( new TH1D("CSCcsc002","",500,-100,100));
  hCSCcsc003 = std::auto_ptr<TH1D>( new TH1D("CSCcsc003","",500,-100,100));

  hCSCcsc011 = std::auto_ptr<TH1D>( new TH1D("CSCcsc011","",500,-100,100));
  hCSCcsc012 = std::auto_ptr<TH1D>( new TH1D("CSCcsc012","",500,-100,100));
  hCSCcsc013 = std::auto_ptr<TH1D>( new TH1D("CSCcsc013","",500,-100,100));

  hCSCcsc021 = std::auto_ptr<TH1D>( new TH1D("CSCcsc021","",500,-100,100));
  hCSCcsc022 = std::auto_ptr<TH1D>( new TH1D("CSCcsc022","",500,-100,100));
  hCSCcsc023 = std::auto_ptr<TH1D>( new TH1D("CSCcsc023","",500,-100,100));

  hCSCcsc031 = std::auto_ptr<TH1D>( new TH1D("CSCcsc031","",500,-100,100));
  hCSCcsc032 = std::auto_ptr<TH1D>( new TH1D("CSCcsc032","",500,-100,100));
  hCSCcsc033 = std::auto_ptr<TH1D>( new TH1D("CSCcsc033","",500,-100,100));

  hCSCcsc041 = std::auto_ptr<TH1D>( new TH1D("CSCcsc041","",500,-100,100));
  hCSCcsc042 = std::auto_ptr<TH1D>( new TH1D("CSCcsc042","",500,-100,100));
  hCSCcsc043 = std::auto_ptr<TH1D>( new TH1D("CSCcsc043","",500,-100,100));

//-----------------------

  hCSCcsc000->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc000","",500,-100,100));
  hCSCcsc004->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc004","",500,-100,100));
  hCSCcsc005->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc005","",500,-100,100));

  hCSCcsc010->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc010","",500,-100,100));
  hCSCcsc014->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc014","",500,-100,100));
  hCSCcsc015->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc015","",500,-100,100));

  hCSCcsc020->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc020","",500,-100,100));
  hCSCcsc024->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc024","",500,-100,100));
  hCSCcsc025->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc025","",500,-100,100));

  hCSCcsc030->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc030","",500,-100,100));
  hCSCcsc034->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc034","",500,-100,100));
  hCSCcsc035->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc035","",500,-100,100));

  hCSCcsc040->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc040","",500,-100,100));
  hCSCcsc044->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc044","",500,-100,100));
  hCSCcsc045->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc045","",500,-100,100));


  hCSCcsc001->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc001","",500,-100,100));
  hCSCcsc002->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc002","",500,-100,100));
  hCSCcsc003->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc003","",500,-100,100));

  hCSCcsc011->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc011","",500,-100,100));
  hCSCcsc012->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc012","",500,-100,100));
  hCSCcsc013->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc013","",500,-100,100));

  hCSCcsc021->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc021","",500,-100,100));
  hCSCcsc022->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc022","",500,-100,100));
  hCSCcsc023->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc023","",500,-100,100));

  hCSCcsc031->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc031","",500,-100,100));
  hCSCcsc032->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc032","",500,-100,100));
  hCSCcsc033->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc033","",500,-100,100));

  hCSCcsc041->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc041","",500,-100,100));
  hCSCcsc042->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc042","",500,-100,100));
  hCSCcsc043->SetDirectory(0); //::auto_ptr<TH2D>( new TH2D("CSCcsc043","",500,-100,100));



}

// ------------ method called once each job just after ending the event loop  ------------
void 
T1RecHitAnalyzer::endJob() {
  theFile = TFile::Open(outputFileName.c_str(), "recreate");
  if(!theFile || !theFile->IsWritable())
  {
    std::cout<<"Output file not opened correctly!!"<<std::endl;
  }
  writeHistograms();
//heFile->Write();
  theFile->Close();
}

void 
T1RecHitAnalyzer::writeHistograms(){
  hCSC000->Write();
  hCSC004->Write();
  hCSC005->Write();

  hCSC010->Write();
  hCSC014->Write();
  hCSC015->Write();

  hCSC020->Write();
  hCSC024->Write();
  hCSC025->Write();

  hCSC030->Write();
  hCSC034->Write();
  hCSC035->Write();

  hCSC040->Write();
  hCSC044->Write();
  hCSC045->Write();

  hCSC001->Write();
  hCSC002->Write();
  hCSC003->Write();

  hCSC011->Write();
  hCSC012->Write();
  hCSC013->Write();

  hCSC021->Write();
  hCSC022->Write();
  hCSC023->Write();

  hCSC031->Write();
  hCSC032->Write();
  hCSC033->Write();

  hCSC041->Write();
  hCSC042->Write();
  hCSC043->Write();
//===========================================================
  hCSCcsc000->Write();
  hCSCcsc004->Write();
  hCSCcsc005->Write();

  hCSCcsc010->Write();
  hCSCcsc014->Write();
  hCSCcsc015->Write();

  hCSCcsc020->Write();
  hCSCcsc024->Write();
  hCSCcsc025->Write();

  hCSCcsc030->Write();
  hCSCcsc034->Write();
  hCSCcsc035->Write();

  hCSCcsc040->Write();
  hCSCcsc044->Write();
  hCSCcsc045->Write();

  hCSCcsc001->Write();
  hCSCcsc002->Write();
  hCSCcsc003->Write();

  hCSCcsc011->Write();
  hCSCcsc012->Write();
  hCSCcsc013->Write();

  hCSCcsc021->Write();
  hCSCcsc022->Write();
  hCSCcsc023->Write();

  hCSCcsc031->Write();
  hCSCcsc032->Write();
  hCSCcsc033->Write();

  hCSCcsc041->Write();
  hCSCcsc042->Write();
  hCSCcsc043->Write();


}


float T1RecHitAnalyzer::Eta(float x,float y,float z){
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
float T1RecHitAnalyzer::Phi(float x,float y){
  float c=0;
  if(x>0 && y>0) c = atan(y/x);
  if(x<0) c = atan(y/x)+3.14159;
  if(x>0 && y<0) c = atan(y/x)+6.28318;
  return c;
}
