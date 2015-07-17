/**
 * Class T2DigiProducer.
 *
 * Author: Erik Br??cken / University of Helsinki
 *         Mirko Berretti / University of Siena & Pisa INFN
 * Email:  brucken@cc.helsinki.fi
 *         mirko.berretti@gmail.com  
 * Date:   2007-11-26
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimTotem/T2Digitizer/interface/T2DigiProducer.h"
#include "SimTotem/T2Digitizer/interface/T2Digitizer.h"
#include "Geometry/TotemGeometry/interface/T2Geometry.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"


#include <iostream>

T2DigiProducer::T2DigiProducer(const edm::ParameterSet& paraSet): DigiProdParSet_(paraSet)
{
  edm::Service<edm::RandomNumberGenerator> rngPR;
  if ( ! rngPR.isAvailable()) {
    throw cms::Exception("Configuration")
      << "This class requires the RandomNumberGeneratorService\n"
      "which is not present in the configuration file.  You must add the service\n"
      "in the configuration file or remove the modules that require it.";
  }
  
  rndEnginePR = &(rngPR->getEngine());

  saveDigiVFAT=paraSet.getParameter<bool>("saveDigiVFAT");

  previousModule = paraSet.getParameter<std::string>("previousModule");
  instanceLabel = paraSet.getParameter<std::string>("instanceLabel");
  simTrackContainerLabel = paraSet.getParameter<edm::InputTag>("SimTrackContainerLabel");
  SimVertexContainerLabel = paraSet.getParameter<edm::InputTag>("SimVertexContainerLabel");

  //previousModule = cms.string('mix'),
  // instanceLabel = cms.string('g4SimHitsTotemHitsT2Gem'),


  if (saveDigiVFAT == true) {
	  produces<T2DigiVfatCollection>("T2Digivfat");
  }

  produces<T2StripDigiCollection>("T2StripDigi");
  produces<T2PadDigiCollection>("T2PadDigi");
} // T2DigiProducer


void T2DigiProducer::produce(edm::Event& ev, const edm::EventSetup& evSet) {
  edm::Handle<CrossingFrame<PSimHit> > cFrame;
  ev.getByLabel(previousModule, instanceLabel, cFrame);

  // get hits from G4Sim
  const std::string nameOfHits("TotemHitsT2Gem");

  std::auto_ptr<MixCollection<PSimHit> > T2SimHits( new MixCollection<PSimHit>( cFrame.product() ) );

  // get tracks from G4Sim
  std::auto_ptr<edm::SimTrackContainer> theSimTracks (new edm::SimTrackContainer());
  edm::Handle<edm::SimTrackContainer> simTrackCollection;
  ev.getByLabel(simTrackContainerLabel, simTrackCollection);
  theSimTracks.get()->clear();
  for(edm::SimTrackContainer::const_iterator ittrack = simTrackCollection->begin();ittrack != simTrackCollection->end(); ++ittrack){ 
    SimTrack thetrack=(*ittrack);
    theSimTracks->push_back(thetrack);
  }

  //LOAD SIMULATED VERTICES
  //Container for all simvertices
  edm::Handle<edm::SimVertexContainer> simVertexCollection;
  ev.getByLabel(SimVertexContainerLabel, simVertexCollection);
  std::auto_ptr<edm::SimVertexContainer> theSimVertices (new edm::SimVertexContainer());
  theSimVertices.get()->clear();
  for(edm::SimVertexContainer::const_iterator itvertex = simVertexCollection->begin();itvertex != simVertexCollection->end(); ++itvertex){     SimVertex thevertex=(*itvertex);   
    theSimVertices->push_back(thevertex);
  }

  edm::Handle<edm::HepMCProduct> HepMCEvtHandle ;
  ev.getByLabel("generator", HepMCEvtHandle ) ;
  const HepMC::GenEvent* hepmcevt = HepMCEvtHandle->GetEvent();
  
 
  // create empty containers for digi
  std::auto_ptr<T2StripDigiCollection> theStripDigis(new T2StripDigiCollection());
  std::auto_ptr<T2PadDigiCollection> thePadDigis(new T2PadDigiCollection());
  


  //Declare Here the T2Digitizer pointer (should be optimized)
  std::auto_ptr<T2Digitizer> theT2Digitizerpt(new T2Digitizer(DigiProdParSet_,evSet,rndEnginePR));
 
  // process the hits
  (*theT2Digitizerpt).process(*T2SimHits, *theStripDigis, *thePadDigis, hepmcevt,theSimVertices,theSimTracks);


 
  //  --------------------------------------------------------------------------------------------------------------
  //Copy the vfat map in a vfatcollection
 
  //std::cout<<"Map size: "<<"Strip: "<<(*theT2Digitizerpt).StripVFatsdig->size()<<" Pad: "<<(*theT2Digitizerpt).PadVFatsdig->size()<<std::endl;
  if(saveDigiVFAT==true)
    {
       
      std::auto_ptr<T2DigiVfatCollection> theVfatsCopy(new T2DigiVfatCollection());
      std::map<unsigned int, T2DigiVfat>::iterator vfinplaneit;
      std::map<int, std::map<unsigned int, T2DigiVfat> >::iterator vfiter;  
      std::map<int, std::map<unsigned int, T2DigiVfat> >::iterator vfbeg=(*theT2Digitizerpt).StripVFatsdig->begin();
      std::map<int, std::map<unsigned int, T2DigiVfat> >::iterator vfend=(*theT2Digitizerpt).StripVFatsdig->end();
      T2GeometryUtil conv;
      T2GeometryUtil::T2DetInfo planeinformation;

      for (vfiter=vfbeg;vfiter!=vfend;vfiter++)
	{
	  planeinformation=conv.GetT2Info(vfiter->first);       
	  //std::cout<<"In digiproducer looks: arm-ht-pl-pls: "<<planeinformation.arm<<"-"<<planeinformation.ht<<"-"<<planeinformation.pl<<"-"<<planeinformation.plside<<std::endl;
	  T2DetId thet2det(planeinformation.cmsswid);
	  for (vfinplaneit=(vfiter->second).begin();vfinplaneit!=(vfiter->second).end();vfinplaneit++)
	    theVfatsCopy->insertDigi(thet2det,(*vfinplaneit).second);
	}
      //std::cout<<"strip saved "<<std::endl;
      
      std::map<int, std::map<unsigned int, T2DigiVfat> >::iterator vfiterp;  
      std::map<int, std::map<unsigned int, T2DigiVfat> >::iterator vfbegp=(*theT2Digitizerpt).PadVFatsdig->begin();
      std::map<int, std::map<unsigned int, T2DigiVfat> >::iterator vfendp=(*theT2Digitizerpt).PadVFatsdig->end();
      std::map<unsigned int, T2DigiVfat>::iterator vfinplaneitpad;
      for (vfiterp=vfbegp;vfiterp!=vfendp;vfiterp++)
	{
	  planeinformation=conv.GetT2Info(vfiterp->first);       
	  // std::cout<<"In digiproducer looks: arm-ht-pl-pls: "<<planeinformation.arm<<"-"<<planeinformation.ht<<"-"<<planeinformation.pl<<"-"<<planeinformation.plside<<std::endl;
	  T2DetId thet2det(planeinformation.cmsswid);
	  //std::cout<<"T2Det initizalized with "<<thet2det.arm()<<std::endl;
	  for (vfinplaneitpad=(vfiterp->second).begin();vfinplaneitpad!=(vfiterp->second).end();vfinplaneitpad++)
	    theVfatsCopy->insertDigi(thet2det,(*vfinplaneitpad).second);
	}
      //  --------------------------------------------------------------------------------------------------------------
      //std::cout<<"DONE!!"<<std::endl;
      ev.put(theVfatsCopy,"T2Digivfat");
    }
  

  // and put the result into the event
  ev.put(theStripDigis, "T2StripDigi");
  ev.put(thePadDigis, "T2PadDigi");
  
} // produce


void T2DigiProducer::beginJob(){

}



