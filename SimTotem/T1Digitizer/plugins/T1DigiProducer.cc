/*
 * modified by Marcin Borratynski (mborratynski@gmail.com)
 */
#include "SimTotem/T1Digitizer/interface/T1DigiProducer.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/PCrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

#include "CLHEP/Matrix/SymMatrix.h"
#include <iostream>
#include <fstream>
#include <iomanip>

//#define _DEBUG_

T1DigiProducer::T1DigiProducer(const edm::ParameterSet& params) :
  _electronics(params.getParameter<std::string>("Electronics"))
{
  theVerbosity = params.getParameter<int>("Verbosity");

  if(theVerbosity >=1 ){std::cout << " Inside T1DigiProducer "<<std::endl;

  std::cout << "Electronics  --- " << _electronics << " ---"<<std::endl;
  }
  _digitizer = new T1Digitizer(params);

  if(_electronics == "Santiard"){
    produces<T1DigiSantiardCollection>("T1DigiSantiard");
  }
  else if(_electronics == "VFAT"){
    produces<T1DigiVfatCollection>("T1DigiVfat");
  }
  else{
    throw cms::Exception("")<<"Electronics for cathodes not set";
  }
  produces<T1DigiWireCollection>("T1DigiWire");

  simulateDeadChannels = false;
  if (params.exists("simulateDeadChannels")) { //check if "simulateDeadChannels" variable is defined in configuration file
	simulateDeadChannels = params.getParameter<bool> ("simulateDeadChannels");
  }

}


T1DigiProducer::~T1DigiProducer() 
{
  delete _digitizer; 
}


// ------------ method called once each job just before starting event loop  ------------
void T1DigiProducer::beginRun(edm::Run&, edm::EventSetup const& es) {
	// get analysis mask to mask channels
	if (simulateDeadChannels) {
		edm::ESHandle<AnalysisMask> analysisMask;
		es.get<TotemDAQMappingRecord> ().get(analysisMask);
		_deadChannelManager = T1DeadChannelManager(analysisMask); //set analysisMask in deadChannelsManager
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void T1DigiProducer::endJob() {
}


void T1DigiProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup) 
{
  evento=event.id().event();
  if(theVerbosity >=2)
  std::cout << " T1DigiProducer -- Event " << evento <<std::endl;
  /*
    int i_myColl=0;
    std::vector<edm::Handle<std::vector<PSimHit> > > resultsim;
    event.getManyByType(resultsim);
    int ss=resultsim.size();
    std::cout << "TAGLIA " << ss << std::endl;
    for (int ii=0;ii<ss;ii++) {
    edm::BranchDescription desc = resultsim[ii].provenance()->product();
    edm::LogInfo("T1DP") <<"For "<<desc.productInstanceName_<<" " <<resultsim[ii].product()->size()<<" Simhits added";
    if(desc.productInstanceName_=="TotemHitsT1"){
    i_myColl=ii;
    std::cout << "DA T1 " << ss << std::endl;
    }     
    }

  */

  edm::Handle<CrossingFrame<PSimHit> > cFrame;
  event.getByLabel("mix", "g4SimHitsTotemHitsT1", cFrame);

  // get hits from G4Sim
  const std::string nameOfHits("TotemHitsT1");

  std::auto_ptr<MixCollection<PSimHit> > 
    T1SimHits( new MixCollection<PSimHit>( cFrame.product() ) );

  std::auto_ptr<T1DigiWireCollection> t1Result(new T1DigiWireCollection);

  if(_electronics == "Santiard"){
    std::auto_ptr<T1DigiSantiardCollection> t1Result2 (new T1DigiSantiardCollection);
    
    _digitizer->doAction(evento,*T1SimHits, *t1Result, *t1Result2, _deadChannelManager);
    event.put(t1Result2,"T1DigiSantiard");

  }
  else if(_electronics == "VFAT"){
    std::auto_ptr<T1DigiVfatCollection> t1Result2 (new T1DigiVfatCollection);
    // run the algorithm
    _digitizer->doAction(evento,*T1SimHits, *t1Result, *t1Result2, _deadChannelManager);
    event.put(t1Result2,"T1DigiVfat");
  }
  event.put(t1Result,"T1DigiWire");

}

