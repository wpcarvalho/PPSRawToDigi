// -*- C++ -*-
//
// Package:    T1Analysis
// Class:      T1Analysis
// 
/**\class T1Analysis T1Analysis.cc TotemAnalysis/T1Analysis/src/T1Analysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  fabrizio ferro
//         Created:  Mon Nov 28 09:33:34 CET 2011
// $Id$
//
//


// system include files
#include <memory>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"

#include <boost/shared_ptr.hpp>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

//
// class declaration
//

class T1Analysis : public edm::EDAnalyzer {
   public:
      explicit T1Analysis(const edm::ParameterSet&);
      ~T1Analysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------


  edm::InputTag rawEventLabel;
  unsigned int BXallT2Trigger;
  unsigned int BXfakeT2Trigger;
  unsigned int BXlostT2Trigger;
  unsigned int BXlostT1orT2Trigger;

  unsigned int BXatleast1T1Track;
  unsigned int BXatleast1T1orT2Track;
  unsigned int BXatleast1T2Track;
  unsigned int BXatleast1T1TrackNoT2Track;
  unsigned int BXatleast1T1TrackNoT2TrackNoT2TBIT;
  unsigned int BXatleast1T1TrackNoT2TrackYesT2TBIT;


  unsigned int BXBBallT2Trigger;
  unsigned int BXBBfakeT2Trigger;
  unsigned int BXBBlostT2Trigger;
  unsigned int BXBBlostT1orT2Trigger;

  unsigned int BXBBatleast1T1Track;
  unsigned int BXBBatleast1T1orT2Track;
  unsigned int BXBBatleast1T2Track;
  unsigned int BXBBatleast1T1TrackNoT2Track;
  unsigned int BXBBatleast1T1TrackNoT2TrackNoT2TBIT;
  unsigned int BXBBatleast1T1TrackNoT2TrackYesT2TBIT;

  unsigned int BXNC1allT2Trigger;
  unsigned int BXNC1fakeT2Trigger;
  unsigned int BXNC1lostT2Trigger;
  unsigned int BXNC1lostT1orT2Trigger;

  unsigned int BXNC1atleast1T1Track;
  unsigned int BXNC1atleast1T1orT2Track;
  unsigned int BXNC1atleast1T2Track;
  unsigned int BXNC1atleast1T1TrackNoT2Track;
  unsigned int BXNC1atleast1T1TrackNoT2TrackNoT2TBIT;
  unsigned int BXNC1atleast1T1TrackNoT2TrackYesT2TBIT;


  unsigned int BXNC2allT2Trigger;
  unsigned int BXNC2fakeT2Trigger;
  unsigned int BXNC2lostT2Trigger;
  unsigned int BXNC2lostT1orT2Trigger;

  unsigned int BXNC2atleast1T1Track;
  unsigned int BXNC2atleast1T1orT2Track;
  unsigned int BXNC2atleast1T2Track;
  unsigned int BXNC2atleast1T1TrackNoT2Track;
  unsigned int BXNC2atleast1T1TrackNoT2TrackNoT2TBIT;
  unsigned int BXNC2atleast1T1TrackNoT2TrackYesT2TBIT;


  unsigned int AllEvents;
  unsigned int AllEventsT2;
  unsigned int AllEventsBigBunch;
  unsigned int AllEventsNC1;
  unsigned int AllEventsNC2;

  std::auto_ptr<TH1D> hT1TracksBXfakeT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXallT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXlostT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXlostT1orT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXatleast1T1TrackNoT2Track;
  std::auto_ptr<TH1D> hT1TracksBXatleast1T1TrackNoT2TrackNoT2TBIT;
  std::auto_ptr<TH1D> hT1TracksBXatleast1T1TrackNoT2TrackYesT2TBIT;
  std::auto_ptr<TH1D> hT1TracksBX;

  std::auto_ptr<TH1D> hT1TracksBXNC1fakeT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXNC1allT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXNC1lostT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXNC1lostT1orT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXNC1atleast1T1TrackNoT2Track;
  std::auto_ptr<TH1D> hT1TracksBXNC1atleast1T1TrackNoT2TrackNoT2TBIT;
  std::auto_ptr<TH1D> hT1TracksBXNC1atleast1T1TrackNoT2TrackYesT2TBIT;
  std::auto_ptr<TH1D> hT1TracksBXNC1;

  std::auto_ptr<TH1D> hT1TracksBXNC2fakeT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXNC2allT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXNC2lostT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXNC2lostT1orT2Trigger;
  std::auto_ptr<TH1D> hT1TracksBXNC2atleast1T1TrackNoT2Track;
  std::auto_ptr<TH1D> hT1TracksBXNC2atleast1T1TrackNoT2TrackNoT2TBIT;
  std::auto_ptr<TH1D> hT1TracksBXNC2atleast1T1TrackNoT2TrackYesT2TBIT;
  std::auto_ptr<TH1D> hT1TracksBXNC2;

  TFile* theFile;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
T1Analysis::T1Analysis(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   rawEventLabel = iConfig.getParameter<edm::InputTag>("RawEventLabel");
}


T1Analysis::~T1Analysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
T1Analysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif



// Trigger data

  unsigned int TriggerData_bunch_num;
  unsigned int TriggerData_input_status_bits;

  Handle< Totem::RawEvent > input;
  iEvent.getByLabel(rawEventLabel, input);

  TriggerData_bunch_num = input->triggerData.bunch_num;
  TriggerData_input_status_bits = input->triggerData.input_status_bits;

// T1 data

    edm::Handle<T1T2TrackCollection> T1trackCollection;
    iEvent.getByLabel("t1tracks2","T1TrackColl",T1trackCollection);

    T1T2TrackCollection::const_iterator T1TC_it;


// T2 data

  Handle<T1T2TrackCollection> T2trackCollection;
  iEvent.getByLabel("T2TrackColl3","T2TrackColl",T2trackCollection);
 
  T1T2TrackCollection::const_iterator T2TC_it;

    unsigned int T2TBIT = TriggerData_input_status_bits & 64;
    unsigned int T2TBIT_HM = TriggerData_input_status_bits & 128;
    unsigned int BXBIT = TriggerData_input_status_bits & 512;

    if(T2TBIT > 0 || T2TBIT_HM > 0)AllEventsT2++;


    if(TriggerData_bunch_num == 0  ){ // only big bunch
    
    if(BXBIT > 0){   // analyse BX triggered events
      
  
      
      if(T2TBIT > 0 || T2TBIT_HM > 0){

	BXBBallT2Trigger++;
      }
      if((T2TBIT > 0 || T2TBIT_HM > 0) && T2trackCollection->size()==0){

	BXBBfakeT2Trigger++;
      }
      if(T2TBIT == 0 && T2TBIT_HM == 0 && T2trackCollection->size()>0){

	BXBBlostT2Trigger++;
      }
      if(T2TBIT == 0 && T2TBIT_HM == 0 && (T2trackCollection->size()>0 || T1trackCollection->size()>0)){

	BXBBlostT1orT2Trigger++;
      }
      
      if(T1trackCollection->size()>0){
	BXBBatleast1T1Track++;
      }
      if(T2trackCollection->size()>0){ 
	BXBBatleast1T2Track++;
      }
       if(T2trackCollection->size()>0 || T1trackCollection->size()>0){ 
	BXBBatleast1T1orT2Track++;
      }
      if(T2trackCollection->size()==0 && T1trackCollection->size()>0){ 

	BXBBatleast1T1TrackNoT2Track++;

      }

      if(T2trackCollection->size()==0 && T2TBIT == 0 && T1trackCollection->size()>0){ 

	BXBBatleast1T1TrackNoT2TrackNoT2TBIT++;
      }
      if(T2trackCollection->size()==0 && T2TBIT > 0 && T1trackCollection->size()>0){ 

	BXBBatleast1T1TrackNoT2TrackYesT2TBIT++;
      }

      AllEventsBigBunch++;
    }// analyse BX triggered events



  }


    
    if(BXBIT > 0){   // analyse BX triggered events
      
      hT1TracksBX->Fill(T1trackCollection->size());
      
      if(T2TBIT > 0 || T2TBIT_HM > 0){
	hT1TracksBXallT2Trigger->Fill(T1trackCollection->size());
	BXallT2Trigger++;
      }
      if((T2TBIT > 0 || T2TBIT_HM > 0) && T2trackCollection->size()==0){
	hT1TracksBXfakeT2Trigger->Fill(T1trackCollection->size());
	BXfakeT2Trigger++;
      }
      if(T2TBIT == 0 && T2TBIT_HM == 0 && T2trackCollection->size()>0){
	hT1TracksBXlostT2Trigger->Fill(T1trackCollection->size());
	BXlostT2Trigger++;
      }
      if(T2TBIT == 0 && T2TBIT_HM == 0 && (T2trackCollection->size()>0 || T1trackCollection->size()>0)){
	hT1TracksBXlostT1orT2Trigger->Fill(T1trackCollection->size());
	BXlostT1orT2Trigger++;
      }
      
      if(T1trackCollection->size()>0){
	BXatleast1T1Track++;
      }
      if(T2trackCollection->size()>0){ 
	BXatleast1T2Track++;
      }
       if(T2trackCollection->size()>0 || T1trackCollection->size()>0){ 
	BXatleast1T1orT2Track++;
      }
      if(T2trackCollection->size()==0 && T1trackCollection->size()>0){ 
	hT1TracksBXatleast1T1TrackNoT2Track->Fill(T1trackCollection->size());
	BXatleast1T1TrackNoT2Track++;

      }

      if(T2trackCollection->size()==0 && T2TBIT == 0 && T1trackCollection->size()>0){ 
	hT1TracksBXatleast1T1TrackNoT2TrackNoT2TBIT->Fill(T1trackCollection->size());
	BXatleast1T1TrackNoT2TrackNoT2TBIT++;
      }
      if(T2trackCollection->size()==0 && T2TBIT > 0 && T1trackCollection->size()>0){ 
	hT1TracksBXatleast1T1TrackNoT2TrackYesT2TBIT->Fill(T1trackCollection->size());
	BXatleast1T1TrackNoT2TrackYesT2TBIT++;
      }


    }// analyse BX triggered events





 if(TriggerData_bunch_num == 1785 ){ // only NON collidng bunches (fill 2232 runs 6872 - 6947)
    
    if(BXBIT > 0){   // analyse BX triggered events
      
      hT1TracksBXNC1->Fill(T1trackCollection->size());
      
      if(T2TBIT > 0 || T2TBIT_HM > 0){
	hT1TracksBXNC1allT2Trigger->Fill(T1trackCollection->size());
	BXNC1allT2Trigger++;
      }
      if((T2TBIT > 0 || T2TBIT_HM > 0) && T2trackCollection->size()==0){
	hT1TracksBXNC1fakeT2Trigger->Fill(T1trackCollection->size());
	BXNC1fakeT2Trigger++;
      }
      if(T2TBIT == 0 && T2TBIT_HM == 0 && T2trackCollection->size()>0){
	hT1TracksBXNC1lostT2Trigger->Fill(T1trackCollection->size());
	BXNC1lostT2Trigger++;
      }
            if(T2TBIT == 0 && T2TBIT_HM == 0 && (T2trackCollection->size()>0 || T1trackCollection->size()>0)){
	hT1TracksBXNC1lostT1orT2Trigger->Fill(T1trackCollection->size());
	BXNC1lostT1orT2Trigger++;
      }
      if(T1trackCollection->size()>0){
	BXNC1atleast1T1Track++;
      }
      if(T2trackCollection->size()>0){ 
	BXNC1atleast1T2Track++;
      }
       if(T2trackCollection->size()>0 || T1trackCollection->size()>0){ 
	BXNC1atleast1T1orT2Track++;
      }
      if(T2trackCollection->size()==0 && T1trackCollection->size()>0){ 
	hT1TracksBXNC1atleast1T1TrackNoT2Track->Fill(T1trackCollection->size());
	BXNC1atleast1T1TrackNoT2Track++;

      }

      if(T2trackCollection->size()==0 && T2TBIT == 0 && T1trackCollection->size()>0){ 
	hT1TracksBXNC1atleast1T1TrackNoT2TrackNoT2TBIT->Fill(T1trackCollection->size());
	BXNC1atleast1T1TrackNoT2TrackNoT2TBIT++;
      }
      if(T2trackCollection->size()==0 && T2TBIT > 0 && T1trackCollection->size()>0){ 
	hT1TracksBXNC1atleast1T1TrackNoT2TrackYesT2TBIT->Fill(T1trackCollection->size());
	BXNC1atleast1T1TrackNoT2TrackYesT2TBIT++;
      }


    }// analyse BXNC1 triggered events

    AllEventsNC1++;

  }


 if( TriggerData_bunch_num == 891 ){ // only NON collidng bunches (fill 2232 runs 6872 - 6947)
    
    if(BXBIT > 0){   // analyse BX triggered events
      
      hT1TracksBXNC2->Fill(T1trackCollection->size());
      
      if(T2TBIT > 0 || T2TBIT_HM > 0){
	hT1TracksBXNC2allT2Trigger->Fill(T1trackCollection->size());
	BXNC2allT2Trigger++;
      }
      if((T2TBIT > 0 || T2TBIT_HM > 0) && T2trackCollection->size()==0){
	hT1TracksBXNC2fakeT2Trigger->Fill(T1trackCollection->size());
	BXNC2fakeT2Trigger++;
      }
      if(T2TBIT == 0 && T2TBIT_HM == 0 && T2trackCollection->size()>0){
	hT1TracksBXNC2lostT2Trigger->Fill(T1trackCollection->size());
	BXNC2lostT2Trigger++;
      }
            if(T2TBIT == 0 && T2TBIT_HM == 0 && (T2trackCollection->size()>0 || T1trackCollection->size()>0)){
	hT1TracksBXNC2lostT1orT2Trigger->Fill(T1trackCollection->size());
	BXNC2lostT1orT2Trigger++;
      }
      if(T1trackCollection->size()>0){
	BXNC2atleast1T1Track++;
      }
      if(T2trackCollection->size()>0){ 
	BXNC2atleast1T2Track++;
      }
       if(T2trackCollection->size()>0 || T1trackCollection->size()>0){ 
	BXNC2atleast1T1orT2Track++;
      }
      if(T2trackCollection->size()==0 && T1trackCollection->size()>0){ 
	hT1TracksBXNC2atleast1T1TrackNoT2Track->Fill(T1trackCollection->size());
	BXNC2atleast1T1TrackNoT2Track++;

      }

      if(T2trackCollection->size()==0 && T2TBIT == 0 && T1trackCollection->size()>0){ 
	hT1TracksBXNC2atleast1T1TrackNoT2TrackNoT2TBIT->Fill(T1trackCollection->size());
	BXNC2atleast1T1TrackNoT2TrackNoT2TBIT++;
      }
      if(T2trackCollection->size()==0 && T2TBIT > 0 && T1trackCollection->size()>0){ 
	hT1TracksBXNC2atleast1T1TrackNoT2TrackYesT2TBIT->Fill(T1trackCollection->size());
	BXNC2atleast1T1TrackNoT2TrackYesT2TBIT++;
      }


    }// analyse BXNC2 triggered events

    AllEventsNC2++;

  }


  AllEvents++;

}


// ------------ method called once each job just before starting event loop  ------------
void 
T1Analysis::beginJob()
{
   BXallT2Trigger=0;
   BXfakeT2Trigger=0;
   BXlostT2Trigger=0;
   BXlostT1orT2Trigger=0;

   BXatleast1T1Track=0;
   BXatleast1T1orT2Track=0;
   BXatleast1T2Track=0;
   BXatleast1T1TrackNoT2Track=0;
   BXatleast1T1TrackNoT2TrackNoT2TBIT=0;
   BXatleast1T1TrackNoT2TrackYesT2TBIT=0;
  
   BXNC1allT2Trigger=0;
   BXNC1fakeT2Trigger=0;
   BXNC1lostT2Trigger=0;
   BXNC1lostT1orT2Trigger=0;

   BXNC1atleast1T1Track=0;
   BXNC1atleast1T1orT2Track=0;
   BXNC1atleast1T2Track=0;
   BXNC1atleast1T1TrackNoT2Track=0;
   BXNC1atleast1T1TrackNoT2TrackNoT2TBIT=0;
   BXNC1atleast1T1TrackNoT2TrackYesT2TBIT=0;
  

   BXBBallT2Trigger=0;
   BXBBfakeT2Trigger=0;
   BXBBlostT2Trigger=0;
   BXBBlostT1orT2Trigger=0;

   BXBBatleast1T1Track=0;
   BXBBatleast1T1orT2Track=0;
   BXBBatleast1T2Track=0;
   BXBBatleast1T1TrackNoT2Track=0;
   BXBBatleast1T1TrackNoT2TrackNoT2TBIT=0;
   BXBBatleast1T1TrackNoT2TrackYesT2TBIT=0;
  

   BXNC2allT2Trigger=0;
   BXNC2fakeT2Trigger=0;
   BXNC2lostT2Trigger=0;
   BXNC2lostT1orT2Trigger=0;

   BXNC2atleast1T1Track=0;
   BXNC2atleast1T1orT2Track=0;
   BXNC2atleast1T2Track=0;
   BXNC2atleast1T1TrackNoT2Track=0;
   BXNC2atleast1T1TrackNoT2TrackNoT2TBIT=0;
   BXNC2atleast1T1TrackNoT2TrackYesT2TBIT=0;
  

   AllEvents=0;
   AllEventsT2=0;
   AllEventsBigBunch=0;
   AllEventsNC1=0;
   AllEventsNC2=0;

 
 
  hT1TracksBX = std::auto_ptr<TH1D>(new TH1D("hT1TracksBX","hT1TracksBX",100,-0.5,99.5));
  hT1TracksBXfakeT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXfakeT2Trigger","hT1TracksBXfakeT2Trigger",100,-0.5,99.5));
  hT1TracksBXlostT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXlostT2Trigger","hT1TracksBXlostT2Trigger",100,-0.5,99.5));
  hT1TracksBXlostT1orT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXlostT1orT2Trigger","hT1TracksBXlostT1orT2Trigger",100,-0.5,99.5));
  hT1TracksBXallT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXallT2Trigger","hT1TracksBXallT2Trigger",100,-0.5,99.5));
  hT1TracksBXatleast1T1TrackNoT2Track = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXatleast1T1TrackNoT2Track","hT1TracksBXatleast1T1TrackNoT2Track",100,-0.5,99.5));
  hT1TracksBXatleast1T1TrackNoT2TrackNoT2TBIT = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXatleast1T1TrackNoT2TrackNoT2TBIT","hT1TracksBXatleast1T1TrackNoT2TrackNoT2TBIT",100,-0.5,99.5));
  hT1TracksBXatleast1T1TrackNoT2TrackYesT2TBIT = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXatleast1T1TrackNoT2TrackYesT2TBIT","hT1TracksBXatleast1T1TrackNoT2TrackYesT2TBIT",100,-0.5,99.5));

  hT1TracksBX->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBX","hT1TracksBX",100,-0.5,99.5));
  hT1TracksBXfakeT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXfakeT2Trigger","hT1TracksBXfakeT2Trigger",100,-0.5,99.5));
  hT1TracksBXlostT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXlostT2Trigger","hT1TracksBXlostT2Trigger",100,-0.5,99.5));
  hT1TracksBXlostT1orT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXlostT2Trigger","hT1TracksBXlostT2Trigger",100,-0.5,99.5));
  hT1TracksBXallT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXallT2Trigger","hT1TracksBXallT2Trigger",100,-0.5,99.5));
  hT1TracksBXatleast1T1TrackNoT2Track->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXatleast1T1TrackNoT2Track","hT1TracksBXatleast1T1TrackNoT2Track",100,-0.5,99.5));
  hT1TracksBXatleast1T1TrackNoT2TrackNoT2TBIT->SetDirectory(0);
  hT1TracksBXatleast1T1TrackNoT2TrackYesT2TBIT->SetDirectory(0);

  hT1TracksBXNC1 = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1","hT1TracksBXNC1",100,-0.5,99.5));
  hT1TracksBXNC1fakeT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1fakeT2Trigger","hT1TracksBXNC1fakeT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1lostT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1lostT2Trigger","hT1TracksBXNC1lostT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1lostT1orT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1lostT1orT2Trigger","hT1TracksBXNC1lostT1orT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1allT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1allT2Trigger","hT1TracksBXNC1allT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1atleast1T1TrackNoT2Track = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1atleast1T1TrackNoT2Track","hT1TracksBXNC1atleast1T1TrackNoT2Track",100,-0.5,99.5));
  hT1TracksBXNC1atleast1T1TrackNoT2TrackNoT2TBIT = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1atleast1T1TrackNoT2TrackNoT2TBIT","hT1TracksBXNC1atleast1T1TrackNoT2TrackNoT2TBIT",100,-0.5,99.5));
  hT1TracksBXNC1atleast1T1TrackNoT2TrackYesT2TBIT = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1atleast1T1TrackNoT2TrackYesT2TBIT","hT1TracksBXNC1atleast1T1TrackNoT2TrackYesT2TBIT",100,-0.5,99.5));

  hT1TracksBXNC1->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1","hT1TracksBXNC1",100,-0.5,99.5));
  hT1TracksBXNC1fakeT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1fakeT2Trigger","hT1TracksBXNC1fakeT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1lostT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1lostT2Trigger","hT1TracksBXNC1lostT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1lostT1orT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1lostT2Trigger","hT1TracksBXNC1lostT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1allT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1allT2Trigger","hT1TracksBXNC1allT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1atleast1T1TrackNoT2Track->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1atleast1T1TrackNoT2Track","hT1TracksBXNC1atleast1T1TrackNoT2Track",100,-0.5,99.5));
  hT1TracksBXNC1atleast1T1TrackNoT2TrackNoT2TBIT->SetDirectory(0);
  hT1TracksBXNC1atleast1T1TrackNoT2TrackYesT2TBIT->SetDirectory(0);

  hT1TracksBXNC2 = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2","hT1TracksBXNC2",100,-0.5,99.5));
  hT1TracksBXNC2fakeT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2fakeT2Trigger","hT1TracksBXNC2fakeT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2lostT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2lostT2Trigger","hT1TracksBXNC2lostT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2lostT1orT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2lostT1orT2Trigger","hT1TracksBXNC2lostT1orT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2allT2Trigger = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2allT2Trigger","hT1TracksBXNC2allT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2atleast1T1TrackNoT2Track = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2atleast1T1TrackNoT2Track","hT1TracksBXNC2atleast1T1TrackNoT2Track",100,-0.5,99.5));
  hT1TracksBXNC2atleast1T1TrackNoT2TrackNoT2TBIT = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2atleast1T1TrackNoT2TrackNoT2TBIT","hT1TracksBXNC2atleast1T1TrackNoT2TrackNoT2TBIT",100,-0.5,99.5));
  hT1TracksBXNC2atleast1T1TrackNoT2TrackYesT2TBIT = std::auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2atleast1T1TrackNoT2TrackYesT2TBIT","hT1TracksBXNC2atleast1T1TrackNoT2TrackYesT2TBIT",100,-0.5,99.5));

  hT1TracksBXNC2->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2","hT1TracksBXNC2",100,-0.5,99.5));
  hT1TracksBXNC2fakeT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2fakeT2Trigger","hT1TracksBXNC2fakeT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2lostT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2lostT2Trigger","hT1TracksBXNC2lostT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2lostT1orT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2lostT2Trigger","hT1TracksBXNC2lostT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2allT2Trigger->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2allT2Trigger","hT1TracksBXNC2allT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2atleast1T1TrackNoT2Track->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2atleast1T1TrackNoT2Track","hT1TracksBXNC2atleast1T1TrackNoT2Track",100,-0.5,99.5));
  hT1TracksBXNC2atleast1T1TrackNoT2TrackNoT2TBIT->SetDirectory(0);
  hT1TracksBXNC2atleast1T1TrackNoT2TrackYesT2TBIT->SetDirectory(0);


}
// ------------ method called once each job just after ending the event loop  ------------
void 
T1Analysis::endJob() 
{

  theFile = TFile::Open("T1Analysis.root", "recreate");
  if(!theFile || !theFile->IsWritable())
    {
      std::cout<<"Output file not opened correctly!!"<<std::endl;
    }

   hT1TracksBX->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBX","hT1TracksBX",100,-0.5,99.5));
  hT1TracksBXfakeT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXfakeT2Trigger","hT1TracksBXfakeT2Trigger",100,-0.5,99.5));
  hT1TracksBXlostT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXlostT2Trigger","hT1TracksBXlostT2Trigger",100,-0.5,99.5));
  hT1TracksBXlostT1orT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXlostT2Trigger","hT1TracksBXlostT2Trigger",100,-0.5,99.5));
  hT1TracksBXallT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXallT2Trigger","hT1TracksBXallT2Trigger",100,-0.5,99.5));
  hT1TracksBXatleast1T1TrackNoT2Track->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXatleast1T1TrackNoT2Track","hT1TracksBXatleast1T1TrackNoT2Track",100,-0.5,99.5));
  hT1TracksBXatleast1T1TrackNoT2TrackNoT2TBIT->Write(); 
  hT1TracksBXatleast1T1TrackNoT2TrackYesT2TBIT->Write(); 





   hT1TracksBXNC1->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1","hT1TracksBXNC1",100,-0.5,99.5));
  hT1TracksBXNC1fakeT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1fakeT2Trigger","hT1TracksBXNC1fakeT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1lostT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1lostT2Trigger","hT1TracksBXNC1lostT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1lostT1orT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1lostT2Trigger","hT1TracksBXNC1lostT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1allT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1allT2Trigger","hT1TracksBXNC1allT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC1atleast1T1TrackNoT2Track->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC1atleast1T1TrackNoT2Track","hT1TracksBXNC1atleast1T1TrackNoT2Track",100,-0.5,99.5));
  hT1TracksBXNC1atleast1T1TrackNoT2TrackNoT2TBIT->Write(); 
  hT1TracksBXNC1atleast1T1TrackNoT2TrackYesT2TBIT->Write(); 









   hT1TracksBXNC2->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2","hT1TracksBXNC2",100,-0.5,99.5));
  hT1TracksBXNC2fakeT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2fakeT2Trigger","hT1TracksBXNC2fakeT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2lostT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2lostT2Trigger","hT1TracksBXNC2lostT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2lostT1orT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2lostT2Trigger","hT1TracksBXNC2lostT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2allT2Trigger->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2allT2Trigger","hT1TracksBXNC2allT2Trigger",100,-0.5,99.5));
  hT1TracksBXNC2atleast1T1TrackNoT2Track->Write(); //auto_ptr<TH1D>(new TH1D("hT1TracksBXNC2atleast1T1TrackNoT2Track","hT1TracksBXNC2atleast1T1TrackNoT2Track",100,-0.5,99.5));
  hT1TracksBXNC2atleast1T1TrackNoT2TrackNoT2TBIT->Write(); 
  hT1TracksBXNC2atleast1T1TrackNoT2TrackYesT2TBIT->Write(); 












  theFile->Close();






  cout << "All read events = " << AllEvents << endl;
  cout << "All events with T2 trigger = " << AllEventsT2 << endl;
  cout << "All events from big bunch = " << AllEventsBigBunch << endl;
  cout << "All events from NC1 bunches= " << AllEventsNC1 << endl;
  cout << "All events from NC2 bunches= " << AllEventsNC2 << endl;
  cout << "BXallT2Trigger="<<BXallT2Trigger << "   BXfakeT2Trigger="<<  BXfakeT2Trigger << "  BXlostT2Trigger="<<BXlostT2Trigger << "  BXlostT1orT2Trigger="<<BXlostT1orT2Trigger<<endl;
  cout << "BXatleast1T1Track="<<BXatleast1T1Track<<"   BXatleast1T2Track="<<BXatleast1T2Track<<"   BXatleast1T1orT2Track="<<BXatleast1T1orT2Track<<"   BXatleast1T1TrackNoT2Track="<<BXatleast1T1TrackNoT2Track<<"   BXatleast1T1TrackNoT2TrackNoT2TBIT="<<BXatleast1T1TrackNoT2TrackNoT2TBIT<<"   BXatleast1T1TrackNoT2TrackYesT2TBIT="<<BXatleast1T1TrackNoT2TrackYesT2TBIT<<endl;
  cout << "BXBBallT2Trigger="<<BXBBallT2Trigger << "   BXBBfakeT2Trigger="<<  BXBBfakeT2Trigger << "  BXBBlostT2Trigger="<<BXBBlostT2Trigger<< "  BXBBlostT1orT2Trigger="<<BXBBlostT1orT2Trigger<<endl;
  cout << "BXBBatleast1T1Track="<<BXBBatleast1T1Track<<"   BXBBatleast1T2Track="<<BXBBatleast1T2Track<<"   BXBBatleast1T1orT2Track="<<BXBBatleast1T1orT2Track<<"   BXBBatleast1T1TrackNoT2Track="<<BXBBatleast1T1TrackNoT2Track<<"   BXBBatleast1T1TrackNoT2TrackNoT2TBIT="<<BXBBatleast1T1TrackNoT2TrackNoT2TBIT<<"   BXBBatleast1T1TrackNoT2TrackYesT2TBIT="<<BXBBatleast1T1TrackNoT2TrackYesT2TBIT<<endl;
  cout << "BXNC1allT2Trigger="<<BXNC1allT2Trigger << "   BXNC1fakeT2Trigger="<<  BXNC1fakeT2Trigger << "  BXNC1lostT2Trigger="<<BXNC1lostT2Trigger<< "  BXNC1lostT1orT2Trigger="<<BXNC1lostT1orT2Trigger<<endl;
  cout << "BXNC1atleast1T1Track="<<BXNC1atleast1T1Track<<"   BXNC1atleast1T2Track="<<BXNC1atleast1T2Track<<"   BXNC1atleast1T1orT2Track="<<BXNC1atleast1T1orT2Track<<"   BXNC1atleast1T1TrackNoT2Track="<<BXNC1atleast1T1TrackNoT2Track<<"   BXNC1atleast1T1TrackNoT2TrackNoT2TBIT="<<BXNC1atleast1T1TrackNoT2TrackNoT2TBIT<<"   BXNC1atleast1T1TrackNoT2TrackYesT2TBIT="<<BXNC1atleast1T1TrackNoT2TrackYesT2TBIT<<endl;

  cout << "BXNC2allT2Trigger="<<BXNC2allT2Trigger << "   BXNC2fakeT2Trigger="<<  BXNC2fakeT2Trigger << "  BXNC2lostT2Trigger="<<BXNC2lostT2Trigger<< "  BXNC2lostT1orT2Trigger="<<BXNC2lostT1orT2Trigger<<endl;
  cout << "BXNC2atleast1T1Track="<<BXNC2atleast1T1Track<<"   BXNC2atleast1T2Track="<<BXNC2atleast1T2Track<<"   BXNC2atleast1T1orT2Track="<<BXNC2atleast1T1orT2Track<<"   BXNC2atleast1T1TrackNoT2Track="<<BXNC2atleast1T1TrackNoT2Track<<"   BXNC2atleast1T1TrackNoT2TrackNoT2TBIT="<<BXNC2atleast1T1TrackNoT2TrackNoT2TBIT<<"   BXNC2atleast1T1TrackNoT2TrackYesT2TBIT="<<BXNC2atleast1T1TrackNoT2TrackYesT2TBIT<<endl;

}

// ------------ method called when starting to processes a run  ------------
void 
T1Analysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
T1Analysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
T1Analysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
T1Analysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
T1Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(T1Analysis);
