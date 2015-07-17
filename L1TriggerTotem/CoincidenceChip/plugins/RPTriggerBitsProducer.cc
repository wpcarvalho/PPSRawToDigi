#include "L1TriggerTotem/CoincidenceChip/interface/RPTriggerBitsProducer.h"

//
// constructors and destructor
//
RPTriggerBitsProducer::RPTriggerBitsProducer(const edm::ParameterSet& iConfig) :
  verbose_(iConfig.getParameter<bool> ("verbose")) {
	stripDigiLabel = iConfig.getParameter<edm::InputTag> ("StripDigiLabel");
  //now do what ever other initialization is needed
  produces<edm::DetSetVector<RPDetTrigger> >();
}

RPTriggerBitsProducer::~RPTriggerBitsProducer() {
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void RPTriggerBitsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  
  if( verbose_ ) LogPrint("RPTriggerBitsProducer") << "RPTriggerBitsProducer " << iEvent.id() <<  ": produce event id = " << iEvent.id();

  Handle<DetSetVector<RPStripDigi> > input;
  iEvent.getByLabel(stripDigiLabel, input );

  DetSetVector<RPStripDigi>::const_iterator inputIterator = input->begin();

  if( verbose_ ) LogPrint("RPTriggerBitsProducer") << "RPTriggerBitsProducer " << iEvent.id() <<  ": total number RPStripDigi objects: " << input->size();

   // Step B: LOOP on hits in event
   theTriggerVector.reserve(240);
   theTriggerVector.clear();


  int triggeredSectorNo, stripNo;
  std::bitset<512> triggerBits;
  triggerBits.reset();
  const int stripsPerSector = 32;

  for (; inputIterator != input->end(); inputIterator++) {
    TotRPDetId detectorId(inputIterator->id);

    if( verbose_ ) LogPrint("RPTriggerBitsProducer") << "RPTriggerBitsProducer " << iEvent.id() <<  ": " << detectorId;
    if( verbose_ ) LogPrint("RPTriggerBitsProducer") << "RPTriggerBitsProducer " << iEvent.id() <<  ": number of triggered strips: " << (inputIterator->data).size();

	edm::DetSet<RPDetTrigger> trigger_collector(inputIterator->id);

    vector<int> trig_cont;

    // inputIterator->data : vector< RPStripDigi>
    // (inputIterator->data)[i] : RPStripDigi
    for (unsigned int i = 0; i < (inputIterator->data).size(); ++i) {
      stripNo = (unsigned int) ((inputIterator->data)[i].GetStripNo());
      triggeredSectorNo = stripNo/stripsPerSector;
  	  if( verbose_ ) LogPrint("RPTriggerBitsProducer") << "RPTriggerBitsProducer " << iEvent.id() <<  ": strip no : " << stripNo;
  	  if( verbose_ ) LogPrint("RPTriggerBitsProducer") << "RPTriggerBitsProducer " << iEvent.id() <<  ": sector no : " << triggeredSectorNo;

	 // add if not present
  	 vector<int>::iterator p = find(trig_cont.begin(), trig_cont.end(), triggeredSectorNo);
	 if (p == trig_cont.end()){
	 	RPDetTrigger dt( inputIterator->id , triggeredSectorNo );   	 
     	trigger_collector.data.push_back(dt);
     	trig_cont.push_back(triggeredSectorNo);
	 }

    }

    if(trigger_collector.data.size()>0)
      {
        theTriggerVector.push_back(trigger_collector);
      }

   }
  
   // Step C: create empty output collection
    std::auto_ptr<edm::DetSetVector<RPDetTrigger> > trigger_output(new edm::DetSetVector<RPDetTrigger>(theTriggerVector) );

    if( verbose_ ) LogPrint("RPTriggerBitsProducer") << "RPTriggerBitsProducer " << iEvent.id() <<  ": number of RPDetTrigger objects: " << trigger_output->size();
    
    // Step D: write output to file
    iEvent.put(trigger_output);
  
}

// ------------ method called once each job just before starting event loop  ------------
void RPTriggerBitsProducer::beginJob() {
	if( verbose_ ) edm::LogPrint("RPTriggerBitsProducer") << "RPTriggerBitsProducer: beginJob";
}

// ------------ method called once each job just after ending the event loop  ------------
void RPTriggerBitsProducer::endJob() {
	if( verbose_ ) edm::LogPrint("RPTriggerBitsProducer") << "RPTriggerBitsProducer: endJob";
}

//define this as a plug-in
DEFINE_FWK_MODULE(RPTriggerBitsProducer);
