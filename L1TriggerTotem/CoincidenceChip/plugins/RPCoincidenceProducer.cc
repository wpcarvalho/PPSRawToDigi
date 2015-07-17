#include "L1TriggerTotem/CoincidenceChip/interface/RPCoincidenceProducer.h"

//
// constructors and destructor
//
RPCoincidenceProducer::RPCoincidenceProducer(const edm::ParameterSet& iConfig) :
  verbose_(iConfig.getParameter<unsigned int> ("verbose")) {
  
  // configure the coincidence chip
  cchip.configure(iConfig.getParameterSet("coincidenceChipConfig"));
  //cchip.setV(3);

  productLabelSimu = iConfig.getParameter<std::string>("productLabelSimu");
  detTriggerLabel =  iConfig.getParameter<edm::InputTag>("DetTriggerLabel");
  produces <std::vector<RPCCBits> > (productLabelSimu);

}

RPCoincidenceProducer::~RPCoincidenceProducer() {
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void RPCoincidenceProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;

  Handle<DetSetVector<RPDetTrigger> > input;
  iEvent.getByLabel(detTriggerLabel, input );

  DetSetVector<RPDetTrigger>::const_iterator inputIterator = input->begin();

  if( verbose_ > 2) LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": number of RPDetTrigger objects:  " << input->size();

  std::map< RPDetId, std::bitset<80> > uBitsMap;
  std::map< RPDetId, std::bitset<80> > vBitsMap;
  int triggeredSectorNo, detectorNo;
  std::bitset<16> triggerBits;
  triggerBits.reset();

  for (; inputIterator != input->end(); inputIterator++) {
    TotRPDetId detectorId(inputIterator->id);

    // inputIterator->data : vector< RPDetTrigger>
    // (inputIterator->data)[i] : RPDetTrigger
    for (unsigned int i = 0; i < (inputIterator->data).size(); ++i) {
      triggeredSectorNo = (unsigned int) ((inputIterator->data)[i].GetSector());
      triggerBits.set(triggeredSectorNo);
      if( verbose_ > 2) LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": trigger bit in " << detectorId << " Sector: " << triggeredSectorNo ;
    }

    // create RAWID for detector 0
    // detector = 1,3,5..., odd, U
    detectorNo = (int) (detectorId.Detector());
    if (detectorId.IsStripsCoordinateUDirection()) {
      // (det-1)/2
      TotRPDetId firstDetectorId(detectorId.Arm(), detectorId.Station(), detectorId.RomanPot(), 1);
      RPDetId rawFirstDetectorId = firstDetectorId.rawId();
      for (unsigned int i = 0; i < 16; ++i) {
        uBitsMap[rawFirstDetectorId][i + 8 * (detectorNo - 1)] = triggerBits[i];
      }
    } else {
      // detector = 0,2,4..., even, V
      // det/2
      TotRPDetId firstDetectorId(detectorId.Arm(), detectorId.Station(), detectorId.RomanPot(), 0);
      RPDetId rawFirstDetectorId = firstDetectorId.rawId();
      for (unsigned int i = 0; i < 16; ++i) {
        vBitsMap[rawFirstDetectorId][i + 8 * detectorNo] = triggerBits[i];
      }
    }
    triggerBits.reset();
  }

  std::auto_ptr<std::vector<RPCCBits> > output(new std::vector<RPCCBits>);

  const int size = uBitsMap.size() + vBitsMap.size();
  output->reserve(size);

  // We print vBits maps
  for (map<RPDetId, std::bitset<80> >::const_iterator vIterator = vBitsMap.begin(); vIterator != vBitsMap.end(); ++vIterator) {
    if( verbose_ ) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": V (even: 0,2,..8) : " << TotRPDetId(vIterator->first);
    if( verbose_ > 2) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": V  input bits: " << vIterator->second;
    if( verbose_) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": input bits";
    if( verbose_ ) printBits(vIterator->second);
    if( verbose_ ) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": V output bits (simulated coincidence): " << cchip.process(vIterator->second);
    const RPCCIdRaw rawChipId = (RPCCIdRaw) (vIterator->first);
    const RPCCBits coinicidenceBits(rawChipId, cchip.process(vIterator->second));
    output->push_back(coinicidenceBits);
  }

  // We print uBits maps
  for (
      map<RPDetId, std::bitset<80> >::const_iterator uIterator = uBitsMap.begin();
      uIterator != uBitsMap.end();
      ++uIterator) {
    if( verbose_ ) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": U (odd: 1,3,..9) : " << TotRPDetId(uIterator->first);
    if( verbose_ > 2) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": U  input bits: " << uIterator->second;
    if( verbose_) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": input bits";
    if( verbose_ ) printBits(uIterator->second);
    if( verbose_ ) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": U output bits (simulated coincidence): " << cchip.process(uIterator->second);
    const RPCCIdRaw rawChipId = (RPCCIdRaw) (uIterator->first);
    const RPCCBits coincidenceBits(rawChipId, cchip.process(uIterator->second));
    output->push_back(coincidenceBits);
  }

  if( verbose_ ) LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer " << iEvent.id() << ": number of RPCCBits objects " << output->size();
  iEvent.put(output, productLabelSimu);
}

// ------------ method called once each job just before starting event loop  ------------
void RPCoincidenceProducer::beginJob() {
  if( verbose_ ) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer: beginJob";
  if( verbose_ ) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer: cfg summary : " << cchip.configSummary();
}

// ------------ method called once each job just after ending the event loop  ------------
void RPCoincidenceProducer::endJob() {
  if( verbose_ ) edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer: endJob";
}

void RPCoincidenceProducer::printBits( std::bitset<80> bs ){
   std::bitset<16> s0,s1,s2,s3,s4;
   for( unsigned int i = 0 ; i < 16 ; i++ ){
		s0[i] = bs[i];
		s1[i] = bs[i+16];
		s2[i] = bs[i+2*16];
		s3[i] = bs[i+3*16];
		s4[i] = bs[i+4*16];
   }
   edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer: " << s0;
   edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer: " << s1;
   edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer: " << s2;
   edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer: " << s3;
   edm::LogPrint("RPCoincidenceProducer") << "RPCoincidenceProducer: " << s4;
}

//define this as a plug-in
DEFINE_FWK_MODULE(RPCoincidenceProducer);
