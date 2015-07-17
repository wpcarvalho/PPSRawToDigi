#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoTotemT1T2/T1RoadProducer/interface/T1RoadProducer.h"


#include "RecoTotemT1T2/T1RoadProducer/interface/T1RoadAnalyzer.h"
#include "RecoTotemT1T2/T1RoadProducer/interface/T1RoadProducerWOaplane.h"

DEFINE_FWK_MODULE(T1RoadProducer);

DEFINE_FWK_MODULE(T1RoadProducerWOaplane);
DEFINE_FWK_MODULE(T1RoadAnalyzer);
