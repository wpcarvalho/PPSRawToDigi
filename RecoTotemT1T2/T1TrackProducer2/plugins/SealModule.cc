#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoTotemT1T2/T1TrackProducer2/interface/T1TrackProducer2.h"
#include "RecoTotemT1T2/T1TrackProducer2/interface/T1TrackProducer3.h"
#include "RecoTotemT1T2/T1TrackProducer2/interface/T1TrackProducer4.h"
#include "RecoTotemT1T2/T1TrackProducer2/interface/T1TrackAnalyzer.h"
#include "RecoTotemT1T2/T1TrackProducer2/interface/T1TrackAnalyzerTB.h"

DEFINE_FWK_MODULE(T1TrackProducer2);
DEFINE_FWK_MODULE(T1TrackProducer3);
DEFINE_FWK_MODULE(T1TrackProducer4);
DEFINE_FWK_MODULE(T1TrackAnalyzer);
DEFINE_FWK_MODULE(T1TrackAnalyzerTB);

