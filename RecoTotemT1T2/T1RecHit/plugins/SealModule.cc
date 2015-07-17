#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoTotemT1T2/T1RecHit/interface/T1RecHitProducer.h"
#include "RecoTotemT1T2/T1RecHit/interface/T1RecHitAnalyzer.h"

DEFINE_FWK_MODULE(T1RecHitProducer);
DEFINE_FWK_MODULE(T1RecHitAnalyzer);
//DEFINE_SEAL_PLUGIN (RPCRecHitAlgoFactory, RPCRecHitStandardAlgo, "RPCRecHitStandardAlgo");

