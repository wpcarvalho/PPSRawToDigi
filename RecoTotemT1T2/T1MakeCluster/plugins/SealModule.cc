
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoTotemT1T2/T1MakeCluster/interface/T1MakeCluster.h"
#include "RecoTotemT1T2/T1MakeCluster/interface/T1MakeClusterMix.h"

DEFINE_FWK_MODULE(T1MakeCluster);
DEFINE_FWK_MODULE(T1MakeClusterMix);
//DEFINE_SEAL_PLUGIN (RPCRecHitAlgoFactory, RPCRecHitStandardAlgo, "RPCRecHitStandardAlgo");

