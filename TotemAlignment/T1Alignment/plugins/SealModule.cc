#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "TotemAlignment/T1Alignment/interface/T1InternalAlignment.h"
#include "TotemAlignment/T1Alignment/interface/T1InternalAlignment2.h"
#include "TotemAlignment/T1Alignment/interface/T1Sex.h"
#include "TotemAlignment/T1Alignment/interface/T1Sex2.h"


DEFINE_FWK_MODULE(T1InternalAlignment);
DEFINE_FWK_MODULE(T1InternalAlignment2);
DEFINE_FWK_MODULE(T1Sex);
DEFINE_FWK_MODULE(T1Sex2);



