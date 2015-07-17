#ifndef COSMICRAYGENERATORFILTER_H
#define COSMICRAYGENERATORFILTER_H

#include "GeneratorInterface/CosmicRayInterface/interface/CosmicRayHadronizer.h"
#include "GeneratorInterface/Core/interface/GeneratorFilter.h"
#include "GeneratorInterface/ExternalDecays/interface/ExternalDecayDriver.h"

namespace gen
{
   typedef edm::GeneratorFilter<gen::CosmicRayHadronizer, gen::ExternalDecayDriver> CosmicRayGeneratorFilter;
}

#endif //#ifndef COSMICRAYGENERATORFILTER_H
