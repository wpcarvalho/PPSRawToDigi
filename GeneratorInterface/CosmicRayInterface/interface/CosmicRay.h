#ifndef COSMICRAY_H
#define COSMICRAY_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "GeneratorInterface/Core/interface/BaseHadronizer.h"
#include "CLHEP/Random/RandomEngine.h"

#include <map>
#include <string>
#include <vector>
#include <math.h>

namespace HepMC {
  class GenEvent;
  class GenParticle;
  class GenVertex;
}


namespace gen
{
  class CosmicRayHadronizer : public BaseHadronizer {
  public:
    CosmicRayHadronizer(const edm::ParameterSet &);
    virtual ~CosmicRayHadronizer();  };

} /*end namespace*/

#endif //ifndef COSMICRAY_H
