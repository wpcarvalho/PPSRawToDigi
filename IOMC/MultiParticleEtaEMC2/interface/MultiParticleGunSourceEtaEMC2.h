#ifndef MultiParticleGunSourceEtaEMC2_H
#define MultiParticleGunSourceEtaEMC2_H

/** \class MultiParticleGuns from FlatRandomEGunSource
 *
 * 
 * 
 ***************************************/

#include "IOMC/MultiParticleEtaEMC2/interface/MultFlatGunSourceEtaEMC2.h"

namespace edm
{

  class MultiParticleGunSourceEtaEMC2 : public MultFlatGunSourceEtaEMC2
  {
  
  public:
    MultiParticleGunSourceEtaEMC2(const ParameterSet &);
    virtual ~MultiParticleGunSourceEtaEMC2();

  private:
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
  protected :
  
    // data members
    std::vector<double> fMinEs;
    std::vector<double> fMaxEs;


  };
} 

#endif
