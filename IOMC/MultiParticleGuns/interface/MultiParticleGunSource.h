#ifndef MultiParticleGunSource_H
#define MultiParticleGunSource_H

/** \class MultiParticleGuns from FlatRandomEGunSource
 *
 * Modified by Mirko Berretti (TOTEM exp)
 * 
 ***************************************/

#include "IOMC/MultiParticleGuns/interface/MultFlatGunSource.h"

namespace edm
{

  class MultiParticleGunSource : public MultFlatGunSource
  {
  
  public:
    MultiParticleGunSource(const ParameterSet&);
    virtual ~MultiParticleGunSource();

  private:
   
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
  protected :
  
    // data members
    std::vector<double> fMinEs;
    std::vector<double> fMaxEs;


  };
} 

#endif
