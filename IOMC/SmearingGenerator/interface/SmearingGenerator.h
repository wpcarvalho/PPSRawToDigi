/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors: 
 *  Jan Kašpar (jan.kaspar@gmail.com) 
 *    
 * $$RCSfile: SmearingGenerator.h,v $: $
 * $Revision: 9977 $
 * $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
 *
 ****************************************************************************/

#ifndef SmearingGenerator_H
#define SmearingGenerator_H

#include "FWCore/Framework/interface/EDProducer.h"

namespace CLHEP {
  class RandGauss;
}


/**
\brief Class for angular, energy and vertex smearing.

All internal variables are in GeV or in mm.

It is an EDProducer, but it doesn't produce anything reasonable.
 **/

class SmearingGenerator : public edm::EDProducer
{
  public:  
    SmearingGenerator(const edm::ParameterSet&);
    virtual ~SmearingGenerator();
    
  private:
    virtual void produce(edm::Event&, const edm::EventSetup&);
//    virtual void beginRun(edm::Run&, edm::EventSetup const&);
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    
  protected:
    unsigned int verbosity;                         ///< verbosity level
    std::string modifyLabel;                        ///< label of module which provided the HepMC to be modified
    std::string originalLabel;                      ///< instance label of the product to which the original HepMC shall be saved
    
    double MeanXi, SiXi;                            ///< Xi distribution
    double Al;                                      ///< Alpha (half crossing angle) d.
    double SiThX, SiThY;                            ///< Theta (beam divergence) d.
    
    double E_CMS;                                   ///< CMS energy of a proton
    
    double MeanX, MeanY, MeanZ;                     ///< [mm] means
    double SiX, SiY, SiZ;                           ///< [mm] variance
    
    void ApplyVertexSmearing(HepMC::GenEvent *evt);
    void ApplyBeamSmearing(HepMC::GenEvent *evt);
    
    CLHEP::RandGauss *rand;                         ///< gaussian random number generator
};

#endif

