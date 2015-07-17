#ifndef MultFlatGunSourceEtaEMC2_H
#define MultFlatGunSourceEtaEMC2_H

/** \class MultFlatGunSourceEtaEMC2
 *
 * Generates single particle gun in HepMC format
 * Julia Yarba 10/2005 
 * Modified By Mirko Berretti
 ***************************************/
#include <string>

#include "HepPDT/defs.h"
// #include "HepPDT/DefaultConfig.hh"
#include "HepPDT/TableBuilder.hh"
//#include "HepPDT/ParticleDataTableT.hh"
#include "HepPDT/ParticleDataTable.hh"

#include "HepMC/GenEvent.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include <memory>
#include "boost/shared_ptr.hpp"
#include <TH2.h>
#include <TFile.h>


//#include <boost/shared_ptr.hpp>
namespace edm
{
  
  class MultFlatGunSourceEtaEMC2 : public EDProducer
  {
  
  public:
    MultFlatGunSourceEtaEMC2(const ParameterSet&);
    // BaseFlatGunSource( const ParameterSet& ) ;
    virtual ~MultFlatGunSourceEtaEMC2();
    virtual void beginRun(Run&, EventSetup const&);
  private:
   
  protected :
  
    // non-virtuals ! this and only way !
    //
    // data members
    
    // gun particle(s) characteristics
    std::vector<int> fPartIDs ;
    std::vector<double>           fMinEtas ;
    std::vector<double>           fMaxEtas ;
    std::vector<double>           fMinPhis ;
    std::vector<double>           fMaxPhis ;

    // the event format itself
    HepMC::GenEvent* fEvt;

    // HepMC/HepPDT related things 
    // (for particle/event construction)
    //std::string      fPDGTablePath ;
    //std::string      fPDGTableName ; 
    // DefaultConfig::ParticleDataTable* fPDGTable;
    // DefaultConfig::ParticleDataTable* fTestTable ;
    // ESHandle<DefaultConfig::ParticleDataTable> fPDGTable ;
    ESHandle<HepPDT::ParticleDataTable> fPDGTable ;
            	    	
    int              fVerbosity ;

    CLHEP::HepRandomEngine* fRandomEngine ;
    CLHEP::RandFlat*        fRandomGenerator;
    
    bool             fAddAntiParticle;
    
    boost::shared_ptr<TFile> MCfile;
    //TFile *MCfile;
    std::string InputMCFileName;
    std::string InputMCHistoName;

    TH2F* MCHisto;



  };
} 

#endif
