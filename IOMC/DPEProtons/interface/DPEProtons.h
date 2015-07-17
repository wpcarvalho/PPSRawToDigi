#ifndef IOMC_DPEProtons_DPEProtons_h
#define IOMC_DPEProtons_DPEProtons_h

#include <string>
#include <vector>
#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"

#include "HepMC/GenEvent.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include <memory>
#include "boost/shared_ptr.hpp"

#include "IOMC/DPEProtons/interface/h101.h"

using namespace edm;
namespace HepMC {
	class GenEvent;
}


class DPEProtons : public edm::EDProducer
{
  public:
    DPEProtons(const edm::ParameterSet &);
    virtual ~DPEProtons();
  
  private:
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void beginRun(Run&, EventSetup const&);
  
  protected:
    // data members
    std::vector<std::string> file_names;
    h101 *dpe_file_access;
    long max_entries;
    
    HepMC::GenEvent* fEvt_;
    edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable;
    int verbosity_;
    CLHEP::HepRandomEngine* fRandomEngine;
    CLHEP::RandFlat* fRandomGenerator;
    
    double m0_;
    int PartID_;
};


#endif
