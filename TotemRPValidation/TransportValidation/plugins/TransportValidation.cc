/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: TransportValidation.cc,v $: $
* $Revision: 1.2 $
* $Date: 2009/02/04 16:23:22 $
*
****************************************************************************/

#include "TotemRPValidation/TransportValidation/interface/TransportValidationLibrary.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

using namespace std;


/**
 *\brief Analyzer for a validation of proton transport functionality (via TransportValidationLibrary)
**/
class TransportValidation : public edm::EDAnalyzer
{
 public:
  TransportValidation(const edm::ParameterSet &ps) : tvl(new TransportValidationLibrary(ps)) {}
  ~TransportValidation() {}

 private:

  auto_ptr<TransportValidationLibrary> tvl;

  virtual void beginRun(edm::Run const&, edm::EventSetup const& es) { tvl->initialize(es); }
  virtual void analyze(const edm::Event &e, const edm::EventSetup &es) {tvl->analyze(e, es); }
  virtual void endJob() { tvl->finalize(); tvl->writeHistogramsToFile(); }
};

DEFINE_FWK_MODULE(TransportValidation);

