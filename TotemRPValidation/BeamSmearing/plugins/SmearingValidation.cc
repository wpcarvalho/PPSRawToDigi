/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: SmearingValidation.cc 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#include <stdio.h>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TotemRPValidation/BeamSmearing/interface/SmearingValLibrary.h"

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

/**
 *\brief An analyzer performing validation of SmearingGenerator via SmearingValLibrary.
**/
class SmearingValidation : public edm::EDAnalyzer
{
 public:
  explicit SmearingValidation(const edm::ParameterSet&);
  ~SmearingValidation();

 private:

   std::auto_ptr<SmearingValLibrary> library;

  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

};


//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

SmearingValidation::SmearingValidation(const edm::ParameterSet& conf)
{
  library = std::auto_ptr<SmearingValLibrary>(new SmearingValLibrary(conf));
}

//----------------------------------------------------------------------------------------------------

SmearingValidation::~SmearingValidation()
{
}

//----------------------------------------------------------------------------------------------------

void SmearingValidation::beginJob()
{
  library->initialize();
}

//----------------------------------------------------------------------------------------------------

void SmearingValidation::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  library->analyze(event, eSetup);
}

//----------------------------------------------------------------------------------------------------

void SmearingValidation::endJob()
{
  library->finalize();
  library->writeHistogramsToFile();
}

DEFINE_FWK_MODULE(SmearingValidation);

