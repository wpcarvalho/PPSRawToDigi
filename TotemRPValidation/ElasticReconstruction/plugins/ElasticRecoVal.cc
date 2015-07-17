/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors: 
*	Leszek Grzanka
*    
* $Id: ElasticRecoVal.cc 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/


#include <stdio.h>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TotemRPValidation/ElasticReconstruction/interface/ElasticRecoVal.h"

//----------------------------------------------------------------------------------------------------

ElasticRecoVal::ElasticRecoVal(const edm::ParameterSet& conf)
{
  libVal = std::auto_ptr<ElasticRecoValLibrary>(new ElasticRecoValLibrary(conf));
  libAcc = std::auto_ptr<ElasticAcceptanceLibrary>(new ElasticAcceptanceLibrary(conf));
}

//----------------------------------------------------------------------------------------------------

ElasticRecoVal::~ElasticRecoVal()
{
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoVal::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  libVal->initialize(es);
  libAcc->initialize(es);
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoVal::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  libVal->analyze(event, eSetup);
  libAcc->analyze(event, eSetup);
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoVal::endJob()
{
  libVal->finalize();
  libVal->writeHistogramsToFile();
  libAcc->finalize();
  libAcc->writeHistogramsToFile();
}

DEFINE_FWK_MODULE(ElasticRecoVal);

