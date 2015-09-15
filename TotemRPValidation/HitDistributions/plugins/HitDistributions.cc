/*
 * HitDistributions.cc
 *
 *  Created on: Oct 29, 2008
 *      Author: grzanka
 */

#include <stdio.h>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TotemRPValidation/HitDistributions/interface/HitDistributions.h"

//----------------------------------------------------------------------------------------------------

HitDistributions::HitDistributions(const edm::ParameterSet& conf)
{
  library = std::auto_ptr<HitDistributionsLibrary>(new HitDistributionsLibrary(conf));
}

//----------------------------------------------------------------------------------------------------

HitDistributions::~HitDistributions()
{
}

//----------------------------------------------------------------------------------------------------

void HitDistributions::beginJob()
{
}

//----------------------------------------------------------------------------------------------------

void HitDistributions::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  library->analyze(event,eSetup);
}

//----------------------------------------------------------------------------------------------------

void HitDistributions::endJob()
{
  library->ExportAllHistograms();
  library->writeToFile();
}

DEFINE_FWK_MODULE(HitDistributions);
