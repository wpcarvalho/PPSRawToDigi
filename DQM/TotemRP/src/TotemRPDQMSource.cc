#include "DQM/TotemRP/interface/TotemRPDQMSource.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>

//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::TotemRPDQMSource(const edm::ParameterSet& ps)
{
  edm::LogInfo("TotemRPDQMSource") <<  "Constructor  TotemRPDQMSource::TotemRPDQMSource " << std::endl;
  
  // Get parameters from configuration file
  //theElectronCollection_   = consumes<reco::GsfElectronCollection>(ps.getParameter<edm::InputTag>("electronCollection"));
  //triggerFilter_           = ps.getParameter<edm::InputTag>("TriggerFilter");
}

//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::~TotemRPDQMSource()
{
  edm::LogInfo("TotemRPDQMSource") <<  "Destructor TotemRPDQMSource::~TotemRPDQMSource " << std::endl;
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::dqmBeginRun(edm::Run const &, edm::EventSetup const &)
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::beginRun" << std::endl;
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::bookHistograms(DQMStore::IBooker & ibooker_, edm::Run const &, edm::EventSetup const &)
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::bookHistograms" << std::endl;
  
  ibooker_.cd();
  ibooker_.setCurrentFolder("TotemRP");

  h_test = ibooker_.book1D("test", "some title", 40, -0.5, 39.5);
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, 
                                            edm::EventSetup const& context) 
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::beginLuminosityBlock" << std::endl;
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::analyze(edm::Event const& e, edm::EventSetup const& eSetup)
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::analyze" << std::endl;

  h_test->Fill(10);
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& eSetup) 
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::endLuminosityBlock" << std::endl;
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::endRun(edm::Run const& run, edm::EventSetup const& eSetup)
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::endRun" << std::endl;
}
