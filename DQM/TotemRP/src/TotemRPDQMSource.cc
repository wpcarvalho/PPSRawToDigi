#include "DQM/TotemRP/interface/TotemRPDQMSource.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>

using namespace edm;

//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::TotemRPDQMSource(const edm::ParameterSet& ps)
  /*
  buildCorrelationPlots(ps.getUntrackedParameter<bool>("buildCorrelationPlots", false)),
  correlationPlotsFilter(ps.getUntrackedParameter<std::string>("correlationPlotsFilter", "")),
  correlationPlotsLimit(ps.getUntrackedParameter<unsigned int>("correlationPlotsLimit", 50)),
  correlationPlotsSelector(correlationPlotsFilter)
  */
{
  edm::LogInfo("TotemRPDQMSource") <<  "Constructor  TotemRPDQMSource::TotemRPDQMSource " << std::endl;
  
  tokenStripDigi = consumes< DetSetVector<RPStripDigi> >(ps.getParameter<edm::InputTag>("tagStripDigi"));
  tokenDigiCluster = consumes< edm::DetSetVector<RPDigCluster> >(ps.getParameter<edm::InputTag>("tagDigiCluster"));
  tokenRecoHit = consumes< edm::DetSetVector<RPRecoHit> >(ps.getParameter<edm::InputTag>("tagRecoHit"));
  tokenPatternColl = consumes< RPRecognizedPatternsCollection >(ps.getParameter<edm::InputTag>("tagPatternColl"));
  tokenTrackCandColl = consumes< RPTrackCandidateCollection >(ps.getParameter<edm::InputTag>("tagTrackCandColl"));
  tokenTrackColl = consumes< RPFittedTrackCollection >(ps.getParameter<edm::InputTag>("tagTrackColl"));
  tokenMultiTrackColl = consumes< RPMulFittedTrackCollection >(ps.getParameter<edm::InputTag>("tagMultiTrackColl"));
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

  //h_test = ibooker_.book1D("test", "some title", 40, -0.5, 39.5);
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, 
                                            edm::EventSetup const& context) 
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::beginLuminosityBlock" << std::endl;
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::analyze(edm::Event const& event, edm::EventSetup const& eSetup)
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::analyze" << std::endl;

  Handle< DetSetVector<RPStripDigi> > digi;
  event.getByToken(tokenStripDigi, digi);

  Handle< DetSetVector<RPDigCluster> > digCluster;
  event.getByToken(tokenDigiCluster, digCluster);

  Handle< DetSetVector<RPRecoHit> > hits;
  event.getByToken(tokenRecoHit, hits);

  Handle<RPRecognizedPatternsCollection> patterns;
  event.getByToken(tokenPatternColl, patterns);

  Handle< RPTrackCandidateCollection > trackCanColl;
  event.getByToken(tokenTrackCandColl, trackCanColl);

  Handle< RPFittedTrackCollection > tracks;
  event.getByToken(tokenTrackColl, tracks);

  Handle< RPMulFittedTrackCollection > multiTracks;
  event.getByToken(tokenMultiTrackColl, multiTracks);
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
