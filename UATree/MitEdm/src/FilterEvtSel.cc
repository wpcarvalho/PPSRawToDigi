//--------------------------------------------------------------------------------------------------
// $Id: FilterEvtSel.cc,v 1.10 2010/06/23 15:01:16 bendavid Exp $
//
// FilterEvtSel
//
// Filter to select events based on event selection object,
// containing calo energy/time, pixel hits, and cluster shape.
//
// Authors: E.Wenger
//--------------------------------------------------------------------------------------------------

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "UATree/MitEdm/interface/EvtSelData.h"
#include <TMath.h>

namespace mitedm
{
  class FilterEvtSel : public edm::EDFilter {
  public:
    explicit FilterEvtSel(const edm::ParameterSet &ps);
    ~FilterEvtSel() {}
    
  protected:
    virtual bool filter (edm::Event &iEvent, const edm::EventSetup &iSetup);

    double              minHfEnergy_;   //minimum hf energy
    double              maxHfTimeDiff_; //maximum hf energy
    std::string         srcEvtSel_;     //event selection data string
    std::vector<double> clusterPars_;   //pixel cluster polynomial pars for vertex compatibility cut
    int                 nhitsTrunc_;    //maximum pixel clusters to apply compatibility check
    int                 nhitsmax_;      //maximum number of pixel clusters
    int                 nHfHits_;       //minimum number of hf coincidence hits
    int                 nHfTowers_;     //minimum number of hf coincidence hits
    double              clusterTrunc_;  //maximum vertex compatibility value for event rejection
  };
}

using namespace mitedm;
using namespace edm;
using namespace std;

//--------------------------------------------------------------------------------------------------
FilterEvtSel::FilterEvtSel(const edm::ParameterSet& iConfig)
  : minHfEnergy_(iConfig.getUntrackedParameter<double>("minHfEnergy",0)),
    maxHfTimeDiff_(iConfig.getUntrackedParameter<double>("maxHfTimeDiff",0)),
    srcEvtSel_(iConfig.getUntrackedParameter<std::string>("srcEvtSel","evtSelData")),
    clusterPars_(iConfig.getUntrackedParameter< std::vector<double> >("clusterPars")),
    nhitsTrunc_(iConfig.getUntrackedParameter<int>("nhitsTrunc",0)),
    nhitsmax_(iConfig.getUntrackedParameter<int>("nhitsmax",0)),
    nHfHits_(iConfig.getUntrackedParameter<int>("nHfHits",0)),
    nHfTowers_(iConfig.getUntrackedParameter<int>("nHfTowers",0)),
    clusterTrunc_(iConfig.getUntrackedParameter<double>("clusterTrunc",0))
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
bool FilterEvtSel::filter( edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  // Filter events based on calo energy/timing and pixel clusters.
  
  Handle<EvtSelData> evtSel;
  iEvent.getByLabel(edm::InputTag(srcEvtSel_),evtSel);

  int nPxlHits       = evtSel->ePxHits();
  double hfTimeDiff  = evtSel->eHfPosTime() - evtSel->eHfNegTime();
  double clusVtxQual = evtSel->eClusVtxQual();
  double hfEnergyMin = min(evtSel->eHfPos(),evtSel->eHfNeg());

  // construct polynomial cut on cluster vertex quality vs. npixelhits
  double polyCut=0;
  for(unsigned int i=0; i < clusterPars_.size(); i++) {
    polyCut += clusterPars_[i]*pow((double)nPxlHits,(int)i);
  }
  if(nPxlHits < nhitsTrunc_) 
    polyCut=0;             // don't use cut below nhitsTrunc_ pixel hits
  if(polyCut > clusterTrunc_ && clusterTrunc_ > 0) 
    polyCut=clusterTrunc_; // no cut above clusterTrunc_

  bool accepted = true;

  if ( (TMath::Abs(hfTimeDiff)>maxHfTimeDiff_ && maxHfTimeDiff_>0) || 
       hfEnergyMin < minHfEnergy_                                  || 
       clusVtxQual < polyCut                                       ||
       (nPxlHits > nhitsmax_ && nhitsmax_>0)                       ||
       (evtSel->nHfNegHits() < nHfHits_)                           ||
       (evtSel->nHfPosHits() < nHfHits_)                           ||
       (evtSel->nHfTowersP() < nHfTowers_)                         ||
       (evtSel->nHfTowersN() < nHfTowers_) )
    accepted = false;

  return accepted;
}

//define this as a plug-in
DEFINE_FWK_MODULE(FilterEvtSel);
