/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "TotemAnalysis/TotemNtuplizer/interface/Ntuplizer.h"
#include "TotemAnalysis/TotemNtuplizer/interface/RawMetaDataNtuplizer.h"
#include "TotemAnalysis/TotemNtuplizer/interface/TriggerDataNtuplizer.h"
#include "TotemAnalysis/TotemNtuplizer/interface/RPNtuplizer.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
#include "DataFormats/TotemDigi/interface/TotemTriggerCounters.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPMulFittedTrackCollection.h"

#include "TTree.h"
#include "TFile.h"

#include <vector>
#include <string>

/**
 * Common Ntuplizer for all TOTEM subdetectors.
 **/
class TotemNtuplizer : public edm::EDAnalyzer
{
  public:
    TotemNtuplizer(const edm::ParameterSet&);
    ~TotemNtuplizer() {}

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    
    virtual void beginJob();
    virtual void endJob();

  protected:
    unsigned int verbosity;

    /// the name of the ROOT file with the final ntuple
    std::string outputFileName;

    /// the file where to save the ntuple
    TFile *file;   

    /// the ntuple
    TTree *tree;

    /// list of ntuple-making objects
    std::vector<Ntuplizer *> workers;
    
    /// internal reference to RPNtuplizer
    RPNtuplizer *rp_ntupl_;
    
    /// whether the tree branches have been declared
    bool branchesCreated; 
};

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

TotemNtuplizer::TotemNtuplizer(const edm::ParameterSet &ps) :
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
  outputFileName(ps.getUntrackedParameter<string>("outputFileName")),
  rp_ntupl_(NULL),
  branchesCreated(false)
{
  // instantiate selected worker modules
  vector<string> modules = ps.getParameter< vector<string> >("modules");

  for (const auto &m : modules)
  {
    if (m.compare("raw") == 0)
      workers.push_back(new RawMetaDataNtuplizer(ps));
      
    if (m.compare("trigger") == 0)
      workers.push_back(new TriggerDataNtuplizer(ps));

    if (m.compare("rp") == 0)
      workers.push_back(rp_ntupl_ = new RPNtuplizer(ps));
  }

  // declare what is consumed
  for (auto &w : workers)
    w->DeclareConsumes(this);

  // TODO: move int DeclareConsumes methods
  auto triggerCountersLabel(ps.getParameter<edm::InputTag>("TriggerCountersLabel"));
  consumes<TotemTriggerCounters>(triggerCountersLabel);
  
  auto rpDigClusterLabel = ps.getParameter<edm::InputTag>("RPDigClusterLabel");
  consumes<DetSetVector<TotemRPCluster>>(rpDigClusterLabel);

  auto rpUVPatternLabel = ps.getParameter<edm::InputTag>("RPUVPatternLabel");
  consumes<DetSetVector<TotemRPUVPattern>>(rpUVPatternLabel);

  auto tagLocalTrack = ps.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
  consumes<DetSetVector<TotemRPLocalTrack>>(tagLocalTrack);

  auto rpMulFittedTrackCollectionLabel = ps.getParameter<edm::InputTag>("RPMulFittedTrackCollectionLabel");
  consumes<RPMulFittedTrackCollection>(rpMulFittedTrackCollectionLabel);
}

//----------------------------------------------------------------------------------------------------

void TotemNtuplizer::beginRun(edm::Run const& r, edm::EventSetup const& es)
{
  // TODO
  /*
  edm::ESHandle<BeamOpticsParams> BOParH;
  es.get<BeamOpticsParamsRcd>().get(BOParH);
  if(!BOParH.isValid())
    throw cms::Exception("TotemNtuplizer::beginRun") << " edm::ESHandle<BeamOpticsParams> is invalid";
  
  rp_ntupl_->SetOpticsConfig(*BOParH);
  */

  // let all workes create their branches
  if (branchesCreated == false)
  {
    for (auto &w : workers)
      w->CreateBranches(es, tree);

    branchesCreated = true;
  }
}

//----------------------------------------------------------------------------------------------------

void TotemNtuplizer::beginJob()
{
  // open dump file and preapre dump tree branches
  file = new TFile(outputFileName.c_str(), "recreate");
  tree = new TTree("TotemNtuple", "TotemNtuple");
}

//----------------------------------------------------------------------------------------------------

void TotemNtuplizer::analyze(const edm::Event &ev, const edm::EventSetup &es)
{
  // let all workes fill their event data
  for (auto &w : workers)
    w->FillEvent(ev, es);

  // commit the data
  tree->Fill();
}

//----------------------------------------------------------------------------------------------------

void TotemNtuplizer::endJob()
{
  file->cd();
  tree->Write();
  delete file;
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemNtuplizer);
