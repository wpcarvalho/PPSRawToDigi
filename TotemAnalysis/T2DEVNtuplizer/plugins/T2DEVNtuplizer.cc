/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Jan Ka??par (jan.kaspar@gmail.com) 
*    
* $Id: TotemRawMetaDataDumper.cc 3290 2010-09-07 11:54:29Z jkaspar $
* $Revision: 2145 $
* $Date: 2010-02-17 15:01:22 +0100 (Wed, 17 Feb 2010) $
*
****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "TotemAnalysis/T2DEVNtuplizer/interface/Ntuplizer.h"
#include "TotemAnalysis/T2DEVNtuplizer/interface/T2Ntuplizer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TTree.h"
#include "TFile.h"

#include <vector>
#include <string>

/**
T2DEVNtuplizer.cc
 * Common Ntuplizer for all TOTEM subdetectors.
 **/
class T2DEVNtuplizer : public edm::EDAnalyzer
{
  public:
    T2DEVNtuplizer(const edm::ParameterSet&);
    ~T2DEVNtuplizer() {}

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
//    RPNtuplizer *rp_ntupl_;
    
    // flag for branch creating
    bool branchCreated; 
};

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

T2DEVNtuplizer::T2DEVNtuplizer(const edm::ParameterSet &ps) :
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
  outputFileName(ps.getUntrackedParameter<string>("outputFileName"))
{
 
  workers.push_back(new T2Ntuplizer(ps));
   
  branchCreated = false;
}

//----------------------------------------------------------------------------------------------------

void T2DEVNtuplizer::beginRun(edm::Run const& r, edm::EventSetup const& es)
{
  /*
  edm::ESHandle<BeamOpticsParams> BOParH;
  es.get<BeamOpticsParamsRcd>().get(BOParH);
  if(!BOParH.isValid())
    throw cms::Exception("T2DEVNtuplizer::beginRun") << " edm::ESHandle<BeamOpticsParams> is invalid";
  
  rp_ntupl_->SetOpticsConfig(*BOParH);
  */

  // let all workes create their branches
  if( branchCreated == false){
    for (vector<Ntuplizer *>::iterator it = workers.begin(); it != workers.end(); ++it)
      (*it)->CreateBranches(es, tree);
    branchCreated = true;
  }
}

//----------------------------------------------------------------------------------------------------

void T2DEVNtuplizer::beginJob()
{
  // open dump file and preapre dump tree branches
  file = new TFile(outputFileName.c_str(), "recreate");
  tree = new TTree("TotemNtuple", "TotemNtuple");
}

//----------------------------------------------------------------------------------------------------

void T2DEVNtuplizer::analyze(const edm::Event &ev, const edm::EventSetup &es)
{
  // let all workes fill their event data
  for (vector<Ntuplizer *>::iterator it = workers.begin(); it != workers.end(); ++it) {
    try {
      (*it)->FillEvent(ev, es);
    }
    catch (...) {
    }
  }

  // commit the data
  tree->Fill();
}

//----------------------------------------------------------------------------------------------------

void T2DEVNtuplizer::endJob()
{
  file->cd();
  tree->Write();
  delete file;
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(T2DEVNtuplizer);

