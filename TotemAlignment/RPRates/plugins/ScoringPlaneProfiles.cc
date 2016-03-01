/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kaspar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: ScoringPlaneProfiles.cc,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"

#include "TotemAlignment/RPRates/interface/ScoringPlaneDistributions.h"

#include <map>

#include "TGraph.h"
#include "TFile.h"
  

/**
 *\brief Shows hit distributions (separately for every particle type) recorded in a scoring plane.
 *
 * Designed to work with PPBckgSource.
 **/
class ScoringPlaneProfiles : public edm::EDAnalyzer
{
  public:
    ScoringPlaneProfiles(const edm::ParameterSet &ps);
    ~ScoringPlaneProfiles() {}

  private:
    edm::InputTag HepMCProductLabel;
    unsigned int verbosity;
    unsigned int offsets_N;
    double offsets_from;
    double offsets_to;    
    std::string outputFile;

    std::map<signed int, ScoringPlaneDistributions> hitDistributions;   ///< map: pdg particle id -> hit distributions

    virtual void analyze(const edm::Event &e, const edm::EventSetup &es);
    virtual void endJob();
};

using namespace std;
using namespace edm;
using namespace HepMC;

//----------------------------------------------------------------------------------------------------

ScoringPlaneProfiles::ScoringPlaneProfiles(const edm::ParameterSet &ps) : 
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
  offsets_N(ps.getParameter<unsigned int>("offsets_N")),
  offsets_from(ps.getParameter<double>("offsets_from")),
  offsets_to(ps.getParameter<double>("offsets_to")),
  outputFile(ps.getParameter<string>("outputFile"))
{
	HepMCProductLabel = ps.getParameter<edm::InputTag>("HepMCProductLabel");

}

//----------------------------------------------------------------------------------------------------

void ScoringPlaneProfiles::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  ESHandle<TotemRPGeometry> geometry;
  eSetup.get<RealGeometryRecord>().get(geometry);

  edm::Handle< HepMCProduct > hepMCProd;
  event.getByLabel(HepMCProductLabel, hepMCProd);

  const GenEvent *gEv = hepMCProd->GetEvent();

  for (GenEvent::vertex_const_iterator it = gEv->vertices_begin(); it != gEv->vertices_end(); ++it) {
    GenVertex::particles_out_const_iterator pit = (*it)->particles_out_const_begin();
    const FourVector &v = (*it)->position();  // in mm
    signed int pId = (*pit)->pdg_id();

    if (verbosity > 5)
      printf("* vertex = %E, %E, %E; pId = %i\n", v.x(), v.y(), v.z(), pId);

    map<signed int, ScoringPlaneDistributions>::iterator git = hitDistributions.find(pId);
    if (git == hitDistributions.end()) {
      git = hitDistributions.insert(pair<signed int, ScoringPlaneDistributions>(pId, ScoringPlaneDistributions())).first;
      git->second.Init(&(*geometry), offsets_N, offsets_from, offsets_to, git->first);
    }
    git->second.Fill(v.x(), v.y());
  }
}

//----------------------------------------------------------------------------------------------------

void ScoringPlaneProfiles::endJob()
{
  if (outputFile.empty()) return;

  if (verbosity)
    printf(">> ScoringPlaneProfiles::endJob\n");

  TFile *f = new TFile(outputFile.c_str(), "recreate");
  
  unsigned long total = 0;
  for (map<signed int, ScoringPlaneDistributions>::iterator it = hitDistributions.begin(); it != hitDistributions.end(); ++it)
    total += it->second.eventCounter;
  if (verbosity)
    printf("\tparticles total: %lu\n", total);

  for (map<signed int, ScoringPlaneDistributions>::iterator it = hitDistributions.begin(); it != hitDistributions.end(); ++it) {
    char buf[20];
    sprintf(buf, "particle %i", it->first);
    gDirectory = f->mkdir(buf);
    it->second.Write();
    if (verbosity)
      printf("\tparticle %5i: hits = %5lu, i.e. %5.1f %% from total\n", it->first, it->second.eventCounter, double(it->second.eventCounter) / total * 100.);
  }

  delete f;
}



DEFINE_FWK_MODULE(ScoringPlaneProfiles);

