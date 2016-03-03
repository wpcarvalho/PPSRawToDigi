/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*  Jan Kašpar (jan.kaspar@gmail.com)
*
* $$RCSfile: MCScoringPlaneProfiles.cc,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

// TODO: clean it
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/ProtonTransportFunctions/interface/ProtonTransportFunctions.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "SimG4CMS/TotemRPProtTranspPar/interface/LHCOpticsApproximator.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "Geometry/TotemRPGeometryBuilder/interface/DetGeomDesc.h"
#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TotemAlignment/RPRates/interface/ScoringPlaneDistributions.h"

#include "TKey.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TFile.h"


/**
 *\brief 
 **/
class MCScoringPlaneProfiles : public edm::EDAnalyzer
{
  public:
    MCScoringPlaneProfiles(const edm::ParameterSet &ps);
    ~MCScoringPlaneProfiles() {}

  private:
    unsigned int verbosity;
    double forwardThLimit;
    double forwardXiLimit;
    unsigned int scoringPlaneAtRP;
    unsigned int offsets_N;
    double offsets_from;
    double offsets_to;
    std::string outputFile;

    ScoringPlaneDistributions distributions;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event &e, const edm::EventSetup &es);
    virtual void endJob();
};


using namespace std;
using namespace edm;
using namespace HepMC;

//----------------------------------------------------------------------------------------------------

MCScoringPlaneProfiles::MCScoringPlaneProfiles(const edm::ParameterSet& conf) :
  verbosity(conf.getUntrackedParameter<unsigned int>("verbosity", 0)),
  forwardThLimit(conf.getParameter<double>("forwardThLimit")),
  forwardXiLimit(conf.getParameter<double>("forwardXiLimit")),
  scoringPlaneAtRP(conf.getParameter<unsigned int>("scoringPlaneAtRP")),
  offsets_N(conf.getParameter<unsigned int>("offsets_N")),
  offsets_from(conf.getParameter<double>("offsets_from")),
  offsets_to(conf.getParameter<double>("offsets_to")),
  outputFile(conf.getParameter<std::string>("outputFile"))
{
}

//----------------------------------------------------------------------------------------------------

void MCScoringPlaneProfiles::beginRun(edm::Run const&, edm::EventSetup const& eSetup)
{
  ESHandle<TotemRPGeometry> geometry;
  eSetup.get<RealGeometryRecord>().get(geometry);

  distributions.Init(&(*geometry), offsets_N, offsets_from, offsets_to);
}

//----------------------------------------------------------------------------------------------------

void MCScoringPlaneProfiles::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{

  //printf("--- event: %i --------------------------------------- \n", event.id().event());

  // get optics and beam parameters
  edm::ESHandle<BeamOpticsParams> optPar;
  eSetup.get<BeamOpticsParamsRcd>().get(optPar);
  if(!optPar.isValid())
    throw cms::Exception("MCScoringPlaneProfiles::analyze") << " edm::ESHandle<BeamOpticsParams> is invalid";
  edm::ESHandle<ProtonTransportFunctions> optFun;
  eSetup.get<BeamOpticsParamsRcd>().get(optFun);
  if(!optFun.isValid())
    throw cms::Exception("MCScoringPlaneProfiles::analyze") << " edm::ESHandle<ProtonTransportFunctions> is invalid";

  // get MC event, the SMEARED one
  Handle< HepMCProduct > mcPr;
  event.getByLabel("source", mcPr);
  const GenEvent *mcEv = mcPr->GetEvent();

  // identify left and right forward particles
  GenParticle *p_l = NULL, *p_r = NULL;
  //printf("* loop over particles\n");
  for (GenEvent::particle_const_iterator it = mcEv->particles_begin(); it != mcEv->particles_end(); ++it) {
    // skip non-final particles and other than protons
    if ((*it)->pdg_id() != 2212 || (*it)->status() != 1) continue;

    //printf("\tid = %i, status = %i, px = %E, py = %E, pz = %E\n", (*it)->pdg_id(), (*it)->status(), (*it)->momentum().x(), (*it)->momentum().y(), (*it)->momentum().z());
    
    // skip high scattering angle particles
    const FourVector &p = (*it)->momentum();
    double th_x = p.x() / p.rho() - optPar->GetCrossingAngleX();
    double th_y = p.y() / p.rho();
    if (fabs(th_x) > forwardThLimit || fabs(th_y) > forwardThLimit) continue;

    // skip high xi particles
    double xi = 1. - p.e() / optPar->GetBeamEnergy();
    if (fabs(xi) > forwardXiLimit) continue;

    //printf("\t\tforward\n");

    if (p.z() > 0) {
      if (!p_r)  p_r = *it;
      else LogWarning("MCScoringPlaneProfiles") << "MCScoringPlaneProfiles::analyze > too many forward particles";
    } else {
      if (!p_l)  p_l = *it;
      else LogWarning("MCScoringPlaneProfiles") << "MCScoringPlaneProfiles::analyze > too many backward particles";
    }
  }

  // check existence of the forward particle in the desired direction
  GenParticle *proton = (scoringPlaneAtRP > 100) ? p_r : p_l;
  if (!proton) {
    if (verbosity)
      printf(">> MCScoringPlaneProfiles::analyze > Failed to identify forward particle.\n");
    return;
  }

  // calculate position at z = 0
  const FourVector &p = proton->momentum();
  const FourVector &vtx = proton->production_vertex()->position();
  double dir = p.z() / fabs(p.z());
  double th_x = p.x() / p.rho();
  double th_y = p.y() / p.rho();
  double vt_x = vtx.x() - dir * vtx.z() * th_x;
  double vt_y = vtx.y() - dir * vtx.z() * th_y;
  double xi = p.rho() / optPar->GetBeamMomentum() - 1.;

  // transport
  double in[5] = {vt_x*1E-3, th_x, vt_y*1E-3, th_y, xi};
  if (verbosity)
    printf("* transport input: %E, %E, %E, %E, %E\n", in[0], in[1], in[2], in[3], in[4]);
  double out[2];  // in m
  LHCOpticsApproximator* opt = optFun->GetFunction(scoringPlaneAtRP);

  bool valid = opt->Transport2D(in, out, true);
  if (!valid) {
    if (verbosity)
      printf(">> MCScoringPlaneProfiles::analyze > Transport2D invalid\n");
    return;
  }

  if (verbosity > 5)
	printf("* x = %.3f mm, y = %.3f mm\n", out[0]*1E3, out[1]*1E3);

  // update graphs
  distributions.Fill(out[0]*1E3, out[1]*1E3);
}


//----------------------------------------------------------------------------------------------------

void MCScoringPlaneProfiles::endJob()
{
  if (outputFile.empty())
    return;

  TFile *of = TFile::Open(outputFile.c_str(), "recreate");
  if (!of || !of->IsWritable())
    throw cms::Exception("MCScoringPlaneProfiles::endJob") << "File `" << outputFile << "' is not writable.";

  distributions.Write();

  /*
  for (map<unsigned int, RPInfo>::iterator it = histograms.begin(); it != histograms.end(); ++it) {
    char buf[10];
    sprintf(buf, "%03u", it->first);
    gDirectory = of->mkdir(buf);
    it->second.reco_diff_x->Write();
    it->second.reco_diff_y->Write();
    it->second.simhit_diff_x->Write();
    it->second.simhit_diff_y->Write();
  }
  */

  delete of;
}

DEFINE_FWK_MODULE(MCScoringPlaneProfiles);

