/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*  Jan Kašpar (jan.kaspar@gmail.com)
*
* $$RCSfile: TransportValidationLibrary.cc,v $: $
* $Revision: 1.7.2.2 $
* $Date: 2009/10/08 16:08:03 $
*
****************************************************************************/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "TotemCondFormats/ProtonTransportFunctions/interface/ProtonTransportFunctions.h"
#include "TotemCondFormats/DataRecord/interface/ProtonTransportRcd.h"
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
#include "TotemRPValidation/TransportValidation/interface/TransportValidationLibrary.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TKey.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TFile.h"

using namespace std;
using namespace edm;
using namespace HepMC;
using namespace CLHEP;
using namespace boost;


//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TransportValidationLibrary::RPInfo::RPInfo() :
  reco_diff_x(new TH1D("reco_diff_x", "reco_diff_x;#Deltax   (mm)", 100, -0.5, 0.5)),
  reco_diff_y(new TH1D("reco_diff_y", "reco_diff_y;#Deltay   (mm)", 100, -0.5, 0.5)),
  simhit_diff_x(new TH1D("simhit_diff_x", "simhit_diff_x;#Deltax   (mm)", 100, -0.5, 0.5)),
  simhit_diff_y(new TH1D("simhit_diff_y", "simhit_diff_y;#Deltay   (mm)", 100, -0.5, 0.5))
{
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TransportValidationLibrary::TransportValidationLibrary(const edm::ParameterSet& conf)
{
  rPFittedTrackCollectionLabel = conf.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
  verbosity = conf.getUntrackedParameter<unsigned int>("verbosity", 0);
  forwardThLimit = conf.getParameter<double>("forwardThLimit");
  outputFile = conf.getParameter<std::string>("outputFile");
}

//----------------------------------------------------------------------------------------------------

TransportValidationLibrary::~TransportValidationLibrary()
{
}

//----------------------------------------------------------------------------------------------------

void TransportValidationLibrary::initialize(const edm::EventSetup&)
{
  TH1::AddDirectory(kFALSE);
}

//----------------------------------------------------------------------------------------------------

void TransportValidationLibrary::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  if (verbosity > 5)
	  printf("\n----------------------------------------------------------------------------------------------------\nEVENT #%i\n", event.id().event());

  // get geometry
  ESHandle<TotemRPGeometry> geometry;
  eSetup.get<RealGeometryRecord>().get(geometry);

  // get optics and beam parameters
  edm::ESHandle<BeamOpticsParams> optPar;
  eSetup.get<BeamOpticsParamsRcd>().get(optPar);
  if(!optPar.isValid())
    throw cms::Exception("TransportValidationLibrary::analyze") << " edm::ESHandle<BeamOpticsParams> is invalid";

  edm::ESHandle<ProtonTransportFunctions> optFun;
  eSetup.get<ProtonTransportRcd>().get(optFun);
  if(!optFun.isValid())
    throw cms::Exception("TransportValidationLibrary::analyze") << " edm::ESHandle<ProtonTransportFunctions> is invalid";

  // get MC event, the SMEARED one
  edm::Handle< HepMCProduct > mcPr;
  HepMC::GenEvent *mcEv = NULL;
  event.getByLabel("generator", "", mcPr);
  mcEv = (HepMC::GenEvent *) mcPr->GetEvent();

  // MC: get vertex
  GenEvent::vertex_iterator vtxIt = mcEv->vertices_begin();
  FourVector vtx = (*vtxIt)->position();    // in mm
  if (++vtxIt != mcEv->vertices_end() && verbosity) printf(">> TransportValidationLibrary::analyze : WARNING: too many vertices in unsmeared MC event.\n");

  // MC: get principal (forward) particles
  HepMC::GenParticle *pt_r = NULL, *pt_l = NULL;
  for (GenEvent::particle_iterator pIt = mcEv->particles_begin(); pIt != mcEv->particles_end(); ++pIt) {
    // skip non-final particles and other than protons
    if ((*pIt)->pdg_id() != 2212 || (*pIt)->status() != 1) continue;

    const FourVector &p = (*pIt)->momentum();
    double th_x = p.x() / p.rho() - optPar->GetCrossingAngleX();
    if (fabs(th_x) > forwardThLimit) continue;
    if (p.z() > 0) {
      if (!pt_r)  pt_r = *pIt;
      else LogWarning("TransportValidationLibrary") << "TransportValidationLibrary::analyze > too many forward particles";
    } else {
      if (!pt_l)  pt_l = *pIt;
      else LogWarning("TransportValidationLibrary") << "TransportValidationLibrary::analyze > too many backward particles";
    }
  }
  if (!pt_r || !pt_l) {
    LogError("TransportValidationLibrary") << "TransportValidationLibrary::analyze > Failed to identify principal particles.";
    return;
  }
  const FourVector &p_r = pt_r->momentum();
  const FourVector &p_l = pt_l->momentum();

  // sim hits to RPIds
  edm::Handle<CrossingFrame<PSimHit> > cf;
  const std::string subdet("g4SimHitsTotemHitsRP");
  event.getByLabel("mix", subdet, cf);
  std::auto_ptr< MixCollection<PSimHit> > simHits(new MixCollection<PSimHit>(cf.product(), std::pair<int,int>(0,0)));
  map< unsigned int, vector<SimHitInfo> > simHitMap;
  MixCollection<PSimHit>::iterator isim;
  for (MixCollection<PSimHit>::iterator isim = simHits->begin(); isim != simHits->end(); ++isim) {
    unsigned int DetId = TotRPDetId::RawToDecId(isim->detUnitId());
    unsigned int RPId = TotRPDetId::RPOfDet(DetId);
    unsigned int det = DetId % 10;
    Hep3Vector loc(isim->localPosition().x(), isim->localPosition().y(), isim->localPosition().z());
    Hep3Vector glb = geometry->LocalToGlobal(isim->detUnitId(), loc);
    simHitMap[RPId].push_back(SimHitInfo(det, glb.x(), glb.y()));
  }

  // traverse all reconstructed hits
  edm::Handle< RPFittedTrackCollection > fitTrCol;
  event.getByLabel(rPFittedTrackCollectionLabel, fitTrCol);
  for (std::map<RPId, RPFittedTrack>::const_iterator it = fitTrCol->begin(); it != fitTrCol->end(); ++it) {
    // skip invalid fits
    if (!it->second.IsValid()) continue;

    // decide which particle to take
    FourVector p;
    if (it->first / 100 == 1)  p = p_r;
    else            p = p_l;

	double dir = p.z() / fabs(p.z());
    double th_x = p.x() / p.rho();
    double th_y = p.y() / p.rho();
	double vt_x = vtx.x() - dir * vtx.z() * th_x;
	double vt_y = vtx.y() - dir * vtx.z() * th_y;
    double xi = p.rho() / optPar->GetBeamMomentum() - 1.;

    double in[5] = {vt_x*1E-3, th_x, vt_y*1E-3, th_y, xi};
    double out[2];  // in m

    LHCOpticsApproximator* opt = optFun->GetFunction(it->first);
    bool valid = opt->Transport2D(in, out, true);
    if (!valid) {
      LogError("TransportValidationLibrary") << "TransportValidationLibrary::analyze > Transport2D invalid";
      continue;
    }

	if (verbosity > 5) {
	  printf("RP %i (everything in mm)\n", it->first);
	  printf("\toptics: x = %.3f, y = %.3f\n", out[0]*1E3, out[1]*1E3);
	  printf("\treco  : x = %.3f, y = %.3f\n", it->second.X0(), it->second.Y0());
	}

    // update histograms
    histograms[it->first].reco_diff_x->Fill(it->second.X0() - out[0]*1E3);
    histograms[it->first].reco_diff_y->Fill(it->second.Y0() - out[1]*1E3);

    if (verbosity > 10)
	  printf("\tsimhits:\n");
    for (unsigned int i = 0; i < simHitMap[it->first].size(); i++) {
      const SimHitInfo &sh = simHitMap[it->first][i];
      histograms[it->first].simhit_diff_x->Fill(sh.x - out[0]*1E3);
      histograms[it->first].simhit_diff_y->Fill(sh.y - out[1]*1E3);
	  if (verbosity > 10) printf("\t\tx = %.3f, y = %.3f\t(%i)\n", sh.x, sh.y, sh.det);
    }
  }
}

//----------------------------------------------------------------------------------------------------

void TransportValidationLibrary::finalize()
{
}

//----------------------------------------------------------------------------------------------------

void TransportValidationLibrary::writeHistogramsToFile()
{
  TFile *of = TFile::Open(outputFile.c_str(), "recreate");
  if(!of || !of->IsWritable())
    {
      std::cout<<"Output file not opened correctly!!"<<std::endl;
    }

  for (map<unsigned int, RPInfo>::iterator it = histograms.begin(); it != histograms.end(); ++it) {
    char buf[10];
    sprintf(buf, "%03u", it->first);
    gDirectory = of->mkdir(buf);
    it->second.reco_diff_x->Write();
    it->second.reco_diff_y->Write();
    it->second.simhit_diff_x->Write();
    it->second.simhit_diff_y->Write();
  }

  delete of;
}

