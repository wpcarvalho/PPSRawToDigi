/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Jan Ka??par (jan.kaspar@gmail.com)
*
* $Id: ElasticAcceptanceLibrary.cc 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecoElasticEvent.h"
#include "TotemRPValidation/ElasticReconstruction/interface/ElasticAcceptanceLibrary.h"
#include "TotemRPValidation/ElasticReconstruction/interface/TrackData.h"

using namespace std;
using namespace edm;
using namespace HepMC;


//----------------------------------------------------------------------------------------------------

  unsigned int t_bins = 100;
  double t_from = 0., t_to = 10.;
  char t_label[] = ";|t|   (GeV^{2})";

  unsigned int logt_bins = 100;
  double logt_from = -4., logt_to = 1.;
  char logt_label[] = ";log_{10}|t / GeV^{2}|";

//----------------------------------------------------------------------------------------------------

ElasticAcceptanceLibrary::ElasticAcceptanceLibrary(const edm::ParameterSet& conf) :
  verbosity(conf.getUntrackedParameter<unsigned int>("verbosity", 0)),
  generatorLabel(conf.getParameter<std::string>("generatorLabel")),
  originalLabel(conf.getParameter<std::string>("originalLabel")),
  outputFile(conf.getParameter<std::string>("acceptanceOutputFile"))
{
	rpFittedTrackCollectionLabel = conf.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
	rpRecoElasticEventLabel = conf.getParameter<edm::InputTag>("RPRecoElasticEventLabel");
}

//----------------------------------------------------------------------------------------------------


ElasticAcceptanceLibrary::~ElasticAcceptanceLibrary()
{
}

//----------------------------------------------------------------------------------------------------

ElasticAcceptanceLibrary::AcceptanceCollection::AcceptanceCollection(const char *prefix, bool reco)
{
  char suffix[] = "_s";
  if (reco)
    suffix[1] = 'r';

  char buf[100];
  sprintf(buf, "%s: acceptance_t%s", prefix, suffix);
  t_s = new TH1D(buf, t_label, t_bins, t_from, t_to);
  sprintf(buf, "%s: acceptance_t_o", prefix);
  t_o = new TH1D(buf, t_label, t_bins, t_from, t_to);

  sprintf(buf, "%s: acceptance_tx%s", prefix, suffix);
  tx_s = new TH1D(buf, t_label, t_bins, t_from, t_to);
  sprintf(buf, "%s: acceptance_tx_o", prefix);
  tx_o = new TH1D(buf, t_label, t_bins, t_from, t_to);

  sprintf(buf, "%s: acceptance_ty%s", prefix, suffix);
  ty_s = new TH1D(buf, t_label, t_bins, t_from, t_to);
  sprintf(buf, "%s: acceptance_ty_o", prefix);
  ty_o = new TH1D(buf, t_label, t_bins, t_from, t_to);

  sprintf(buf, "%s: acceptance_logt%s", prefix, suffix);
  logt_s = new TH1D(buf, logt_label, logt_bins, logt_from, logt_to);
  sprintf(buf, "%s: acceptance_logt_o", prefix);
  logt_o = new TH1D(buf, logt_label, logt_bins, logt_from, logt_to);
}

//----------------------------------------------------------------------------------------------------

void ElasticAcceptanceLibrary::initialize(const edm::EventSetup &eSetup)
{
  TH1::AddDirectory(kFALSE);

  total_l = AcceptanceCollection("total_left");
  total_r = AcceptanceCollection("total_right");
  total_reco = AcceptanceCollection("total_reco");

  char buf[100];
  unsigned int id = 0;
  for (unsigned int a = 0; a < 2; a++) {
    id = a;
    sprintf(buf, "arm %u", id); 
    accepted_arm[id] = AcceptanceCollection(buf);

    for (unsigned int s = 0; s < 3; s++) {
      if (s == 1)
        continue;

      id = 10*a + s;
      sprintf(buf, "station %02u", id);
      accepted_station[id] = AcceptanceCollection(buf);

      for (unsigned int u = 0; u < 2; u++) {
        id = 100*a + 10*s + u;
        sprintf(buf, "unit %03u", id);
        accepted_unit[id] = AcceptanceCollection(buf);
      }

      for (unsigned int r = 0; r < 6; r++) {
        id = 100*a + 10*s + r;
        sprintf(buf, "rp %03u", id);
        accepted_rp[id] = AcceptanceCollection(buf);
      }
    }
  }

  accepted_subsystem = AcceptanceCollection("subsystem", true);
}

//----------------------------------------------------------------------------------------------------

void ElasticAcceptanceLibrary::AcceptanceCollection::Fill(const TrackData &tr_o, const TrackData &tr_s)
{
  t_o->Fill(tr_o.t);
  t_s->Fill(tr_s.t);
  tx_o->Fill(tr_o.t_x);
  tx_s->Fill(tr_s.t_x);
  ty_o->Fill(tr_o.t_y);
  ty_s->Fill(tr_s.t_y);
  logt_o->Fill(log10(tr_o.t));
  logt_s->Fill(log10(tr_s.t));
}

//----------------------------------------------------------------------------------------------------

void ElasticAcceptanceLibrary::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  // get optics and beam parameters
  edm::ESHandle<BeamOpticsParams> optPar;
  eSetup.get<BeamOpticsParamsRcd>().get(optPar);
  if(!optPar.isValid())
    throw cms::Exception("ElasticAcceptanceLibrary::analyze") << " edm::ESHandle<BeamOpticsParams> is invalid";

  // get MC event
  if (verbosity > 4) cout << "[0;34mMONTE CARLO[0m" << endl;
  edm::Handle< HepMCProduct > mcPr;
  HepMC::GenEvent *mcEv = NULL, *smearMCEv = NULL, *origMCEv = NULL;

  event.getByLabel(generatorLabel, mcPr);
  mcEv = (HepMC::GenEvent *) mcPr->GetEvent();

  try {
    event.getByLabel("SmearingGenerator", originalLabel, mcPr);
    origMCEv = (HepMC::GenEvent *) mcPr->GetEvent();
    smearMCEv = mcEv;
  }
  catch (...) {
    if (verbosity) printf(">> ElasticAcceptanceLibrary::analyze : WARNING: this event is not smeared\n");
    smearMCEv = origMCEv = mcEv;
  }

  // MC: get principal particles (usualy protons)
  HepMC::GenParticle *part_o_l = NULL, *part_o_r = NULL, *part_s_l = NULL, *part_s_r = NULL;
  if (origMCEv->signal_process_id() == 91) {  // elastic scattering
    part_o_r = origMCEv->barcode_to_particle(3);
    part_o_l = origMCEv->barcode_to_particle(4);
  }
  if (smearMCEv->signal_process_id() == 91) {
    part_s_r = smearMCEv->barcode_to_particle(3);
    part_s_l = smearMCEv->barcode_to_particle(4);
  }
  if (!part_o_l || !part_o_r || !part_s_l || !part_s_r) {
    printf(">> ElasticAcceptanceLibrary::analyze : ERROR: Failed to identify principal particles.\n");
    return;
  }

  // extract kinematical information (original event is ASSUMED to be left-right symmetric)
  TrackData tr_o(part_o_r->momentum(), *optPar, false), tr_s_r(part_s_r->momentum(), *optPar), tr_s_l(part_s_l->momentum(), *optPar);

  // update total histograms
  total_l.Fill(tr_o, tr_s_l);
  total_r.Fill(tr_o, tr_s_r);
  total_reco.Fill(tr_o, tr_o);

  // get list of active elements
  edm::Handle< RPFittedTrackCollection > fitTrColl; 
  event.getByLabel(rpFittedTrackCollectionLabel, fitTrColl);
  set<unsigned int> activeRPs;
  set<unsigned int> activeUnits;
  set<unsigned int> activeStations;
  set<unsigned int> activeArms;
  for (RPFittedTrackCollection::const_iterator it = fitTrColl->begin(); it != fitTrColl->end(); ++it) {
    // skip invalid fits
    if (!it->second.IsValid())
      continue;

    unsigned int rp = it->first;
    unsigned int unit = rp - rp % 10;
    if ((rp%10) > 2)
      unit += 1;
    activeRPs.insert(rp);
    activeUnits.insert(unit);
    activeStations.insert(rp / 10);
    activeArms.insert(rp / 100);
  }
  
  // update accepted histograms
  for (set<unsigned int>::iterator it = activeRPs.begin(); it != activeRPs.end(); ++it)
    accepted_rp[*it].Fill(tr_o, (*it >= 100) ? tr_s_r : tr_s_l);
  for (set<unsigned int>::iterator it = activeUnits.begin(); it != activeUnits.end(); ++it)
    accepted_unit[*it].Fill(tr_o, (*it >= 100) ? tr_s_r : tr_s_l);
  for (set<unsigned int>::iterator it = activeStations.begin(); it != activeStations.end(); ++it)
    accepted_station[*it].Fill(tr_o, (*it >= 10) ? tr_s_r : tr_s_l);
  for (set<unsigned int>::iterator it = activeArms.begin(); it != activeArms.end(); ++it)
    accepted_arm[*it].Fill(tr_o, (*it >= 1) ? tr_s_r : tr_s_l);

  // get elastic reconstruction result
  edm::Handle< RPRecoElasticEvent > elasticReco;
  event.getByLabel(rpRecoElasticEventLabel, elasticReco);
  if (elasticReco->isValid()) {
    TrackData tr_r(elasticReco->result, *optPar);
    accepted_subsystem.Fill(tr_o, tr_r);
  }
}

//----------------------------------------------------------------------------------------------------

void ElasticAcceptanceLibrary::AcceptanceCollection::Normalize(const AcceptanceCollection &total)
{
  t_o->Divide(total.t_o);
  t_s->Divide(total.t_s);
  tx_o->Divide(total.tx_o);
  tx_s->Divide(total.tx_s);
  ty_o->Divide(total.ty_o);
  ty_s->Divide(total.ty_s);
  logt_o->Divide(total.logt_o);
  logt_s->Divide(total.logt_s);
}

//----------------------------------------------------------------------------------------------------

void ElasticAcceptanceLibrary::finalize()
{
  for (map<unsigned int, AcceptanceCollection>::iterator it = accepted_rp.begin(); it != accepted_rp.end(); ++it)
    it->second.Normalize((it->first >= 100) ? total_r : total_l);

  for (map<unsigned int, AcceptanceCollection>::iterator it = accepted_unit.begin(); it != accepted_unit.end(); ++it)
    it->second.Normalize((it->first >= 100) ? total_r : total_l);

  for (map<unsigned int, AcceptanceCollection>::iterator it = accepted_station.begin(); it != accepted_station.end(); ++it)
    it->second.Normalize((it->first >= 10) ? total_r : total_l);

  for (map<unsigned int, AcceptanceCollection>::iterator it = accepted_arm.begin(); it != accepted_arm.end(); ++it)
    it->second.Normalize((it->first >= 1) ? total_r : total_l);

  accepted_subsystem.Normalize(total_reco);
}

//----------------------------------------------------------------------------------------------------

void ElasticAcceptanceLibrary::AcceptanceCollection::Write()
{
  t_o->Write();
  t_s->Write();
  tx_o->Write();
  tx_s->Write();
  ty_o->Write();
  ty_s->Write();
  logt_o->Write();
  logt_s->Write();
}

//----------------------------------------------------------------------------------------------------

void ElasticAcceptanceLibrary::writeHistogramsToFile()
{
  if (outputFile.empty())
    return;

  TFile *of = TFile::Open(outputFile.c_str(), "recreate");
  if(!of || !of->IsWritable()) {
    std::cout << "Output file not opened correctly!!" << std::endl;
    return;
  }

  gDirectory = of->mkdir("total");
  total_l.Write();
  total_r.Write();

  char buf[100];
  
  gDirectory = of->mkdir("rp");
  for (map<unsigned int, AcceptanceCollection>::iterator it = accepted_rp.begin(); it != accepted_rp.end(); ++it) {
    sprintf(buf, "rp %03u", it->first);
    gDirectory = gDirectory->mkdir(buf);
    it->second.Write();
    gDirectory->cd("..");
  }

  gDirectory = of->mkdir("unit");
  for (map<unsigned int, AcceptanceCollection>::iterator it = accepted_unit.begin(); it != accepted_unit.end(); ++it) {
    sprintf(buf, "unit %03u", it->first);
    gDirectory = gDirectory->mkdir(buf);
    it->second.Write();
    gDirectory->cd("..");
  }

  gDirectory = of->mkdir("station");
  for (map<unsigned int, AcceptanceCollection>::iterator it = accepted_station.begin(); it != accepted_station.end(); ++it) {
    sprintf(buf, "station %02u", it->first);
    gDirectory = gDirectory->mkdir(buf);
    it->second.Write();
    gDirectory->cd("..");
  }

  gDirectory = of->mkdir("arm");
  for (map<unsigned int, AcceptanceCollection>::iterator it = accepted_arm.begin(); it != accepted_arm.end(); ++it) {
    sprintf(buf, "arm %01u", it->first);
    gDirectory = gDirectory->mkdir(buf);
    it->second.Write();
    gDirectory->cd("..");
  }

  gDirectory = of->mkdir("subsystem");
  accepted_subsystem.Write();

  delete of;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void ElasticAcceptanceLibrary::ExportAllHistograms()
{
  // TODO
  //ExportHistogram(vtx_x_rs_diff);
}

//----------------------------------------------------------------------------------------------------

void ElasticAcceptanceLibrary::ExportHistogram(TH1D *h)
{
  TH1::AddDirectory(kFALSE);
  TCanvas *c1 = new TCanvas(h->GetName(),h->GetName());
  c1->cd();
  h->Draw();
  elValHists.push_back(c1);
  //  delete c1;
}

