/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: RPAcceptanceBuilder.cc,v $: $
* $Revision: 1.1.2.1 $
* $Date: 2009/10/15 10:03:39 $
*
****************************************************************************/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TotemRecords/interface/MisalignedGeometryRecord.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrack.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "TH2D.h"
#include "TFile.h"

//----------------------------------------------------------------------------------------------------

/**
 * \brief Builds log t vs. log xi acceptance histograms for each RP.
 *
 * Designed to work with FlatProtonLogKsiLogTGun
 */
class RPAcceptanceBuilder : public edm::EDAnalyzer
{
  public:
    RPAcceptanceBuilder(const edm::ParameterSet &);

  protected:
    edm::InputTag rPFittedTrackCollectionLabel;
    double m;                                 ///< proton mass
    double p_nom;                             ///< nominal beam proton momentum

    unsigned int verbosity;
    std::string hepMCLabel;
    std::string outputFileName;

    struct AcceptancePlots {
      TH2D *t_xi;
      TH1D *t;
      AcceptancePlots() : t_xi(NULL), t(NULL) {}
      AcceptancePlots(const std::string &name);
      void Fill(double vt, double vxi) { t_xi->Fill(vt, vxi); t->Fill(vt); }
    };

    std::map<unsigned int, AcceptancePlots> armHits;
    std::map<unsigned int, AcceptancePlots> rpHits;
    std::map<unsigned int, AcceptancePlots> stHits;

    double CalculateT(double th, double xi);  ///< calculates -t from theta and xi (negative)
    TH2D* MakeStdHist(const char* name);

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event &, const edm::EventSetup&);
    virtual void endJob();
};


using namespace edm;
using namespace std;
using namespace HepMC;

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

RPAcceptanceBuilder::AcceptancePlots::AcceptancePlots(const std::string &name)
{
  t_xi = new TH2D((name + "_t_xi").c_str(), ";log_{10} |t|;log_{10} |#xi|", 53, -4., +1.3, 40, -4., 0.);
  t = new TH1D((name + "_t").c_str(), ";log_{10} |t|", 53, -4., +1.3);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

RPAcceptanceBuilder::RPAcceptanceBuilder(const ParameterSet &iConfig) :
    verbosity(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0)),
    hepMCLabel(iConfig.getUntrackedParameter<string>("hepMCLabel", "source")),
    outputFileName(iConfig.getParameter<string>("outputFileName"))
{
  rPFittedTrackCollectionLabel = iConfig.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
  armHits[0] = AcceptancePlots("left event distribution");
  armHits[1] = AcceptancePlots("right event distribution");
}

//----------------------------------------------------------------------------------------------------

double RPAcceptanceBuilder::CalculateT(double th, double xi)
{
  double omx = 1. + xi;
  double p_sq = p_nom * p_nom;
  double m_sq = m*m;
  double E_sq =  p_sq + m_sq;
  double t0 = 2.*p_sq*xi - 2.*E_sq*(sqrt((p_sq*omx*omx + m_sq) / E_sq) - 1.);
  double S = sin(th/2.);

  return -t0 + 4.*p_sq*omx*S*S;
}

//----------------------------------------------------------------------------------------------------

void RPAcceptanceBuilder::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  ESHandle<TotemRPGeometry> geom;
  es.get<MisalignedGeometryRecord>().get(geom);

  // allocate histograms
  for (TotemRPGeometry::RPDeviceMapType::const_iterator it = geom->beginRP(); it != geom->endRP(); ++it) {
    // rp hists
    char buf[30];
    sprintf(buf, "RP %i: hits", it->first);
    rpHits.insert(pair<unsigned int, AcceptancePlots>(it->first, AcceptancePlots(buf)));

    // station
    unsigned id = it->first / 10;
    if (stHits.find(id) == stHits.end()) {
      sprintf(buf, "station %i: hits", it->first / 10);
      stHits.insert(pair<unsigned int, AcceptancePlots>(id, AcceptancePlots(buf)));
    }
  }

  // get proton mass
  ESHandle<HepPDT::ParticleDataTable> fPDGTable;
  es.getData(fPDGTable);
  const HepPDT::ParticleData* PData = fPDGTable->particle(HepPDT::ParticleID(2212));
  m = PData->mass().value();

  // get nominal momentum
  ESHandle<BeamOpticsParams> optPar;
  es.get<BeamOpticsParamsRcd>().get(optPar);
  p_nom = optPar->GetBeamMomentum();
}

//----------------------------------------------------------------------------------------------------

void RPAcceptanceBuilder::analyze(const edm::Event &event, const EventSetup &es)
{
  // get hepmc event
  Handle<HepMCProduct> mcProd;
  event.getByLabel(hepMCLabel, mcProd);
  const GenEvent *mcEvent = mcProd->GetEvent();

  // extract left and right particles
  GenParticle *partRight = mcEvent->barcode_to_particle(1);
  GenParticle *partLeft = mcEvent->barcode_to_particle(2);

  // calculate right kinematics and fill distribution
  FourVector p = partRight->momentum();
  double xi_r = p.rho() / p_nom - 1.;
  double th_r = p.perp()/p.rho();
  double t_r = CalculateT(th_r, xi_r);
  double ts_r = th_r * p_nom; ts_r *= ts_r;
  //printf("t = %E, xi = %E, ts = %E\n", t_r, xi_r, ts_r);
  //printf("log10 of t = %E, xi = %E\n", log10(t_r), log10(-xi_r));
  armHits[1].Fill(log10(t_r), log10(-xi_r));

  // calculate left kinematics and fill distribution
  p = partLeft->momentum();
  double xi_l = p.rho() / p_nom - 1.;
  double th_l = p.perp()/p.rho();
  double t_l = CalculateT(th_l, xi_l);
  armHits[0].Fill(log10(t_l), log10(-xi_l));

  // get track fits
  Handle<RPFittedTrackCollection> trFits;
  event.getByLabel(rPFittedTrackCollectionLabel, trFits);

  //printf("hit coll size %u\n", trFits->size());
  
  set<unsigned int> hitSt;
  // fill RP hit histograms
  for (RPFittedTrackCollection::const_iterator it = trFits->begin(); it != trFits->end(); ++it) {
    if (!it->second.IsValid())
      continue;

    if (it->first < 100)
      rpHits[it->first].Fill(log10(t_l), log10(-xi_l));
    else
      rpHits[it->first].Fill(log10(t_r), log10(-xi_r));

    hitSt.insert(it->first / 10);
  }

  // fill station hit histograms
  for (set<unsigned int>::iterator it = hitSt.begin(); it != hitSt.end(); ++it) {
    if (*it < 10)
      stHits[*it].Fill(log10(t_l), log10(-xi_l));
    else
      stHits[*it].Fill(log10(t_r), log10(-xi_r));
  }
}

//----------------------------------------------------------------------------------------------------

void RPAcceptanceBuilder::endJob()
{
  TFile *outF = new TFile(outputFileName.c_str(), "recreate");

  // write event distribution
  armHits[0].t->Write();
  armHits[0].t_xi->Write();
  armHits[1].t->Write();
  armHits[1].t_xi->Write();

  for (unsigned int c = 0; c < 2; ++c) {
    map<unsigned int, AcceptancePlots> &collection = (c == 0) ? stHits : rpHits;
    string tag = (c == 0) ? "station" : "RP";
    unsigned int threshold = (c == 0) ? 10 : 100;

    for (map<unsigned int, AcceptancePlots>::const_iterator it = collection.begin(); it != collection.end(); ++it) {
      char buf[30];
      sprintf(buf, "%s %i", tag.c_str(), it->first);
      gDirectory = gDirectory->mkdir(buf);
  
      // save hit distributions
      it->second.t->Write();
      it->second.t_xi->Write();
      
      const AcceptancePlots &reference = armHits[(it->first > threshold) ? 1 : 0];
      
      // build t acceptance
      const TH1D *tnh = reference.t;
      TH1D *tah = new TH1D(* (it->second.t));
      sprintf(buf, "%s %i: acceptance_t", tag.c_str(), it->first);
      tah->SetName(buf);
      for (signed int i = 1; i <= tah->GetNbinsX(); i++) {
        const double norm = tnh->GetBinContent(i);
        if (norm < 1)
          tah->SetBinContent(i, 0.);
        else
          tah->SetBinContent(i, tah->GetBinContent(i) / norm);
      }
      tah->Write();
  
      // build t_xi acceptance
      const TH2D *nh = reference.t_xi;
      TH2D *ah = new TH2D(* (it->second.t_xi));
      sprintf(buf, "%s %i: acceptance_t_xi", tag.c_str(), it->first);
      ah->SetName(buf);
      for (signed int i = 1; i <= ah->GetNbinsX(); i++)
        for (signed int j = 1; j <= ah->GetNbinsY(); j++) {
          const double norm = nh->GetBinContent(i, j);
          if (norm < 2)
            ah->SetBinContent(i, j, 0.);
          else
            ah->SetBinContent(i, j, ah->GetBinContent(i, j) / norm);
        }
      ah->Write();
  
      gDirectory->cd("..");
    }
  }

  delete outF;
}

DEFINE_FWK_MODULE(RPAcceptanceBuilder);

