/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
* $Id: ElasticRecoValLibrary.cc 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/

#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"

#include "TotemRPValidation/ElasticReconstruction/interface/ElasticRecoValLibrary.h"
#include "TotemRPValidation/ElasticReconstruction/interface/TrackData.h"

#include "TF1.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TFile.h"

using namespace std;
using namespace edm;
using namespace HepMC;

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TResolution::TResolution(const char* _name, const char* _title, unsigned int _N, double from,
  double to, double range) :
  name(_name), title(_title), N(_N), x0(from), dx((to - from) / N), entries(0), overfill(0),
  underfill(0)
{
  for (unsigned int i = 0; i < N; i++)
    bins.push_back(TMoments(range));
}

//----------------------------------------------------------------------------------------------------

void TResolution::Fill(double x, double y)
{
  if (x < x0) {
    underfill++;
    return;
  }

  unsigned int bin = int(floor((x - x0) / dx));
  if (bin >= N) {
    overfill++;
    return;
  }

  bins[bin].Fill(y);
  entries++;
}

//----------------------------------------------------------------------------------------------------

void TResolution::Write()
{
  printf(">> TResolution::Write > %s\n\tentries = %u\n\tunderfill = %u\n\toverfill = %u\n\ttotal = %u\n", 
      name.c_str(), entries, underfill, overfill, entries+underfill+overfill);
  gDirectory = gDirectory->mkdir(name.c_str());

  // write raw data
  gDirectory = gDirectory->mkdir("profiles");
  for (unsigned int i = 0; i < bins.size(); i++) {
    char str[100];
    sprintf(str, "bin %i (t=%.3f)", i, x0 + (0.5 + i)*dx);
    bins[i].hist->Write(str);
  }
  gDirectory->cd("../");

  // build graphs
  char buf[100], qn[20];
  strcpy(qn, title.c_str());
  TGraphErrors *gr = new TGraphErrors();
  sprintf(buf, ";%s;#sigma(%s^{reco glob} - %s^{orig}) / %s^{orig}", qn, qn, qn, qn);
  gr->SetTitle(buf); gr->SetMarkerStyle(20); gr->SetMarkerSize(0.4);

  TGraphErrors *gm = new TGraphErrors();
  sprintf(buf, ";%s;<%s^{reco glob} - %s^{orig}>", qn, qn, qn);
  gm->SetTitle(buf); gm->SetMarkerStyle(20); gm->SetMarkerSize(0.4);

  TGraphErrors *gs = new TGraphErrors();
  sprintf(buf, ";%s;#sigma(%s^{reco glob} - %s^{orig})", qn, qn, qn);
  gs->SetTitle(buf); gs->SetMarkerStyle(20); gs->SetMarkerSize(0.4);

  double x = x0 + dx/2.;
  for (unsigned int b = 0; b < bins.size(); ++b, x += dx) {
    TMoments& m = bins[b];
    
    // need at least two hits
    if (m.y0 < 2)
      continue;

    double M = m.y1/m.y0;
    double S = sqrt((m.y2 - m.y1 * m.y1 / m.y0) / (m.y0 - 1)); // barlow (5.14)
    double R = S / x;

    double M_err = S / sqrt(m.y0); // barlow (5.11)
    double S_err = S / sqrt(2.*(m.y0 - 1)); // barlow (5.23)
    double x_err = dx/sqrt(12.);
    double R_err = sqrt(pow(S_err/x, 2.) + pow(S/x/x*x_err, 2.));

    int point = gr->GetN();

    gm->SetPoint(point, x, M);
    gm->SetPointError(point, 0., M_err);

    gs->SetPoint(point, x, S);
    gs->SetPointError(point, 0., S_err);

    gr->SetPoint(point, x, R);
    gr->SetPointError(point, 0., R_err);

    printf("\t%i\tx = %E\tN = %.0f\tmean = %E\tsigma = %E\n", b, x, m.y0, M, S);
  }

  // fit resolution curve
  if (gr->GetN() > 5) {
    TF1 *fitFun = new TF1("user", "[0]/sqrt(x)");
    gr->Fit(fitFun, "Q");
    printf("\tA = %.2E\n", fitFun->GetParameter(0));
  }

  // write the graphs
  gm->Write("mean profile");
  gs->Write("sigma profile");
  gr->Write("resolution");
  gDirectory->cd("../");
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

ElasticRecoValLibrary::RecoCounter::RecoCounter() : tot(0), trig(0), rOK(0), rNoRoad(0), rNoGoodRoad(0),
  rOnlyOneArm(0), rBadVertexX(0), rBadVertexY(0), rBadAngleX(0), rBadAngleY(0)
{
  th_x = new TH1D("", ";#Delta#vartheta_{x}", 200, -3E-5, 3E-5);
  th_y = new TH1D("", ";#Delta#vartheta_{y}", 200, -3E-6, 3E-6);
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::RecoCounter::Update(bool triggered, const RPRecoElasticEvent &ee)
{
  tot++;
  if (triggered) trig++;
  if (ee.status == ee.sOK) rOK++;
  if (ee.status == ee.sNoRoad) rNoRoad++;
  if (ee.status == ee.sNoGoodRoad) rNoGoodRoad++;
  if (ee.rejectReason & ee.rVertexX) rBadVertexX++;
  if (ee.rejectReason & ee.rVertexY) rBadVertexY++;
  if (ee.rejectReason & ee.rAngleX) { rBadAngleX++; th_x->Fill(ee.rightFit.th_x - ee.leftFit.th_x); }
  if (ee.rejectReason & ee.rAngleY) { rBadAngleY++; th_y->Fill(ee.rightFit.th_y - ee.leftFit.th_y); }

  if (ee.status == ee.sNoGoodRoad) {
    bool left = false, right = false;
    for (unsigned int i = 0; i < ee.roads.size(); i++) {
      for (unsigned int j = 0; j < ee.roads[i].members.size(); j++) {
        if (ee.roads[i].members[j] / 100 == 1)  right = true;
        else                  left = true;
      }
    }
    if (!left || !right) rOnlyOneArm++;
  }
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

ElasticRecoValLibrary::FitStatistic::FitStatistic(const string &_name) : name(_name)
{
  string label;
  label = "ndf_x_" + name;
  ndf_x = new TH1D(label.c_str(), ";number of degrees of freedom, x fit", 23, -0.5, 22.5);
  label = "ndf_y_" + name;
  ndf_y = new TH1D(label.c_str(), ";number of degrees of freedom, y fit", 23, -0.5, 22.5);
  label = "p_x_" + name;
  p_x = new TH1D(label.c_str(), "; p value, x fit", 100, 0., 1.); 
  label = "p_y_" + name;
  p_y = new TH1D(label.c_str(), "; p value, y fit", 100, 0., 1.); 
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::FitStatistic::Fill(unsigned int _ndf_x, unsigned int _ndf_y, double s2_min_x, double s2_min_y)
{
  ndf_x->Fill(_ndf_x);
  ndf_y->Fill(_ndf_y);

  map<unsigned int, TH1D*>::iterator it;
  if (_ndf_x > 0) {
    p_x->Fill(TMath::Prob(s2_min_x, _ndf_x));
    it = p_ndf_x.find(_ndf_x);
    if (it == p_ndf_x.end()) {
      char buf[50];
      sprintf(buf, "p_x_%s_ndf=%i", name.c_str(), _ndf_x);
      it = p_ndf_x.insert(pair<unsigned int, TH1D*> (_ndf_x, new TH1D(buf, "; p value, x fit", 100, 0., 1.))).first;
    }
    it->second->Fill(TMath::Prob(s2_min_x, _ndf_x));
  }

  if (_ndf_y > 0) {
    p_y->Fill(TMath::Prob(s2_min_y, _ndf_y));
    it = p_ndf_y.find(_ndf_y);
    if (it == p_ndf_y.end()) {
      char buf[50];
      sprintf(buf, "p_y_%s_ndf=%i", name.c_str(), _ndf_y);
      it = p_ndf_y.insert(pair<unsigned int, TH1D*> (_ndf_y, new TH1D(buf, "; p value, y fit", 100, 0., 1.))).first;
    }
    it->second->Fill(TMath::Prob(s2_min_y, _ndf_y));
  }
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::FitStatistic::Write()
{
  ndf_x->Write();
  ndf_y->Write();
  p_x->Write();
  p_y->Write();
  for (map<unsigned int, TH1D*>::iterator it = p_ndf_x.begin(); it != p_ndf_x.end(); ++it)
    it->second->Write();
  for (map<unsigned int, TH1D*>::iterator it = p_ndf_y.begin(); it != p_ndf_y.end(); ++it)
    it->second->Write();
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

ElasticRecoValLibrary::ElasticRecoValLibrary(const edm::ParameterSet& conf) :
  verbosity(conf.getUntrackedParameter<unsigned int>("verbosity", 0)),
  generatorLabel(conf.getParameter<std::string>("generatorLabel")),
  originalLabel(conf.getParameter<std::string>("originalLabel")),
  outputFile(conf.getParameter<std::string>("validationOutputFile")),
  



  lim_vtx_dx(conf.getUntrackedParameter<double>("lim_vtx_dx", 0.)),
  lim_vtx_dy(conf.getUntrackedParameter<double>("lim_vtx_dy", 0.)),
  lim_th_dx(conf.getUntrackedParameter<double>("lim_th_dx", 0.)),
  lim_th_dy(conf.getUntrackedParameter<double>("lim_th", 0.)),
  lim_t_min(conf.getUntrackedParameter<double>("lim_t_min", 0.)),
  lim_t_max(conf.getUntrackedParameter<double>("lim_t_max", 0.)),

  fStat_left("left"),
  fStat_right("right"),
  fStat_global("global")
{
	  rpDetTriggerSetLabel = conf.getParameter<edm::InputTag>("RPDetTriggerSetLabel");
	  rpFittedTrackCollectionLabel = conf.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
	  rpRecoElasticEventLabel = conf.getParameter<edm::InputTag>("RPRecoElasticEventLabel");
}

//----------------------------------------------------------------------------------------------------


ElasticRecoValLibrary::~ElasticRecoValLibrary()
{
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::initialize(const edm::EventSetup &eSetup)
{
  // get optics and beam parameters
  edm::ESHandle<BeamOpticsParams> optPar;
  eSetup.get<BeamOpticsParamsRcd>().get(optPar);
  if(!optPar.isValid())
    throw cms::Exception("ElasticRecoValLibrary::analyze") << " edm::ESHandle<BeamOpticsParams> is invalid";

  TH1::AddDirectory(kFALSE);

  double range_x, range_y;

  // init vertex plots
  range_x = (lim_vtx_dx == 0.) ? 1E-3 : lim_vtx_dx;
  range_y = (lim_vtx_dy == 0.) ? 1E-3 : lim_vtx_dy;
  vtx_x_rs_diff = new TH1D("vtx_x_rs_diff", ";x*^{reco} - x*^{smear}   (m)", 100, -range_x, range_x);
  vtx_y_rs_diff = new TH1D("vtx_y_rs_diff", ";y*^{reco} - y*^{smear}   (m)", 100, -range_y, range_y);
  vtx_x_ro_diff = new TH1D("vtx_x_ro_diff", ";x*^{reco} - x*^{orig}   (m)", 100, -range_x, range_x);
  vtx_y_ro_diff = new TH1D("vtx_y_ro_diff", ";y*^{reco} - y*^{orig}   (m)", 100, -range_y, range_y);

  vtx_x_rs_corr = new TGraph(); vtx_x_rs_corr->SetName("vtx_x_rs_corr"); vtx_x_rs_corr->SetTitle(";x*^{smear}   (m);x*^{reco}   (m)");
  vtx_y_rs_corr = new TGraph(); vtx_y_rs_corr->SetName("vtx_y_rs_corr"); vtx_y_rs_corr->SetTitle(";y*^{smear}   (m);y*^{reco}   (m)");
  vtx_x_ro_corr = new TGraph(); vtx_x_ro_corr->SetName("vtx_x_ro_corr"); vtx_x_ro_corr->SetTitle(";x*^{orig}   (m);x*^{reco}   (m)");
  vtx_y_ro_corr = new TGraph(); vtx_y_ro_corr->SetName("vtx_y_ro_corr"); vtx_y_ro_corr->SetTitle(";y*^{orig}   (m);y*^{reco}   (m)");

  // init theta plots
  range_x = (lim_th_dx == 0.) ? 5*optPar->GetBeamDivergenceX() : lim_th_dx;
  range_y = (lim_th_dy == 0.) ? 5*optPar->GetBeamDivergenceY() : lim_th_dy;
  th_rlsl_diff = new TH1D("th_rlsl_diff", ";#vartheta^{reco left} - #vartheta^{smear left}   (rad)", 500, -range_x, range_x);
  th_x_rlsl_diff = new TH1D("th_x_rlsl_diff", ";#vartheta_{x}^{reco left} - #vartheta_{x}^{smear left}   (rad)", 500, -range_x, range_x);
  th_y_rlsl_diff = new TH1D("th_y_rlsl_diff", ";#vartheta_{y}^{reco left} - #vartheta_{y}^{smear left}   (rad)", 500, -range_y, range_y);
  th_x_rlo_diff = new TH1D("th_x_rlo_diff", ";#vartheta_{x}^{reco left} - #vartheta_{x}^{orig left}   (rad)", 500, -range_x, range_x);
  th_y_rlo_diff = new TH1D("th_y_rlo_diff", ";#vartheta_{y}^{reco left} - #vartheta_{y}^{orig left}   (rad)", 500, -range_y, range_y);
  th_x_slo_diff = new TH1D("th_x_slo_diff", ";#vartheta_{x}^{smear left} - #vartheta_{x}^{orig left}   (rad)", 500, -range_x, range_x);
  th_y_slo_diff = new TH1D("th_y_slo_diff", ";#vartheta_{y}^{smear left} - #vartheta_{y}^{orig left}   (rad)", 500, -range_y, range_y);
  th_rrsr_diff = new TH1D("th_rrsr_diff", ";#vartheta^{reco right} - #vartheta^{smear right}   (rad)", 500, -range_x, range_x);
  th_x_rrsr_diff = new TH1D("th_x_rrsr_diff", ";#vartheta_{x}^{reco right} - #vartheta_{x}^{smear right}   (rad)", 500, -range_x, range_x);
  th_y_rrsr_diff = new TH1D("th_y_rrsr_diff", ";#vartheta_{y}^{reco right} - #vartheta_{y}^{smear right}   (rad)", 500, -range_y, range_y);
  th_rro_diff = new TH1D("th_rro_diff", ";#vartheta^{reco right} - #vartheta^{orig right}   (rad)", 500, -range_x, range_x);
  th_rlo_diff = new TH1D("th_rlo_diff", ";#vartheta^{reco left} - #vartheta^{orig left}   (rad)", 500, -range_x, range_x);
  th_x_rro_diff = new TH1D("th_x_rro_diff", ";#vartheta_{x}^{reco right} - #vartheta_{x}^{orig right}   (rad)", 500, -range_x, range_x);
  th_y_rro_diff = new TH1D("th_y_rro_diff", ";#vartheta_{y}^{reco right} - #vartheta_{y}^{orig right}   (rad)", 500, -range_y, range_y);
  th_x_sro_diff = new TH1D("th_x_sro_diff", ";#vartheta_{x}^{smear right} - #vartheta_{x}^{orig right}   (rad)", 500, -range_x, range_x);
  th_y_sro_diff = new TH1D("th_y_sro_diff", ";#vartheta_{y}^{smear right} - #vartheta_{y}^{orig right}   (rad)", 500, -range_y, range_y);

  th_rgo_diff = new TH1D("th_rgo_diff", ";#vartheta^{reco glob} - #vartheta^{orig}   (rad)", 500, -range_x, range_x);
  th_x_rgo_diff = new TH1D("th_x_rgo_diff", ";#vartheta_{x}^{reco glob} - #vartheta_{x}^{orig}   (rad)", 500, -range_x, range_x);
  th_y_rgo_diff = new TH1D("th_y_rgo_diff", ";#vartheta_{y}^{reco glob} - #vartheta_{y}^{orig}   (rad)", 500, -range_y, range_y);

  th_rgo_corr = new TGraph(); th_rgo_corr->SetTitle(";#theta^{orig};#theta^{reco}");
  th_x_rgo_corr = new TGraph(); th_x_rgo_corr->SetTitle(";#theta_{x}^{orig};#theta_{x}^{reco}");
  th_y_rgo_corr = new TGraph(); th_y_rgo_corr->SetTitle(";#theta_{y}^{orig};#theta_{y}^{reco}");

  th_rgo_res = new TResolution("th_rgo_res", "#vartheta", 200, -3E-4, +3E-4, range_x);
  th_x_rgo_res = new TResolution("th_x_rgo_res", "#vartheta_{x}", 200, -3E-4, +3E-4, range_x);
  th_y_rgo_res = new TResolution("th_y_rgo_res", "#vartheta_{y}", 200, -3E-4, +3E-4, range_y);

  // init t plots
  double range = 2*optPar->GetBeamMomentum()*optPar->GetBeamDivergenceX() * 1 * 5; // typical |t| = 1 Gev^2
  t_rlsl_diff = new TH1D("t_rlsl_diff", ";t^{reco left} - t^{smear left}   (GeV^{2})", 500, -range, range);
  t_x_rlsl_diff = new TH1D("t_x_rlsl_diff", ";t^{reco left}_{x} - t^{smear left}_{x}   (GeV^{2})", 500, -range, range);
  t_y_rlsl_diff = new TH1D("t_y_rlsl_diff", ";t^{reco left}_{y} - t^{smear left}_{y}   (GeV^{2})", 500, -range, range);

  t_rrsr_diff = new TH1D("t_rrsr_diff", ";t^{reco right} - t^{smear right}   (GeV^{2})", 500, -range, range);
  t_x_rrsr_diff = new TH1D("t_x_rrsr_diff", ";t^{reco right}_{x} - t^{smear right}_{x}   (GeV^{2})", 500, -range, range);
  t_y_rrsr_diff = new TH1D("t_y_rrsr_diff", ";t^{reco right}_{y} - t^{smear right}_{y}   (GeV^{2})", 500, -range, range);

  t_rgo_diff = new TH1D("t_rgo_diff", ";t^{reco glob} - t^{orig}   (GeV^{2})", 500, -range, range);
  t_x_rgo_diff = new TH1D("t_x_rgo_diff", ";t^{reco glob}_{x} - t^{orig}_{x}   (GeV^{2})", 500, -range, range);
  t_y_rgo_diff = new TH1D("t_y_rgo_diff", ";t^{reco glob}_{y} - t^{orig}_{y}   (GeV^{2})", 500, -range, range);

  int N = 50;
  t_rgo_res = new TResolution("t_rgo_res", "|t|", N, lim_t_min, lim_t_max, range);
  t_x_rgo_res = new TResolution("t_x_rgo_res", "|t_{y}|", N, lim_t_min, lim_t_max, range);
  t_y_rgo_res = new TResolution("t_y_rgo_res", "|t_{x}|", N, lim_t_min, lim_t_max, range);

  // init distribution plots
  t_rg_dist = new TH1D("t_rg_dist", ";t^{reco}   (GeV^{2})", 100, 0, 3);
  t_sl_dist = new TH1D("t_sl_dist", ";t^{smear left}   (GeV^{2})", 100, 0, 3);
  t_sr_dist = new TH1D("t_sr_dist", ";t^{smear right}   (GeV^{2})", 100, 0, 3);
  t_o_dist = new TH1D("t_o_dist", ";t^{orig}   (GeV^{2})", 100, 0, 3);

  th_rg_dist = new TH1D("th_rg_dist",  ";#vartheta^{reco}   (rad)", 100, -1E-3, 1E-3);
  th_sl_dist = new TH1D("th_sl_dist",  ";#vartheta^{smear left}   (rad)", 100, -1E-3, 1E-3);
  th_sr_dist = new TH1D("th_sr_dist",  ";#vartheta^{smear right}   (rad)", 100, -1E-3, 1E-3);
  th_o_dist = new TH1D("th_o_dist",  ";#vartheta^{orig}   (rad)", 100, -1E-3, 1E-3);

  phi_rg_dist = new TH1D("phi_rg_dist", ";#phi^{reco}   (rad)", 100, -M_PI, M_PI);
  phi_sl_dist = new TH1D("phi_sl_dist", ";#phi^{smear left}   (rad)", 100, -M_PI, M_PI);
  phi_sr_dist = new TH1D("phi_sr_dist", ";#phi^{smear right}   (rad)", 100, -M_PI, M_PI);
  phi_o_dist = new TH1D("phi_o_dist", ";#phi^{orig}   (rad)", 100, -M_PI, M_PI);

  // init left-right diffrence plots
  range = 9*optPar->GetBeamDivergenceX();
  th_x_r_rldiff = new TH1D("th_x_r_rldiff", ";#Delta_{R-L}#vartheta_{x}^{reco}   (rad)", 500, -range, range);
  th_y_r_rldiff = new TH1D("th_y_r_rldiff", ";#Delta_{R-L}#vartheta_{y}^{reco}   (rad)", 500, -range, range);
  vtx_x_r_rldiff = new TH1D("vtx_x_r_rldiff", ";#Delta_{R-L}x*^{reco}   (m)", 100, -1E-2, 1E-2);
  vtx_y_r_rldiff = new TH1D("vtx_y_r_rldiff", ";#Delta_{R-L}y*^{reco}   (m)", 100, -1E-2, 1E-2);
  th_x_s_rldiff = new TH1D("th_x_s_rldiff", ";#Delta_{R-L}#vartheta_{x}^{smear}   (rad)", 500, -range, range);
  th_y_s_rldiff = new TH1D("th_y_s_rldiff", ";#Delta_{R-L}#vartheta_{y}^{smear}   (rad)", 500, -range, range);

  // statistics
  aRPCount_left = new TH1D("aRPCount_left", ";count of active left RP", 13, -0.5, 12.5);
  aRPCount_right = new TH1D("aRPCount_right", ";count of active right RP", 13, -0.5, 12.5);
  aRPCount_global = new TH1D("aRPCount_global", ";count of active global RP", 25, -0.5, 24.5);

  // road-size histograms
  rs_x = new TH1D("rs_x", ";road-size x", 500, 0, 20E-6);
  rs_y = new TH1D("rs_y", ";road-size y", 500, 0, 20E-6);

#ifdef _ER_DEBUG_
  debug1 = new TGraph(); debug1->SetName("DeThR vx DeThL");
  range = 1e-3;
  debug2 = new TH1D("th_x_detec", ";#theta_x when detected   (rad)", 500, -range, range);
  debug3 = new TH1D("th_y_detec", ";#theta_y when detected   (rad)", 500, -range, range);
#endif
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  if (verbosity > 4) cout << endl << "[0;32mEVENT[0m " << event.id().event() << endl;

  // get optics and beam parameters
  edm::ESHandle<BeamOpticsParams> optPar;
  eSetup.get<BeamOpticsParamsRcd>().get(optPar);
  if(!optPar.isValid())
    throw cms::Exception("ElasticRecoValLibrary::analyze") << " edm::ESHandle<BeamOpticsParams> is invalid";

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
    if (verbosity) printf(">> ElasticRecoValLibrary::analyze : WARNING: this event is not smeared\n");
    smearMCEv = origMCEv = mcEv;
  }

  // MC: get vertices
  GenEvent::vertex_iterator vtxIt = origMCEv->vertices_begin();
  FourVector vtx_o = (*vtxIt)->position();    // in mm
  if (++vtxIt != origMCEv->vertices_end()) printf(">> ElasticRecoValLibrary::analyze : WARNING: too many vertices in unsmeared MC event.\n");
  vtxIt = smearMCEv->vertices_begin();
  FourVector vtx_s = (*vtxIt)->position();    // in mm
  if (++vtxIt != smearMCEv->vertices_end()) printf(">> ElasticRecoValLibrary::analyze : WARNING: too many vertices in smeared MC event.\n");
  if (verbosity > 4) {
    printf("> vertex:\n\toriginal: %E, %E, %E (mm)\n\tsmeared: %E, %E, %E (mm)\n", vtx_o.x(), vtx_o.y(), vtx_o.z(), vtx_s.x(), vtx_s.y(), vtx_s.z());
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
    printf(">> ElasticRecoValLibrary::analyze : ERROR: Failed to identify principal particles.\n");
    return;
  }

  // extract kinematical information (original event is ASSUMED to be left-right symmetric)
  TrackData tr_o(part_o_r->momentum(), *optPar, false), tr_s_r(part_s_r->momentum(), *optPar), tr_s_l(part_s_l->momentum(), *optPar);
  if (verbosity > 4) {
    printf("> tracks:\n\tsource\t\tt (GeV^2)\tphi\tth (murad)\tth_x\tth_y\n");
    printf("\toriginal\t%.4f\t\t%.3f\t%.1f\t\t%.1f\t%.1f\n", tr_o.t, tr_o.phi, tr_o.th*1E6, tr_o.th_x*1E6, tr_o.th_y*1E6);
    printf("\tsmeared-right\t%.4f\t\t%.3f\t%.1f\t\t%.1f\t%.1f\n", tr_s_r.t, tr_s_r.phi, tr_s_r.th*1E6, tr_s_r.th_x*1E6, tr_s_r.th_y*1E6);
    printf("\tsmeared-left\t%.4f\t\t%.3f\t%.1f\t\t%.1f\t%.1f\n", tr_s_l.t, tr_s_l.phi, tr_s_l.th*1E6, tr_s_l.th_x*1E6, tr_s_l.th_y*1E6);
  }

  // update reco independent plots
  t_sl_dist->Fill(tr_s_l.t);
  t_sr_dist->Fill(tr_s_r.t);
  t_o_dist->Fill(tr_o.t);

  th_x_sro_diff->Fill(tr_s_r.th_x - tr_o.th_x);
  th_y_sro_diff->Fill(tr_s_r.th_y - tr_o.th_y);
  th_x_slo_diff->Fill(tr_s_l.th_x - tr_o.th_x);
  th_y_slo_diff->Fill(tr_s_l.th_y - tr_o.th_y);

  th_sl_dist->Fill(tr_s_l.th);
  th_sr_dist->Fill(tr_s_r.th);
  th_o_dist->Fill(tr_o.th);

  phi_sl_dist->Fill(tr_s_l.phi);
  phi_sr_dist->Fill(tr_s_r.phi);
  phi_o_dist->Fill(tr_o.phi);

  th_x_s_rldiff->Fill(tr_s_r.th_x - tr_s_l.th_x);
  th_y_s_rldiff->Fill(tr_s_r.th_y - tr_s_l.th_y);

  // event triggered?
  edm::Handle< DetSetVector<RPDetTrigger> > trigCol;
  event.getByLabel(rpDetTriggerSetLabel, trigCol);
  bool leftHit = false, rightHit = false;
  for (DetSetVector<RPDetTrigger>::const_iterator it = trigCol->begin(); it != trigCol->end(); ++it) {
    unsigned int arm = TotRPDetId::ArmOfDet(TotRPDetId::RawToDecId(it->detId()));
    if (it->size() && arm == 0) leftHit = true;
    if (it->size() && arm == 1) rightHit = true;
    if (leftHit && rightHit) break;
  }
  bool triggered = leftHit && rightHit;

  // count active RPs
  edm::Handle< RPFittedTrackCollection > fitTrColl; 
  event.getByLabel(rpFittedTrackCollectionLabel, fitTrColl);
  unsigned int acLeftRPs = 0, acRightRPs = 0;
  for (RPFittedTrackCollection::const_iterator it = fitTrColl->begin(); it != fitTrColl->end(); ++it) {
    unsigned int arm = TotRPDetId::ArmOfRP(it->first);
    if (arm == 0) acLeftRPs++;
    if (arm == 1) acRightRPs++;
  }
  aRPCount_left->Fill(acLeftRPs);
  aRPCount_right->Fill(acRightRPs);
  aRPCount_global->Fill(acLeftRPs + acRightRPs);

  // read/print elastic report
  edm::Handle< RPRecoElasticEvent > elasticReco;
  event.getByLabel(rpRecoElasticEventLabel, elasticReco);
  if (verbosity > 4) {
    cout << "[0;34mELASTIC RECONSTRUCTION[0m" << endl;
    printf("status = ");
    switch (elasticReco->status) {
      case RPRecoElasticEvent::sOK: printf("OK\n"); break;
      case RPRecoElasticEvent::sNoRoad: printf("no road\n"); break;
      case RPRecoElasticEvent::sNoGoodRoad: printf("no good road\n"); break;
      case RPRecoElasticEvent::sRejected:
        printf("rejected because of");
        if (elasticReco->rejectReason & RPRecoElasticEvent::rVertexX) printf("vertex x, ");
        if (elasticReco->rejectReason & RPRecoElasticEvent::rVertexY) printf("vertex y, ");
        if (elasticReco->rejectReason & RPRecoElasticEvent::rAngleX) printf("theta x, ");
        if (elasticReco->rejectReason & RPRecoElasticEvent::rAngleY) printf("theta y");
        printf("\n");
        break;
      default: printf("undefined\n");
    }
  }
  if (verbosity > 10) PrintReport(* elasticReco);

  // increment counters
  counters[part_s_r->pdg_id()].Update(triggered, *elasticReco);

  // stop if no valid elastic track
  if (!elasticReco->isValid())
    return;

  // update vertex plots
  vtx_x_rs_diff->Fill(elasticReco->result.x - vtx_s.x()*1E-3);  // vtx: from mm to m
  vtx_y_rs_diff->Fill(elasticReco->result.y - vtx_s.y()*1E-3);
  vtx_x_ro_diff->Fill(elasticReco->result.x - vtx_o.x()*1E-3);
  vtx_y_ro_diff->Fill(elasticReco->result.y - vtx_o.y()*1E-3);

  vtx_x_rs_corr->SetPoint(vtx_x_rs_corr->GetN(), vtx_s.x()*1E-3, elasticReco->result.x);
  vtx_y_rs_corr->SetPoint(vtx_y_rs_corr->GetN(), vtx_s.y()*1E-3, elasticReco->result.y);
  vtx_x_ro_corr->SetPoint(vtx_x_ro_corr->GetN(), vtx_o.x()*1E-3, elasticReco->result.x);
  vtx_y_ro_corr->SetPoint(vtx_y_ro_corr->GetN(), vtx_o.t()*1E-3, elasticReco->result.y);

  // update theta plots
  TrackData tr_r_g(elasticReco->result, *optPar), tr_r_l(elasticReco->leftFit, *optPar), tr_r_r(elasticReco->rightFit, *optPar);

  th_rlsl_diff->Fill(tr_r_l.th - tr_s_l.th);
  th_x_rlsl_diff->Fill(tr_r_l.th_x - tr_s_l.th_x);
  th_y_rlsl_diff->Fill(tr_r_l.th_y - tr_s_l.th_y);
  th_x_rlo_diff->Fill(tr_r_l.th_x - tr_o.th_x);
  th_y_rlo_diff->Fill(tr_r_l.th_y - tr_o.th_y);
  th_rlo_diff->Fill(tr_r_l.th - tr_o.th);

  th_rrsr_diff->Fill(tr_r_r.th - tr_s_r.th);
  th_x_rrsr_diff->Fill(tr_r_r.th_x - tr_s_r.th_x);
  th_y_rrsr_diff->Fill(tr_r_r.th_y - tr_s_r.th_y);
  th_x_rro_diff->Fill(tr_r_r.th_x - tr_o.th_x);
  th_y_rro_diff->Fill(tr_r_r.th_y - tr_o.th_y);
  th_rro_diff->Fill(tr_r_r.th - tr_o.th);

  th_rgo_diff->Fill(tr_r_g.th - tr_o.th);
  th_x_rgo_diff->Fill(tr_r_g.th_x - tr_o.th_x);
  th_y_rgo_diff->Fill(tr_r_g.th_y - tr_o.th_y);

  th_rgo_corr->SetPoint(th_rgo_corr->GetN(), tr_r_g.th, tr_o.th);
  th_x_rgo_corr->SetPoint(th_x_rgo_corr->GetN(), tr_r_g.th_x, tr_o.th_x);
  th_y_rgo_corr->SetPoint(th_y_rgo_corr->GetN(), tr_r_g.th_y, tr_o.th_y);

  th_rgo_res->Fill(tr_o.th, tr_r_g.th - tr_o.th);
  th_x_rgo_res->Fill(tr_o.th_x, tr_r_g.th_x - tr_o.th_x);
  th_y_rgo_res->Fill(tr_o.th_y, tr_r_g.th_y - tr_o.th_y);

  // update t plots
  t_rlsl_diff->Fill(tr_r_l.t - tr_s_l.t);
  t_x_rlsl_diff->Fill(tr_r_l.t_x - tr_s_l.t_x);
  t_y_rlsl_diff->Fill(tr_r_l.t_y - tr_s_l.t_y);

  t_rrsr_diff->Fill(tr_r_r.t - tr_s_r.t);
  t_x_rrsr_diff->Fill(tr_r_r.t_x - tr_s_r.t_x);
  t_y_rrsr_diff->Fill(tr_r_r.t_y - tr_s_r.t_y);

  t_rgo_diff->Fill(tr_r_g.t - tr_o.t);
  t_x_rgo_diff->Fill(tr_r_g.t_x - tr_o.t_x);
  t_y_rgo_diff->Fill(tr_r_g.t_y - tr_o.t_y);

  t_rgo_res->Fill(tr_o.t, tr_r_g.t - tr_o.t);
  t_x_rgo_res->Fill(tr_o.t_x, tr_r_g.t_x - tr_o.t_x);
  t_y_rgo_res->Fill(tr_o.t_y, tr_r_g.t_y - tr_o.t_y);

  // update (the remaining) distribution plots
  t_rg_dist->Fill(tr_r_g.t);
  th_rg_dist->Fill(tr_r_g.th);
  phi_rg_dist->Fill(tr_r_g.phi);

  // update (the remaining) left-right diffrence plots
  th_x_r_rldiff->Fill(tr_r_r.th_x - tr_r_l.th_x);
  th_y_r_rldiff->Fill(tr_r_r.th_y - tr_r_l.th_y);
  vtx_x_r_rldiff->Fill(elasticReco->rightFit.x - elasticReco->leftFit.x);
  vtx_y_r_rldiff->Fill(elasticReco->rightFit.y - elasticReco->leftFit.y);

  // update statistics plots
  fStat_left.Fill(elasticReco->leftFit.ndf_x, elasticReco->leftFit.ndf_y, elasticReco->leftFit.s2min_x, elasticReco->leftFit.s2min_y);
  fStat_right.Fill(elasticReco->rightFit.ndf_x, elasticReco->rightFit.ndf_y, elasticReco->rightFit.s2min_x, elasticReco->rightFit.s2min_y);
  fStat_global.Fill(elasticReco->globalFit.ndf_x, elasticReco->globalFit.ndf_y, elasticReco->globalFit.s2min_x, elasticReco->globalFit.s2min_y);

  // road-size histograms
  rs_x->Fill(elasticReco->roads[elasticReco->preferredRoad].SizeX());
  rs_y->Fill(elasticReco->roads[elasticReco->preferredRoad].SizeY());

#ifdef _ER_DEBUG_
  debug1->SetPoint(debug1->GetN(), tr_r_r.th_x - tr_o.th_x, tr_r_l.th_x - tr_o.th_x);
  debug2->Fill(tr_o.th_x);
  debug3->Fill(tr_o.th_y);
#endif
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::PrintReport(const RPRecoElasticEvent &elasticReco)
{
  printf("number of roads = %lu\n", elasticReco.roads.size());
  for (unsigned int i = 0; i < elasticReco.roads.size(); i++) {
    printf("\t[%i] at th_x = %8.1f, th_y = %8.1f (mu rad) : ", i, elasticReco.roads[i].centerX()*1E6, elasticReco.roads[i].centerY()*1E6);
    for (unsigned int j = 0; j < elasticReco.roads[i].members.size(); j++)
    printf("%i, ", elasticReco.roads[i].members[j]);
    printf("\n");
  }

  printf("preferred road = %i\n", elasticReco.preferredRoad);

  if (elasticReco.preferredRoad >= 0) {
    printf("---------------------------------------------------------------------------------------------------------\n");
    printf("   fit | projection | theta (rad) | theta error |  vertex (m) |  vertex err |     ndf |  S2_min/ndf |\n");
    printf("---------------------------------------------------------------------------------------------------------\n");
    RPRecoElasticEvent::fit_type* f = (RPRecoElasticEvent::fit_type *) (&elasticReco.leftFit);
    printf("  left |      x |%12.1f |%12.1f |%12.2E |%12.2E |%12i |%12.2E |\n", f->th_x*1E6, f->si_th_x*1E6, f->x, f->si_x, f->ndf_x, f->s2minPerDf_x());
    printf("          y |%12.1f |%12.1f |%12.2E |%12.2E |%12i |%12.2E |\n", f->th_y*1E6, f->si_th_y*1E6, f->y, f->si_y, f->ndf_y, f->s2minPerDf_y());
    printf("---------------------------------------------------------------------------------------------------------\n");
    f = (RPRecoElasticEvent::fit_type *) (&elasticReco.rightFit);
    printf(" right |      x |%12.1f |%12.1f |%12.2E |%12.2E |%12i |%12.2E |\n", f->th_x*1E6, f->si_th_x*1E6, f->x, f->si_x, f->ndf_x, f->s2minPerDf_x());
    printf("          y |%12.1f |%12.1f |%12.2E |%12.2E |%12i |%12.2E |\n", f->th_y*1E6, f->si_th_y*1E6, f->y, f->si_y, f->ndf_y, f->s2minPerDf_y());
    printf("---------------------------------------------------------------------------------------------------------\n");
    f = (RPRecoElasticEvent::fit_type *) (&elasticReco.globalFit);
    printf("global |      x |%12.1f |%12.1f |%12.2E |%12.2E |%12i |%12.2E |\n", f->th_x*1E6, f->si_th_x*1E6, f->x, f->si_x, f->ndf_x, f->s2minPerDf_x());
    printf("          y |%12.1f |%12.1f |%12.2E |%12.2E |%12i |%12.2E |\n", f->th_y*1E6, f->si_th_y*1E6, f->y, f->si_y, f->ndf_y, f->s2minPerDf_y());
    printf("---------------------------------------------------------------------------------------------------------\n");
  }
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::finalize()
{
  // print summary
  map<signed int, string> pNames;
  pNames[2212] = "p+";
  pNames[2112] = "n";
  pNames[211] = "pi+";
  pNames[-211] = "pi-";
  pNames[11] = "e-";
  pNames[-11] = "e+";
  pNames[22] = "gamma";

  printf("\n\nSUMMARY:\n");
  printf("------------------------------------------------------------------------------------------------------------------------------------\n");
  printf("  particle |   total   | triggered | OK recon. |   no road | no good r | one-sided | incon. x* | incon. y* | inc. th_x | inc. th_y |\n");
  printf("------------------------------------------------------------------------------------------------------------------------------------\n");
  for (map<signed int, RecoCounter>::iterator it = counters.begin(); it != counters.end(); ++it) {
    RecoCounter& c = it->second;
    printf("%10s |%10i |%10i |%10i |%10i |%10i |%10i |%10i |%10i |%10i |%10i |\n", pNames[it->first].c_str(), c.tot, c.trig, c.rOK, c.rNoRoad,
      c.rNoGoodRoad, c.rOnlyOneArm, c.rBadVertexX, c.rBadVertexY, c.rBadAngleX, c.rBadAngleY);
  }
  printf("------------------------------------------------------------------------------------------------------------------------------------\n");
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::writeHistogramsToFile()
{
  if (outputFile.empty())
    return;

  TFile *of = TFile::Open(outputFile.c_str(), "recreate");
  if (!of || !of->IsWritable()) {
    std::cout << "Output file not opened correctly!!" << std::endl;
    return;
  }

  TDirectory *d1;

  // write vertex plots
  d1 = of->mkdir("vertex");
  gDirectory = d1->mkdir("reco global vs. smeared");
  vtx_x_rs_diff->Write();
  vtx_y_rs_diff->Write();
  vtx_x_rs_corr->Write();
  vtx_y_rs_corr->Write();
  gDirectory = d1->mkdir("reco global vs. original");
  vtx_x_ro_diff->Write();
  vtx_y_ro_diff->Write();
  vtx_x_ro_corr->Write();
  vtx_y_ro_corr->Write();

  // write theta plots
  d1 = of->mkdir("theta");
  gDirectory = d1->mkdir("smear right vs. original");
  th_x_sro_diff->Write();
  th_y_sro_diff->Write();
  gDirectory = d1->mkdir("smear left vs. original");
  th_x_slo_diff->Write();
  th_y_slo_diff->Write();
  gDirectory = d1->mkdir("reco right vs. original");
  th_rro_diff->Write();
  th_x_rro_diff->Write();
  th_y_rro_diff->Write();
  gDirectory = d1->mkdir("reco left vs. original");
  th_rlo_diff->Write();
  th_x_rlo_diff->Write();
  th_y_rlo_diff->Write();
  gDirectory = d1->mkdir("reco right vs. smeared right");
  th_rrsr_diff->Write();
  th_x_rrsr_diff->Write();
  th_y_rrsr_diff->Write();
  gDirectory = d1->mkdir("reco left vs. smeared left");
  th_rlsl_diff->Write();
  th_x_rlsl_diff->Write();
  th_y_rlsl_diff->Write();
  gDirectory = d1->mkdir("reco global vs. original");
  th_rgo_diff->Write();
  th_x_rgo_diff->Write();
  th_y_rgo_diff->Write();

  th_rgo_corr->Write();
  th_x_rgo_corr->Write();
  th_y_rgo_corr->Write();

  th_rgo_res->Write();
  th_x_rgo_res->Write();
  th_y_rgo_res->Write();

  // write t plots
  d1 = of->mkdir("t");
  gDirectory = d1->mkdir("reco left vs. smeared left");
  t_rlsl_diff->Write();
  t_x_rlsl_diff->Write();
  t_y_rlsl_diff->Write();
  gDirectory = d1->mkdir("reco right vs. smeared right");
  t_rrsr_diff->Write();
  t_x_rrsr_diff->Write();
  t_y_rrsr_diff->Write();
  gDirectory = d1->mkdir("reco global vs. original");
  t_rgo_diff->Write();
  t_x_rgo_diff->Write();
  t_y_rgo_diff->Write();

  printf("\n\nFITTING resolution curves:\n");
  t_rgo_res->Write();
  t_x_rgo_res->Write();
  t_y_rgo_res->Write();

  // write distribution plots
  d1 = of->mkdir("distributions");
  gDirectory = d1->mkdir("reco global");
  t_rg_dist->Write();
  th_rg_dist->Write();
  phi_rg_dist->Write();
  gDirectory = d1->mkdir("smeared left");
  th_sl_dist->Write();
  t_sl_dist->Write();
  phi_sl_dist->Write();
  gDirectory = d1->mkdir("smeared right");
  t_sr_dist->Write();
  th_sr_dist->Write();
  phi_sr_dist->Write();
  gDirectory = d1->mkdir("original");
  phi_o_dist->Write();
  t_o_dist->Write();
  th_o_dist->Write();

  // init left-right diffrence plots
  d1 = of->mkdir("right-left differences");
  gDirectory = d1->mkdir("reco");
  th_x_r_rldiff->Write();
  th_y_r_rldiff->Write();
  vtx_x_r_rldiff->Write();
  vtx_y_r_rldiff->Write();
  gDirectory = d1->mkdir("smeared");
  th_x_s_rldiff->Write();
  th_y_s_rldiff->Write();

  // statistics
  d1 = of->mkdir("statistics");
  gDirectory = d1->mkdir("left");
  aRPCount_left->Write();
  fStat_left.Write();
  gDirectory = d1->mkdir("right");
  aRPCount_right->Write();
  fStat_right.Write();
  gDirectory = d1->mkdir("global");
  aRPCount_global->Write();
  fStat_global.Write();
  
  // road-size histograms
  gDirectory = of->mkdir("road-size");
  rs_x->Write();
  rs_y->Write();
  
#ifdef _ER_DEBUG_
  gDirectory = of->mkdir("debug");
  debug1->Write();
  debug2->Write();
  debug3->Write();
#endif

  delete of;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::ExportAllHistograms()
{
  ExportHistogram(vtx_x_rs_diff);
  ExportHistogram(vtx_y_rs_diff);
  ExportHistogram(vtx_x_rs_corr);
  ExportHistogram(vtx_y_rs_corr);
  ExportHistogram(vtx_x_ro_diff);
  ExportHistogram(vtx_y_ro_diff);
  ExportHistogram(vtx_x_ro_corr);
  ExportHistogram(vtx_y_ro_corr);

  ExportHistogram(th_rlsl_diff);
  ExportHistogram(th_rlo_diff);
  ExportHistogram(th_x_rlo_diff);
  ExportHistogram(th_y_rlo_diff);
  ExportHistogram(th_x_rlsl_diff);
  ExportHistogram(th_y_rlsl_diff);
  ExportHistogram(th_rrsr_diff);
  ExportHistogram(th_x_rrsr_diff);
  ExportHistogram(th_y_rrsr_diff);
  ExportHistogram(th_rro_diff);
  ExportHistogram(th_x_rro_diff);
  ExportHistogram(th_y_rro_diff);
  ExportHistogram(th_x_sro_diff);
  ExportHistogram(th_y_sro_diff);
  ExportHistogram(th_x_slo_diff);
  ExportHistogram(th_y_slo_diff);
  ExportHistogram(th_rgo_diff);
  ExportHistogram(th_x_rgo_diff);
  ExportHistogram(th_y_rgo_diff);
  ExportHistogram(th_rgo_corr);
  ExportHistogram(th_x_rgo_corr);
  ExportHistogram(th_y_rgo_corr);

  ExportHistogram(t_rlsl_diff);
  ExportHistogram(t_x_rlsl_diff);
  ExportHistogram(t_y_rlsl_diff);
  ExportHistogram(t_rrsr_diff);
  ExportHistogram(t_x_rrsr_diff);
  ExportHistogram(t_y_rrsr_diff);
  ExportHistogram(t_rgo_diff);
  ExportHistogram(t_x_rgo_diff);
  ExportHistogram(t_y_rgo_diff);

  ExportHistogram(t_rg_dist);
  ExportHistogram(th_rg_dist);
  ExportHistogram(phi_rg_dist);
  ExportHistogram(th_sl_dist);
  ExportHistogram(t_sl_dist);
  ExportHistogram(phi_sl_dist);
  ExportHistogram(t_sr_dist);
  ExportHistogram(th_sr_dist);
  ExportHistogram(phi_sr_dist);
  ExportHistogram(phi_o_dist);
  ExportHistogram(t_o_dist);
  ExportHistogram(th_o_dist);

  ExportHistogram(th_x_r_rldiff);
  ExportHistogram(th_y_r_rldiff);
  ExportHistogram(vtx_x_r_rldiff);
  ExportHistogram(vtx_y_r_rldiff);
  ExportHistogram(th_x_s_rldiff);
  ExportHistogram(th_y_s_rldiff);

  ExportHistogram(aRPCount_left);
  ExportHistogram(aRPCount_right);
  ExportHistogram(aRPCount_global);
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::ExportHistogram(TH1D *h)
{
  TH1::AddDirectory(kFALSE);
  TCanvas *c1 = new TCanvas(h->GetName(),h->GetName());
  c1->cd();
  h->Draw();
  elValHists.push_back(c1);
  //  delete c1;
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::ExportHistogram(TProfile *h)
{
  TCanvas *c2 = new TCanvas(h->GetName(),h->GetName());
  c2->cd();
  h->Draw();
  elValHists.push_back(c2);
  //  delete c2;
}

//----------------------------------------------------------------------------------------------------

void ElasticRecoValLibrary::ExportHistogram(TGraph *h)
{
  TCanvas *c3 = new TCanvas(h->GetName(),h->GetName());
  c3->cd();
  h->Draw("P");
  elValHists.push_back(c3);
  //  delete c3;
}

