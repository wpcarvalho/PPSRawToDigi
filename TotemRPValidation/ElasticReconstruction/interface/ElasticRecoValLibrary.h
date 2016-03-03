/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors: 
*	Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: ElasticRecoValLibrary.h 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/


#ifndef _TotemRPValidationElasticReconstructionElasticRecoValLibrary_H_
#define _TotemRPValidationElasticReconstructionElasticRecoValLibrary_H_

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}

class TGraph;
class TGraphErrors;
class TF1;
class TCanvas;
class TProfile;

class BeamOpticsParams;

#include <vector>
#include <map>

#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecoElasticEvent.h"
#include "TH1D.h"


//----------------------------------------------------------------------------------------------------

/**
 *\brief Resolution graph.
 **/
class TResolution
{
  public:
    struct TMoments {
      double y0, y1, y2, y3, y4;
      double range;
      TH1F* hist;

      TMoments(double r) : y0(0.), y1(0.), y2(0.), y3(0.), y4(0.),
        range(r), hist(new TH1F("", "", 100, -r, +r)) {}

      void Fill(double y) {
        if (fabs(y) > range)
          return;
        y0 += 1.;
        y1 += y;
        y2 += y*y;
        y3 += y*y*y;
        y4 += y*y*y*y;
        hist->Fill(y);
      }
    };

  public:
    TResolution(const char* _name, const char* _title, unsigned int _N, double from, double to, double range);
    void Fill(double x, double y);
    void Write();

  private:
    std::string name, title;
    unsigned int N;
    double x0, dx;
    unsigned int entries, overfill, underfill;
    std::vector<TMoments> bins;
};


//----------------------------------------------------------------------------------------------------


/**
 *\brief The code for elastic reconstruction validation. Called from ElasticRecoVal analyzer.
**/
class ElasticRecoValLibrary
{
 public:
  ElasticRecoValLibrary(const edm::ParameterSet&);
  ~ElasticRecoValLibrary();

  void initialize(const edm::EventSetup&);
  void analyze(const edm::Event&, const edm::EventSetup&);
  void finalize();
  void writeHistogramsToFile();

  /// reconstruction status counters for a give particle type, see counters variable
  struct RecoCounter {
    unsigned int tot, trig;
    unsigned int rOK, rNoRoad, rNoGoodRoad, rOnlyOneArm, rBadVertexX, rBadVertexY, rBadAngleX, rBadAngleY;

    TH1D *th_x, *th_y;

    RecoCounter();
    void Update(bool triggered, const RPRecoElasticEvent &ee);
  };

  /// statistics properties of a fit
  struct FitStatistic {
    std::string name;
    TH1D *ndf_x, *ndf_y;  /// ndf = number of degrees of freedom
    TH1D *p_x, *p_y;      /// histogram of p values, see TMath::Prob
    std::map<unsigned int, TH1D *> p_ndf_x, p_ndf_y;
    FitStatistic(const std::string &_name);
    void Fill(unsigned int _ndf_x, unsigned int _ndf_y, double s2_min_x, double s2_min_y);
    void Write();
  };

 private:
  edm::InputTag rpDetTriggerSetLabel;
  edm::InputTag rpFittedTrackCollectionLabel;
  edm::InputTag rpRecoElasticEventLabel;


  unsigned char verbosity;														///< verbosity level
  std::string generatorLabel;                                                   ///< label of the HepMC product to be loaded as smeared event
  std::string originalLabel;                                                    ///< label of the HepMC product to be loaded as event before smaering
  std::string outputFile;														///< file to store output plots

  double lim_vtx_dx, lim_vtx_dy;
  double lim_th_dx, lim_th_dy;
  double lim_t_min, lim_t_max;

  std::map<signed int, RecoCounter> counters;									///< reconstruction status counters for a given particle type

  TH1D *aRPCount_left, *aRPCount_right, *aRPCount_global;						///< active RP count histograms

  TH1D *vtx_x_rs_diff, *vtx_y_rs_diff, *vtx_x_ro_diff, *vtx_y_ro_diff;			///< vertex difference plots: rs = reco vs. smeared, ro = reco vs. original
  TGraph *vtx_x_rs_corr, *vtx_y_rs_corr, *vtx_x_ro_corr, *vtx_y_ro_corr;		///< vertex correlation plots

  TH1D *th_rlsl_diff, *th_x_rlsl_diff, *th_y_rlsl_diff;							///< theta difference plots; rlsl = Left track from Reco vs. Left track from Smeared MC event
  TH1D *th_rrsr_diff, *th_x_rrsr_diff, *th_y_rrsr_diff;							
  TH1D *th_rlo_diff, *th_x_rlo_diff, *th_y_rlo_diff;
  TH1D *th_rro_diff, *th_x_rro_diff, *th_y_rro_diff;							
  TH1D *th_x_slo_diff, *th_y_slo_diff;
  TH1D *th_x_sro_diff, *th_y_sro_diff;							
  TH1D *th_rgo_diff, *th_x_rgo_diff, *th_y_rgo_diff;							///< rgo means Global Reco vs. Original MC event
  TGraph *th_rgo_corr, *th_x_rgo_corr, *th_y_rgo_corr;
  TResolution *th_rgo_res, *th_x_rgo_res, *th_y_rgo_res;						///< theta resolution plots

  TH1D *t_rlsl_diff, *t_x_rlsl_diff, *t_y_rlsl_diff;							///< t difference plots
  TH1D *t_rrsr_diff, *t_x_rrsr_diff, *t_y_rrsr_diff;
  TH1D *t_rgo_diff, *t_x_rgo_diff, *t_y_rgo_diff;
  TResolution *t_rgo_res, *t_x_rgo_res, *t_y_rgo_res;							///< t resolution plots
  TGraphErrors *t_rgo_resg, *t_x_rgo_resg, *t_y_rgo_resg; 

  TH1D *t_rg_dist, *t_sl_dist, *t_sr_dist, *t_o_dist;							///< t distribution plots
  TH1D *th_rg_dist, *th_sl_dist, *th_sr_dist, *th_o_dist;						///< theta distribution plots
  TH1D *phi_rg_dist, *phi_sl_dist, *phi_sr_dist, *phi_o_dist;					///< phi distribution plots

  TH1D *th_x_r_rldiff, *th_y_r_rldiff, *vtx_x_r_rldiff, *vtx_y_r_rldiff;		///< right-left differences; s = smeared MC event, r = reco
  TH1D *th_x_s_rldiff, *th_y_s_rldiff;

  FitStatistic fStat_left, fStat_right, fStat_global;                           ///< fit statistics

  /// road-size histograms
  TH1D *rs_x, *rs_y;

  static void PrintReport(const RPRecoElasticEvent&);

  void AddFunction(TH1D* h, TF1 *f);
  TGraphErrors* FitAndDrawErrProfile(TProfile *p, bool doFit = true, bool logScale = false);

#ifdef _ER_DEBUG_
  TGraph *debug1;
  TH1D *debug2, *debug3;
#endif

  /* interface to Leszek's validation kit */
 public:
  void ExportAllHistograms();
  std::vector<TCanvas*> getHistograms() { return elValHists; }

 private:
  std::vector<TCanvas*> elValHists;

  void ExportHistogram(TH1D *h);
  void ExportHistogram(TProfile *h);
  void ExportHistogram(TGraph *h);
};

#endif 

