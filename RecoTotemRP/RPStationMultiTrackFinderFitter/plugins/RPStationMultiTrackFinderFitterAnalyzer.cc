/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.sawicki@cern.ch, jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/RPStationMultiTrackFinderFitterAnalyzer.h"

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMath.h"

namespace RPStationMultiTrackFinderFitter {

RPStationMultiTrackFinderFitterAnalyzer::RPStationMultiTrackFinderFitterAnalyzer(const edm::ParameterSet &ps) :
    geometry(ps),
    verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),

    simTracksProducer(ps.getParameter<string>("simTracksProducer")),
    recoTracksProducer(ps.getParameter<string>("recoTracksProducer")),
    oneRPPatternsProducer(ps.getParameter<string>("oneRPPatternsProducer")),

    outputFileName(ps.getParameter<string>("outputFileName")),
    rootFileName(ps.getParameter<string>("rootFileName")),

    missing_max(ps.getParameter<unsigned int>("missing_max")),
    fake_max(ps.getParameter<unsigned int>("fake_max")),
    si_de_ax(ps.getParameter<double>("si_de_ax")),
    si_de_ay(ps.getParameter<double>("si_de_ay")),
    si_de_bx(ps.getParameter<double>("si_de_bx")),
    si_de_by(ps.getParameter<double>("si_de_by")),
    acceptance(ps.getParameter<double>("acceptance")),

    f_out_enable(outputFileName.size() > 0),
    tf_out_enable(rootFileName.size() > 0),

    stats({0,0,0,0,0})
{
    if (f_out_enable)
      f_out = fopen(outputFileName.c_str(), "w");

    if (tf_out_enable)
    {
        missing_fake = new TH2D("missing_fake", "missing_fake",
                missing_max, 0, missing_max,
                fake_max, 0, fake_max);

        missing_fake->GetXaxis()->SetTitle("missing");
        missing_fake->GetYaxis()->SetTitle("false");

        h_min_pat_dist_all = new TH1D("h_min_pat_dist_all", ";min distance between U/V patterns   (mm)", 100, 0., 6.6);
        h_min_pat_dist_all->SetLineColor(1);

        h_min_pat_dist_fail = new TH1D("h_min_pat_dist_fail", ";min distance between U/V patterns   (mm)", 100, 0, 6.6);
        h_min_pat_dist_fail->SetLineColor(2);

        h_min_pat_dist_fake_comb = new TH1D("h_min_pat_dist_fake_comb", ";min distance between U/V patterns   (mm)", 100, 0, 6.6);
        h_min_pat_dist_fake_comb->SetLineColor(4);

        // TODO: units ??
        track_distances = new TH1D("track_distances", ";track distance squared", 100, 0, 0);

        dist_x = new TH2D("dist_x", "dist_x", 100, 0, 0, 100, 0, 0);
        dist_y = new TH2D("dist_y", "dist_y", 100, 0, 0, 100, 0, 0);
        fail_dist_x = new TH2D("fail_dist_x", "fail_dist_x", 100, 0, 0, 100, 0, 0);
        fail_dist_y = new TH2D("fail_dist_y", "fail_dist_y", 100, 0, 0, 100, 0, 0);

        dist_u = new TH2D("dist_u", "dist_u", 100, 0, 0, 100, 0, 0);
        dist_v = new TH2D("dist_v", "dist_v", 100, 0, 0, 100, 0, 0);
        fail_dist_u = new TH2D("fail_dist_u", "fail_dist_u", 100, 0, 0, 100, 0, 0);
        fail_dist_v = new TH2D("fail_dist_v", "fail_dist_v", 100, 0, 0, 100, 0, 0);

        h_num_reas_fits = new TH1D("h_num_reas_fits", ";reasonable fits;fraction", 11, -0.5, 10.5);
    }
}

//----------------------------------------------------------------------------------------------------

RPStationMultiTrackFinderFitterAnalyzer::~RPStationMultiTrackFinderFitterAnalyzer()
{
    // finalize stats file
    WriteSummary();

    if (f_out_enable)
      fclose(f_out);

    if (tf_out_enable)
    {
        TFile *tf_out = new TFile(rootFileName.c_str(), "recreate");

        // finalize root file
        Double_t norm = missing_fake->GetEntries();
        if (norm != 0.)
          missing_fake->Scale(1/norm);
        missing_fake->Write();

        //norm = h_min_pat_dist_all->GetEntries();
        //if (norm != 0.) h_min_pat_dist_all->Scale(1/norm);
        h_min_pat_dist_all->Write();

        //norm = h_min_pat_dist_fail->GetEntries();
        //if (norm != 0.) h_min_pat_dist_fail->Scale(1/norm);
        h_min_pat_dist_fail->Write();

        h_min_pat_dist_fake_comb->Write();

        norm = track_distances->GetEntries();
        if (norm != 0.)
          track_distances->Scale(1./norm);
        track_distances->Write();

        dist_x->Write();
        dist_y->Write();
        fail_dist_x->Write();
        fail_dist_y->Write();

        dist_u->Write();
        dist_v->Write();
        fail_dist_u->Write();
        fail_dist_v->Write();

        h_num_reas_fits->Scale(1./h_num_reas_fits->GetEntries());
        h_num_reas_fits->Write();

        delete tf_out;
    }
}

//----------------------------------------------------------------------------------------------------

void RPStationMultiTrackFinderFitterAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &es)
{
    using namespace edm;

    if (verbosity > 0)
    {
        printf("\n>> RPStationMultiTrackFinderFitterAnalyzer::analyze(%u:%u)\n",
                event.id().run(), event.id().event());
    }

    // get geometry
    if (watcherRealGeometry.check(es))
    {
        ESHandle<TotemRPGeometry> newGeometry;
        es.get<RealGeometryRecord>().get(newGeometry);
        geometry.setGeometry(&(*newGeometry));
    }

    edm::Handle<RPStationTrackFitCollection> station_simu;
    event.getByLabel(simTracksProducer, station_simu);

    edm::Handle<RPStationTrackFitCollection> station_reco;
    event.getByLabel(recoTracksProducer, station_reco);

    edm::Handle<RPRecognizedPatternsCollection> patterns;
    event.getByLabel(oneRPPatternsProducer, patterns);

    const vector<RPStationTrackFit> *tracks_simu = NULL, *tracks_reco = NULL;

    unsigned int valid = 0;

    if (verbosity > 0)
        printf("* SIMULATED TRACKS\n");

    for (RPStationTrackFitCollection::const_iterator ait = station_simu->begin(); ait != station_simu->end(); ++ait)
    {
        tracks_simu = &ait->second;
        for (vector<RPStationTrackFit>::const_iterator fit = ait->second.begin(); fit != ait->second.end(); ++fit)
        {
            const RPStationTrackFit& tf = *fit;
            if (tf.valid)
              valid++;

            if (verbosity > 0)
            {
                printf("\tax=%+.1f urad, ay=%+.1f urad, bx=%+.3f mm, by=%+.3f mm, z0=%+.3f mm, reconstructable=%u\n",
                        tf.ax*1E6, tf.ay*1E6, tf.bx, tf.by, tf.z0, tf.valid);
            }
        }
    }

    if (verbosity > 0)
        printf("reconstructable: %u\n", valid);

    if (verbosity > 0)
        printf("* RECONSTRUCTED TRACKS\n");

    for (RPStationTrackFitCollection::const_iterator ait = station_reco->begin(); ait != station_reco->end(); ++ait)
    {
        tracks_reco = &ait->second;
        if (verbosity > 0)
        {
            for (vector<RPStationTrackFit>::const_iterator fit = ait->second.begin(); fit != ait->second.end(); ++fit)
            {
                const RPStationTrackFit& tf = *fit;
                printf("\tax=%+.1f urad, ay=%+.1f urad, bx=%+.3f mm, by=%+.3f mm, z0=%+.3f\n",
                        tf.ax*1E6, tf.ay*1E6, tf.bx, tf.by, tf.z0);
            }
        }
    }

    vector<RPStationTrackFit> dummy;
    if (tracks_reco == 0)
        tracks_reco = &dummy;

    if (tracks_simu == 0)
        tracks_simu = &dummy;

    missing_fake_pair res = AnalyzeTracks(event, *tracks_simu, *tracks_reco);

    unsigned int n_reas_fits = AnalyzeReasonableInputFits(*patterns);
    bool fake_comb = (n_reas_fits > tracks_simu->size());

    AnalyzePatterns(*patterns, res, fake_comb);

    AnalyzeSimuTracks(*tracks_simu, res);

    if ((res.missing + res.fake == 0) && fake_comb)
      printf("THIS\n");
}

//----------------------------------------------------------------------------------------------------

RPStationMultiTrackFinderFitterAnalyzer::missing_fake_pair
RPStationMultiTrackFinderFitterAnalyzer::AnalyzeTracks(const edm::Event &event,
        const vector<RPStationTrackFit> &tracks_simu,
        const vector<RPStationTrackFit> &tracks_reco)
{
    unsigned int missing = 0,
                 fake = 0,
                 total_simu_valid = 0;

    set<unsigned int> found;

    for (unsigned int i = 0; i < tracks_simu.size() && tracks_simu[i].valid; i++)
    {
        total_simu_valid++;

        bool found_match = false;
        double dist_sq;

        for (unsigned int j = 0; j < tracks_reco.size() && !found_match; j++)
        {
            if (found.find(j) != found.end())
              continue;

            const RPStationTrackFit &t_s = tracks_simu[i],
                                    &t_r = tracks_reco[j];

            double de_ax = (t_r.ax - t_s.ax)*1E6;
            double de_bx = (t_r.bx - t_s.bx - t_s.ax*(t_r.z0 - t_s.z0));
            double de_ay = (t_r.ay - t_s.ay)*1E6;
            double de_by = (t_r.by - t_s.by - t_s.ay*(t_r.z0 - t_s.z0));

            double dist_de_ax = de_ax / si_de_ax;
            double dist_de_ay = de_ay / si_de_ay;
            double dist_de_bx = de_bx / si_de_bx;
            double dist_de_by = de_by / si_de_by;

            // normalisation ??
            dist_sq = dist_de_ax*dist_de_ax + dist_de_ay*dist_de_ay + dist_de_bx*dist_de_bx + dist_de_by*dist_de_by;

            if (dist_sq < acceptance)
            {
                found.insert(j);
                found_match = true;
            }
        }

        if (!found_match)
        {
            missing++;
        } else {
          if (tf_out_enable)
            track_distances->Fill(dist_sq);
        }
    }

    if (total_simu_valid == 0)
        return {0, 0};

    fake = tracks_reco.size() + missing - total_simu_valid;

    if (f_out_enable)
    {
        unsigned int run = event.id().run(),
                     ev  = event.id().event();

        fprintf(f_out, "%u:%u %u/%u\n", run, ev, missing, fake);
    }

    stats.total++;
 
    if (fake + missing > 0)
        stats.failed++;

    if (fake > 0)
        stats.fake++;

    if (missing > 0)
        stats.missing++;

    if (fake * missing > 0)
        stats.fakemissing++;

    if (tf_out_enable)
      missing_fake->Fill(missing, fake);

    return {missing, fake};
}

//----------------------------------------------------------------------------------------------------

void RPStationMultiTrackFinderFitterAnalyzer::AnalyzePatterns(
        const RPRecognizedPatternsCollection &patterns,
        const missing_fake_pair &miss_fake,
        bool fake_comb)
{
    if
      (!tf_out_enable) return;

    double minDist = 1e100;

    // loop over all RPs
    for (RPRecognizedPatternsCollection::const_iterator rpIt = patterns.begin(); rpIt != patterns.end(); rpIt++)
    {
        // calculate distances between any two U patterns
        if (rpIt->second.uLines.size() > 1)
        {
            for (vector<RPRecognizedPatterns::Line>::const_iterator lineIt = rpIt->second.uLines.begin();
                    lineIt != rpIt->second.uLines.end(); lineIt++)
            {
                for (vector<RPRecognizedPatterns::Line>::const_iterator lineIt2 = lineIt + 1;
                        lineIt2 != rpIt->second.uLines.end(); lineIt2++)
                {
                    double dist = fabs(lineIt->b - lineIt2->b);
                    if (dist < minDist)
                        minDist = dist;
                }
            }
        }

        // calculate distances between any two V patterns
        if (rpIt->second.vLines.size() > 1)
        {
            for (vector<RPRecognizedPatterns::Line>::const_iterator lineIt = rpIt->second.vLines.begin();
                    lineIt != rpIt->second.vLines.end(); lineIt++)
            {
                for (vector<RPRecognizedPatterns::Line>::const_iterator lineIt2 = lineIt + 1;
                        lineIt2 != rpIt->second.vLines.end(); lineIt2++)
                {
                    double dist = fabs(lineIt->b - lineIt2->b);
                    if (dist < minDist)
                        minDist = dist;
                }
            }
        }
    }

    h_min_pat_dist_all->Fill(minDist);

    if (miss_fake.missing + miss_fake.fake > 0)
      h_min_pat_dist_fail->Fill(minDist);

    if (fake_comb)
      h_min_pat_dist_fake_comb->Fill(minDist);
}

//----------------------------------------------------------------------------------------------------

void RPStationMultiTrackFinderFitterAnalyzer::AnalyzeSimuTracks(
        const vector<RPStationTrackFit> &tracks_simu,
        const missing_fake_pair &miss_fake)
{
    if (!tf_out_enable)
      return;

    for (vector<RPStationTrackFit>::const_iterator track1It = tracks_simu.begin();
            track1It != tracks_simu.end(); track1It++)
    {
        for (vector<RPStationTrackFit>::const_iterator track2It = track1It + 1;
                track2It != tracks_simu.end(); track2It++)
        {
            const RPStationTrackFit &t_r = *track1It,
                                    &t_s = *track2It;

            //double si_de_ax = 3.,
                   //si_de_ay = 3.,
                   //si_de_bx = 0.02,
                   //si_de_by = 0.02;

            double de_ax = (t_r.ax - t_s.ax)*1E6;
            double de_bx = (t_r.bx - t_s.bx - t_s.ax*(t_r.z0 - t_s.z0));
            double de_ay = (t_r.ay - t_s.ay)*1E6;
            double de_by = (t_r.by - t_s.by - t_s.ay*(t_r.z0 - t_s.z0));

            double alpha = 45. / 180. * M_PI;
            double d_x = cos(alpha);
            double d_y = sin(alpha);

            double de_au = de_ax*d_x + de_ay*d_y,
                   de_bu = de_bx*d_x + de_by*d_y,
                   de_av = de_ax*d_y - de_ay*d_x,
                   de_bv = de_bx*d_y - de_by*d_x;

            //double dist_de_ax = de_ax / si_de_ax;
            //double dist_de_ay = de_ay / si_de_ay;
            //double dist_de_bx = de_bx / si_de_bx;
            //double dist_de_by = de_by / si_de_by;

            //double dist_sq = dist_de_ax*dist_de_ax + dist_de_ay*dist_de_ay + dist_de_bx*dist_de_bx + dist_de_by*dist_de_by;

            if (miss_fake.fake == 0 && miss_fake.missing == 2)
            {
                fail_dist_x->Fill(fabs(de_ax), fabs(de_bx));
                fail_dist_y->Fill(fabs(de_ay), fabs(de_by));
                fail_dist_u->Fill(fabs(de_au), fabs(de_bu));
                fail_dist_v->Fill(fabs(de_av), fabs(de_bv));
            } else if (fabs(de_av) < 110 && fabs(de_bv) < 1.4 &&
                    fabs(de_au) < 60 && fabs(de_bu) < 3)
            {
            //} else {
                dist_x->Fill(fabs(de_ax), fabs(de_bx));
                dist_y->Fill(fabs(de_ay), fabs(de_by));
                dist_u->Fill(fabs(de_au), fabs(de_bu));
                dist_v->Fill(fabs(de_av), fabs(de_bv));
            }
        }
    }

}

//----------------------------------------------------------------------------------------------------

void RPStationMultiTrackFinderFitterAnalyzer::WriteSummary()
{
    if (!f_out_enable)
      return;

    fprintf(f_out, "\nSUMMARY:\n");
    fprintf(f_out, "total: %u\n", stats.total);
    fprintf(f_out, "failed: %u, i.e. %.2f%%\n", stats.failed, 100.0*stats.failed/stats.total);
    fprintf(f_out, "fake: %u, i.e. %.2f%%\n", stats.fake, 100.0*stats.fake/stats.total);
    fprintf(f_out, "missing: %u, i.e. %.2f%%\n", stats.missing, 100.0*stats.missing/stats.total);
    fprintf(f_out, "both: %u, i.e. %.2f%%\n", stats.fakemissing, 100.0*stats.fakemissing/stats.total);
}

//----------------------------------------------------------------------------------------------------
  
unsigned int RPStationMultiTrackFinderFitterAnalyzer::AnalyzeReasonableInputFits(const RPRecognizedPatternsCollection &patterns)
{
  // fit quality limit to be considered as reasonable
  double sig_limit = 3.;  // sigmas

  // prepare indeces
  vector<unsigned int> indeces;
  vector<unsigned int> indeces_lim;
  vector<unsigned int> rp_ids;
  for (RPRecognizedPatternsCollection::const_iterator rpit = patterns.begin(); rpit != patterns.end(); ++rpit)
  {
    rp_ids.push_back(rpit->first);

    indeces.push_back(0);
    indeces.push_back(0);

    indeces_lim.push_back(rpit->second.uLines.size());
    indeces_lim.push_back(rpit->second.vLines.size());
  }

  // build fit matrix
  unsigned int n_RP = patterns.size();
  TMatrixD A(2*n_RP, 4);
  unsigned int rpi = 0;
  for (RPRecognizedPatternsCollection::const_iterator rpit = patterns.begin(); rpit != patterns.end(); ++rpit, ++rpi)
  {
    const RPInfo &g = geometry.getRPInfo(rpit->first);
    double z_eff = g.cz - 213E3;
    A(2*rpi + 0, 0) = z_eff * g.u_dx; A(2*rpi + 0, 1) = g.u_dx; A(2*rpi + 0, 2) = z_eff * g.u_dy; A(2*rpi + 0, 3) = g.u_dy;
    A(2*rpi + 1, 0) = z_eff * g.v_dx; A(2*rpi + 1, 1) = g.v_dx; A(2*rpi + 1, 2) = z_eff * g.v_dy; A(2*rpi + 1, 3) = g.v_dy;
  }

  //printf(">> A\n");
  //A.Print();
  
  TMatrixD AT(TMatrixD::kTransposed, A);
  TMatrixD ATAi(TMatrixD::kInverted, AT * A);
  TMatrixD ATAi_AT = ATAi * AT;

  TVectorD M(2*n_RP);

  // loop over all possible index configurations
  bool stop = false;
  unsigned int N_reas_fits = 0;
  while (!stop)
  {
    // print indeces
    if (verbosity > 5)
    {
      printf("* ");
      for (unsigned int ii = 0; ii < indeces.size(); ++ii)
        printf("%u, ", indeces[ii]);
      printf("\n");
    }

    // build data vector
    rpi = 0;
    for (RPRecognizedPatternsCollection::const_iterator rpit = patterns.begin(); rpit != patterns.end(); ++rpit, ++rpi)
    {
      unsigned int idx_u = indeces[2*rpi + 0];
      unsigned int idx_v = indeces[2*rpi + 1];

      M(2*rpi + 0) = rpit->second.uLines[idx_u].b;
      M(2*rpi + 1) = rpit->second.vLines[idx_v].b;
    }

    // TODO: remove
    //printf(">> M\n");
    //M.Print();

    // calculate chi square
    TVectorD par = ATAi_AT * M;
    TVectorD diff = M - A * par;

    // TODO: remove
    //printf(">> par\n");
    //par.Print();
    //printf(">> diff\n");
    //diff.Print();

    double chi_sq = 0;
    for (unsigned int i = 0; i < 2*n_RP; i++)
      chi_sq += diff(i) * diff(i);

    double meas_unc = 19E-3;  // mm
    chi_sq /= meas_unc * meas_unc;

	unsigned int ndf = 2*n_RP - 4;
	double prob = TMath::Prob(chi_sq, ndf);
	double sig = (prob <= 0.) ? 999. : sqrt(2.) * TMath::ErfcInverse(prob);

    if (verbosity > 5)
      printf("----> chi^2 = %E, prob = %E, sig = %.1f\n", chi_sq, prob, sig);

    if (sig < sig_limit)
      N_reas_fits++;

    // advance indeces
    bool advanced = false;
    for (unsigned int ii = 0; ii < indeces.size(); ++ii)
    {
      if (indeces[ii] + 1 < indeces_lim[ii])
      {
        indeces[ii]++;
        advanced = true;
        break;
      } else {
        indeces[ii] = 0;
      }
    }
    stop = !advanced;
  }

  printf("\n* reasonable hit combinations: %u\n", N_reas_fits);

  // fill results
  h_num_reas_fits->Fill(N_reas_fits);

  return N_reas_fits;
}

// define as plug-in
DEFINE_FWK_MODULE(RPStationMultiTrackFinderFitterAnalyzer);

} // namespace RPStationMultiTrackFinderFitter
