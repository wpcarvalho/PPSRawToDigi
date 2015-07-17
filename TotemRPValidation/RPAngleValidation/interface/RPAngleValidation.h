#ifndef TotemRPValidation_RPAngleValidation_RPAngleValidation_h
#define TotemRPValidation_RPAngleValidation_RPAngleValidation_h

#include <boost/shared_ptr.hpp>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "TH2D.h"

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}

struct PrimaryProton
{
  HepMC::FourVector vertex;
  HepMC::FourVector momentum;
  double thetaX, thetaY;
  bool found;
};

struct TrackHit
{
	double x, y, z;
};

class RPAngleValidation : public edm::EDAnalyzer
{
  public:
    explicit RPAngleValidation(const edm::ParameterSet&);
    ~RPAngleValidation();
    
  private:
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    enum { rp_150_l, rp_220_l, rp_150_r, rp_220_r, NSTATIONS };
    enum { near_top, near_bottom, near_horiz, far_horiz, far_top, far_bottom, NPOTS };
    typedef std::map<RPId, TrackHit> rec_tracks_collection;
    typedef std::vector<int> station_rp_ids_type;

    bool FindProtons(const edm::Handle<edm::HepMCProduct> &HepMCEvt, int smeared);

    int SelectHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll, const station_rp_ids_type& st_ids);

    int calcTrackTheta(std::vector<double> &theta, const rec_tracks_collection &coll, int arm);
    std::string StIdToName(int st_id);

    edm::ParameterSet conf_;
    BeamOpticsParams BOPar_;
    struct PrimaryProton prot_[2][2]; /// left, right : primary, secondary
    std::string OriginalHepMCProductLabel_;
    std::string OriginalHepMCModuleName_;
    std::string SmearedHepMCProductLabel_;
    std::string SmearedHepMCModuleName_;
    std::string hist_file_name_;
    station_rp_ids_type station_ids_[NSTATIONS]; /// station
    std::auto_ptr<TH2D> angle_dists_[2][NSTATIONS]; /// X, Y : station
    std::auto_ptr<TH1D> angle_dists_1d_[2][2][NSTATIONS]; /// IP, RP : X, Y : station
    int verbosity_;
    edm::InputTag rpFittedTrackCollectionLabel;
};


#endif
