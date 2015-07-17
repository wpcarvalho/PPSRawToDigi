#ifndef Validation_RecoTrackRP_RecoTrackRPVal_h
#define Validation_RecoTrackRP_RecoTrackRPVal_h

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "HepMC/SimpleVector.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDigCluster.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RP2DHit.h"
#include "RecoTotemRP/RPRomanPotResolutionService/interface/RPFitResolution.h"

#include "TotemRPValidation/BaseValidationClasses/interface/BaseCollectionManager.h"
#include "TotemRPValidation/ParamMADRefTransport/interface/ParamMADRefTransport.h"
#include "TotemRPValidation/RecoTrackRP/interface/DetInfo.h"
#include "TotemRPValidation/RecoTrackRP/interface/RPTrackInfo.h"

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


class RecoTrackRPVal : public edm::EDAnalyzer
{
  public:
    explicit RecoTrackRPVal(const edm::ParameterSet&);
    ~RecoTrackRPVal();
    
  private:
    edm::InputTag rPFittedTrackCollectionLabel;
    edm::InputTag clusterSetLabel;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    
    enum st_id { rp_150_l, rp_220_l, rp_150_r, rp_220_r, NSTATIONS };
    enum pot_pos { near_top, near_bottom, near_horiz, far_horiz, far_top, far_bottom, NPOTS };

    typedef edm::DetSetVector<RPDigCluster> cluster_set;
    typedef std::map<RPId, RP2DHit> rec_tracks_collection;
    typedef std::vector<RPId> station_rp_ids_type;
    
    void InitBeamSmearingHistograms();
    void InitOverlapsHistograms();
    void InitHitHistograms();
    void InitHitDistHistograms();
    void InitAngleHistograms();
    bool FindProtons(const edm::Handle<edm::HepMCProduct> &HepMCEvt, int smeared);
    void WriteBeamSmearingHistograms(TFile *f);
    void WriteOverlapsHistograms(TFile *f);
    void WriteHitHistorgrams(TFile *f);
    void WriteHitDistHistograms(TFile*);
    void WriteAngleHistograms(TFile *f);
    void FillReferenceHistograms(RPId id, const PrimaryProton &proton);
    void FillResidualHistograms(const rec_tracks_collection &tracks, 
        const PrimaryProton &prim_proton, 
        const PrimaryProton &smeared_proton, bool proton_at_150_station);
    void FillDetectorHistograms(const RPFittedTrack &track, const cluster_set &clusters);
    void FillDetectorHistogramsForReconstructedTracks( 
          const RPFittedTrackCollection &tracks, const cluster_set &clusters);
    void FillEventSmearingInformation();
    void FillTrackPositions(edm::Handle< RPFittedTrackCollection > &input);
    void FillOverlappingAreas(const RPFittedTrackCollection &tracks);
    void FillHits(const RPFittedTrackCollection &tracks);
    void FillHitDists(const RPFittedTrackCollection &tracks);
    void FillAngleHistograms(const rec_tracks_collection &tracks, int sid);
    
    void setupCanvas();
    void DrawRPContour(RPId rpid);
    void WriteHistogramAsCanvas(TH2F *h, const std::vector<RPId> &rpids);
          
    int SelectHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll, const station_rp_ids_type& st_ids);
    int calcTrackTheta(std::vector<double> &theta, const rec_tracks_collection &coll, int arm);
    std::string StIdToName(int st_id);
    
    TotemRPGeometry Totem_RP_geometry_;
    BaseCollectionManager<RPTrackInfo, RPId, edm::ParameterSet> rp_histograms_;
    BaseCollectionManager<DetInfo, RPDetId, edm::ParameterSet> det_histograms_;
    
    std::auto_ptr<ParamMADRefTransport> parameterized_madx_transport_;
    struct PrimaryProton prot_[2][2]; /// left, right : primary, smeared
    
    std::string hist_file_name_;
    std::string OriginalHepMCProductLabel_;
    std::string OriginalHepMCModuleName_;
    std::string SmearedHepMCProductLabel_;
    std::string SmearedHepMCModuleName_;
    int verbosity_;
    int validTracksOnly;
    
    station_rp_ids_type station_ids_[NSTATIONS];
    std::auto_ptr<TH2D> angle_dists_[2][NSTATIONS]; /// X, Y : station
    std::auto_ptr<TH1D> angle_dists_1d_[2][2][NSTATIONS]; /// IP, RP : X, Y : station
    
    RPFitResolution resol_degrad_service_;
    BeamOpticsParams BOPar_;
    edm::ParameterSet conf_;
    
    std::auto_ptr<TH1D> vertex_x_smearing_dist_;
    std::auto_ptr<TH1D> vertex_y_smearing_dist_;
    std::auto_ptr<TH1D> vertex_z_smearing_dist_;
    
    std::auto_ptr<TH1D> px_smearing_dist_[2];
    std::auto_ptr<TH1D> py_smearing_dist_[2];
    std::auto_ptr<TH1D> pz_smearing_dist_[2];
    
    std::auto_ptr<TH1D> thx_smearing_dist_[2];
    std::auto_ptr<TH1D> thy_smearing_dist_[2];
    
    std::auto_ptr<TH1D> xi_smearing_dist_[2];
    
    //Overlapping tracks' histograms
    std::auto_ptr<TH2F> overlaps_RP_120_121_122_;
    std::auto_ptr<TH2F> overlaps_RP_100_101_102_;
    std::auto_ptr<TH2F> overlaps_RP_020_021_022_;
    std::auto_ptr<TH2F> overlaps_RP_000_001_002_;
    
    std::auto_ptr<TH2F> overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_;
    std::auto_ptr<TH2F> overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_;
    std::auto_ptr<TH2F> overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_;
    std::auto_ptr<TH2F> overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_;
    
    //Particle flux per station
    std::auto_ptr<TH2F> hits_rp_[NSTATIONS];
    std::map< unsigned int, TGraph *> hitDists;
};


#endif
