#ifndef TotemRPValidation_RPReconstructedTracksValidation_RPReconstructedTracksValidation_h
#define TotemRPValidation_RPReconstructedTracksValidation_RPReconstructedTracksValidation_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TotemRPValidation/BaseValidationClasses/interface/BaseCollectionManager.h"
#include "TotemRPValidation/RPReconstructedTracksValidation/interface/RPTrackInfo.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "TotemRPValidation/ParamMADRefTransport/interface/ParamMADRefTransport.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProton.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RP2DHit.h"
#include "RecoTotemRP/RPRomanPotResolutionService/interface/RPFitResolution.h"
#include "TotemRPValidation/RPReconstructedTracksValidation/interface/DetInfo.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDigCluster.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrack.h"
#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "TotemRPValidation/RPReconstructedTracksValidation/interface/RPStationInfo.h"

#include <memory>
#include "TH1D.h"
#include "TH2D.h"
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TH2F.h>
#include <TH1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TKey.h>
#include <TLine.h>
#include "TF1.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TRint.h"
#include "TStyle.h"
#include "TMath.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "HepMC/GenEvent.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
#include "Geometry/TotemRPGeometryBuilder/interface/DetGeomDesc.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}

class PrimaryProton
{
  public: 
    HepMC::FourVector vertex;
    HepMC::FourVector momentum;
    bool found;
};


class RPReconstructedTracksValidation : public edm::EDAnalyzer
{
  public:
    explicit RPReconstructedTracksValidation(const edm::ParameterSet&);
    ~RPReconstructedTracksValidation();
    
  private:
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    
    typedef edm::DetSetVector<RPDigCluster> cluster_set;
    typedef std::map<RPId, RP2DHit> rec_tracks_collection;
    typedef std::vector<int> station_rp_ids_type;
    
    void InitBeamSmearingHistograms();
    void InitOverlapsHistograms();
    void InitStationFlux();
    void InitGeneratorDebugHistograms();
    bool FindPrimaryProtons(const edm::Handle<edm::HepMCProduct> &HepMCEvt);
    bool FindSmearedProtons(const edm::Handle<edm::HepMCProduct> &HepMCEvt);
          
    int SelectHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll, const station_rp_ids_type& st_ids);
    int Select220RightHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll);
    int Select220LeftHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll);
    int Select150RightHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll);
    int Select150LeftHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll);
    int SelectArmHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll_l, rec_tracks_collection & coll_r);
    
    bool IsRightRPId(RPId id);
    bool IsLeftRPId(RPId id);
    
    void FillReferenceHistograms(RPId id, const PrimaryProton &proton);
    void FillReferenceStationHistograms(station_rp_ids_type stations, const PrimaryProton &prim_prot_);
    void FillResidualHistograms(const rec_tracks_collection &tracks, 
        const PrimaryProton &prim_proton, 
        const PrimaryProton &smeared_proton, bool proton_at_150_station);
    
    void FillDetectorHistograms(const RPFittedTrack &track, const edm::Handle< cluster_set > &clusters);
    void FillDetectorHistogramsForReconstructedTracks( 
          const edm::Handle<RPFittedTrackCollection> &tracks, const edm::Handle<cluster_set> &clusters, const edm::Handle< edm::DetSetVector<RPDetTrigger> > &det_triggers);
    void FillEventSmearingInformation();
    void FillTrackPositions(edm::Handle< RPFittedTrackCollection > &input, edm::Handle<std::vector<RPCCBits> > &cc_chip_bits);
    void FillOverlappingAreas(const RPFittedTrackCollection &tracks);
    void FillStationFlux(const RPFittedTrackCollection &tracks);
    void FillStationProtonTracks(rec_tracks_collection & coll);
    void FillStationAcceptance(RPStationId id, double t, double xi);
    void FillStationHistograms(const rec_tracks_collection &tracks, const PrimaryProton &prim_proton, 
        const PrimaryProton &smeared_proton, bool proton_at_150_station);
    void FillGeneratorDebugHistogram_InserProtPair(const PrimaryProton &prim_prot_left_, const PrimaryProton &prim_prot_right_);

    void DrawRPContour(RPId rpid);
    void WriteHistogramAsCanvas(TH1 *h, std::vector<RPId>);
    void WriteHistograms(const std::string &root_file_name);
    void WriteBeamSmearingHistograms(TFile *f);
    void WriteOverlapsHistograms(TFile *f);
    void WriteStationFlux(TFile *f);
    void WriteGeneratorDebugHistograms(TFile *f);
    
    TotemRPGeometry Totem_RP_geometry_;
    BaseCollectionManager<RPTrackInfo, RPId, edm::ParameterSet> rp_histograms_;
    BaseCollectionManager<DetInfo, RPDetId, edm::ParameterSet> det_histograms_;
    BaseCollectionManager<RPStationInfo, RPStationId, edm::ParameterSet> station_histograms_;
    
    bool org_HepMCEvt_found_, smeared_HepMCEvt_found_;
    std::auto_ptr<ParamMADRefTransport> parameterized_madx_transport_;
    
    PrimaryProton right_prim_prot_;
    PrimaryProton left_prim_prot_;
    PrimaryProton right_smeared_prot_;
    PrimaryProton left_smeared_prot_;
    
    std::string hist_file_name_;
    std::string OriginalHepMCProductLabel_;
    std::string OriginalHepMCModuleName_;
    std::string SmearedHepMCProductLabel_;
    std::string SmearedHepMCModuleName_;
    bool SDValidation_;
    
    int verbosity_;
    
    station_rp_ids_type right_rp_ids_; 
    station_rp_ids_type left_rp_ids_;
    
    station_rp_ids_type station_150_right_ids_; 
    station_rp_ids_type station_150_left_ids_; 
    station_rp_ids_type station_220_right_ids_; 
    station_rp_ids_type station_220_left_ids_;
    station_rp_ids_type stations_left_;  //station ids in left arm
    station_rp_ids_type stations_right_;  //station ids in right arm
    
    RPFitResolution resol_degrad_service_;
    BeamOpticsParams BOPar_;
    edm::ParameterSet conf_;

    std::string modulLabelSimu_;
    std::string productLabelSimu_;
    
    std::auto_ptr<TH1D> vertex_x_smearing_dist_;
    std::auto_ptr<TH1D> vertex_y_smearing_dist_;
    std::auto_ptr<TH1D> vertex_z_smearing_dist_;
    
    std::auto_ptr<TH1D> px_right_smearing_dist_;
    std::auto_ptr<TH1D> py_right_smearing_dist_;
    std::auto_ptr<TH1D> pz_right_smearing_dist_;
    
    std::auto_ptr<TH1D> px_left_smearing_dist_;
    std::auto_ptr<TH1D> py_left_smearing_dist_;
    std::auto_ptr<TH1D> pz_left_smearing_dist_;
    
    std::auto_ptr<TH1D> thx_right_smearing_dist_;
    std::auto_ptr<TH1D> thy_right_smearing_dist_;
    
    std::auto_ptr<TH1D> thx_left_smearing_dist_;
    std::auto_ptr<TH1D> thy_left_smearing_dist_;
    
    std::auto_ptr<TH1D> xi_right_smearing_dist_;
    std::auto_ptr<TH1D> xi_left_smearing_dist_;
    
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
    std::auto_ptr<TH2F> station_flux_RP_120_121_122_;
    std::auto_ptr<TH2F> station_flux_RP_100_101_102_;
    std::auto_ptr<TH2F> station_flux_RP_020_021_022_;
    std::auto_ptr<TH2F> station_flux_RP_000_001_002_;
    
    //Generator debuging
    std::auto_ptr<TH2F> fractional_momentum_left_vs_right_;

    edm::InputTag rpFittedTrackCollectionLabel;
    edm::InputTag rpDigClusterSetLabel;
    edm::InputTag rpDetTriggerSetLabel;


};


#endif
