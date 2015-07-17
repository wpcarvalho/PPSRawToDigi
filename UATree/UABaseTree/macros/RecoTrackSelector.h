#ifndef RecoSelectors_RecoTrackSelector_h
#define RecoSelectors_RecoTrackSelector_h
/* \class RecoTrackSelector
 *
 * \author Giuseppe Cerati, INFN
 *
 *  $Date: 2009/03/04 13:11:28 $
 *  $Revision: 1.1 $
 *
 */
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <iostream>
using namespace std;
class RecoTrackSelector {
 public:
  typedef reco::TrackCollection collection;
  typedef std::vector<const reco::Track *> container;
  typedef container::const_iterator const_iterator;

  /// Constructors
  RecoTrackSelector() {}
  RecoTrackSelector ( const edm::ParameterSet & cfg ) :
    ptMin_(cfg.getParameter<double>("ptMin")),
    minRapidity_(cfg.getParameter<double>("minRapidity")),
    maxRapidity_(cfg.getParameter<double>("maxRapidity")),
    tip_(cfg.getParameter<double>("tip")),
    lip_(cfg.getParameter<double>("lip")),
    minHit_(cfg.getParameter<int>("minHit")),
    min3DHit_(cfg.getParameter<int>("min3DHit")),
    maxChi2_(cfg.getParameter<double>("maxChi2")),
    maxVtxZ_(cfg.getParameter<double>("maxVtxZ")),
    ptErr_pt_(cfg.getParameter<double>("ptErr_pt")),
    ndof_(cfg.getParameter<double>("ndof")),
    bsSrc_(cfg.getParameter<edm::InputTag>("beamSpot")), 
    vertexCollection_(cfg.getParameter<edm::InputTag>("vertexCollection"))
    {
      std::vector<std::string> quality = cfg.getParameter<std::vector<std::string> >("quality");
      for (unsigned int j=0;j<quality.size();j++) quality_.push_back(reco::TrackBase::qualityByName(quality[j]));
      std::vector<std::string> algorithm = cfg.getParameter<std::vector<std::string> >("algorithm");
      for (unsigned int j=0;j<algorithm.size();j++) algorithm_.push_back(reco::TrackBase::algoByName(algorithm[j]));
    }
  
    RecoTrackSelector ( double ptMin, double minRapidity, double maxRapidity,
		      double tip, double lip, int minHit, int min3DHit, double maxChi2,
		      double maxVtxZ , double ptErr_pt, double ndof,
		      std::vector<std::string> quality , std::vector<std::string> algorithm) :
       ptMin_( ptMin ),
       minRapidity_( minRapidity ),
       maxRapidity_( maxRapidity ),
       tip_( tip ),
       lip_( lip ),
       minHit_( minHit ),
       min3DHit_( min3DHit),
       maxChi2_( maxChi2 ),
       maxVtxZ_( maxVtxZ ),
       ptErr_pt_( ptErr_pt ),
       ndof_( ndof )
    { 
      for (unsigned int j=0;j<quality.size();j++) quality_.push_back(reco::TrackBase::qualityByName(quality[j]));
      for (unsigned int j=0;j<algorithm.size();j++) algorithm_.push_back(reco::TrackBase::algoByName(algorithm[j]));
    }

  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  
  void select( const edm::Handle<collection>& c, const edm::Event & event, const edm::EventSetup&) {
    selected_.clear();
       
    //std::cout << "Hello ! I'm going to select your best tracks ... Please stay tuned ! " << endl;
    
    //cout << "Konichiwa ! Boku wa RecoTrackSelector::select() desu. Yoroshiku !" << endl;
    edm::Handle<reco::BeamSpot> beamSpot;
    event.getByLabel(bsSrc_,beamSpot); 
    
    edm::Handle< reco::VertexCollection > vertexColl;
    event.getByLabel(vertexCollection_, vertexColl);
    reco::VertexCollection theVertices = *(vertexColl.product());
    
    for( reco::TrackCollection::const_iterator trk = c->begin(); trk != c->end(); ++trk ){
      if ( operator()(*trk , beamSpot.product() , theVertices) ) {
	selected_.push_back( & * trk );
      }
    }
    
    
  }
  
  
  /// Operator() performs the selection: e.g. if (recoTrackSelector(track)) {...}
  inline  bool operator()( const reco::Track & t, const reco::BeamSpot* bs , reco::VertexCollection& vtxColl) {
  
    if(vtxColl.size() == 0) return 0;
    reco::Vertex& vtx = vtxColl.at(0);

    double dzMin = 9999.;
    for (reco::VertexCollection::const_iterator itVtx = vtxColl.begin(); itVtx != vtxColl.end(); ++itVtx) {
      if (itVtx->isFake() || itVtx->ndof() < ndof_ || fabs(itVtx->z()) > maxVtxZ_) continue;
      
      // and check this track is closer to any other vertex (if more than 1 vertex considered)
      if ( t.dz(itVtx->position()) < dzMin ) {
        vtx = *itVtx;
	dzMin = t.dz(itVtx->position());
      }
    }
    if(dzMin > 9998.) return 0;
  
  
    //------- Minimum number of hits on the track
    bool b_layerMinCut = false;
    if ( t.hitPattern().trackerLayersWithMeasurement() >= minHit_ ) b_layerMinCut = true;


    //------- pT_error / pT sel
    bool b_ptError = false;
    if (t.ptError()/t.pt() < ptErr_pt_ ) b_ptError = true;

    
    //------- Quality of the track
    bool quality_ok = true;
    if (quality_.size()!=0) {
      quality_ok = false;
      for (unsigned int i = 0; i<quality_.size();++i) {
	if (t.quality(quality_[i])){
	  quality_ok = true;
	  break;	  
	}
      }
    }


    //------- Algo Type : NOT USED
    bool algo_ok = true;
    if (algorithm_.size()!=0) {
      if (std::find(algorithm_.begin(),algorithm_.end(),t.algo())==algorithm_.end()) algo_ok = false;
    }
    
    bool keep = false;
    
    keep = ( b_layerMinCut &&
      t.hitPattern().pixelLayersWithMeasurement() + t.hitPattern().numberOfValidStripLayersWithMonoAndStereo() >= min3DHit_ &&
      b_ptError &&
      fabs(t.pt()) >= ptMin_ &&
      t.eta() >= minRapidity_ && t.eta() <= maxRapidity_ &&
      fabs(t.dxy(vtx.position())/sqrt(pow(t.dxyError(),2)+pow(vtx.xError(),2)+pow(vtx.yError(),2))) < tip_ &&
      fabs(t.dz(vtx.position())/sqrt(pow(t.dzError(),2)+pow(vtx.zError(),2))) < lip_  &&
      t.normalizedChi2()<=maxChi2_ &&
      quality_ok &&
      algo_ok);
    
    /*cout << bool( b_layerMinCut) << "  " 
         << bool( t.hitPattern().pixelLayersWithMeasurement() + t.hitPattern().numberOfValidStripLayersWithMonoAndStereo() >= min3DHit_) << "  " 
         << bool( b_ptError) << "  " 
         << bool( fabs(t.pt()) >= ptMin_) << "  " 
         << bool( t.eta() >= minRapidity_ && t.eta() <= maxRapidity_) << "  " 
         <<  fabs(t.dxy(vtx.position())/sqrt(pow(t.dxyError(),2)+pow(vtx.xError(),2)+pow(vtx.yError(),2))) << "  " 
         <<  fabs(t.dz(vtx.position())/sqrt(pow(t.dzError(),2)+pow(vtx.zError(),2))) << "  " 
         << bool( t.normalizedChi2()<=maxChi2_) << "  " 
         << bool( quality_ok) << "  " 
         << bool( algo_ok) << "  " 
	 << endl;
    
    cout << keep << "  " << t.pt() << endl;*/
    
    return keep;
  }
	
	

  size_t size() const { return selected_.size(); }
  
 protected:
  double ptMin_;
  double minRapidity_;
  double maxRapidity_;
  double tip_;
  double lip_;
  int    minHit_;
  int    min3DHit_;
  double maxChi2_;
  double maxVtxZ_;
  double ptErr_pt_;
  double ndof_;
  std::vector<reco::TrackBase::TrackQuality> quality_;
  std::vector<reco::TrackBase::TrackAlgorithm> algorithm_;
  edm::InputTag bsSrc_;
  edm::InputTag vertexCollection_;
  container selected_;
};

#endif
