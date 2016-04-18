
#include "RecoTotemRP/RPRecoDataFormats/interface/RPMulTrackCandidateCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateDistinctCollectionsSet.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPMulFittedTrackCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
#include <map>
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProton.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RP2DHit.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RP2DHitDebug.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonPair.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonPairCollection.h"

#include <vector>
#include "RecoTotemRP/RPRecoDataFormats/interface/CentralMassInfo.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPStationTrackFit.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPStationTrackFitCollection.h"



namespace {
  namespace {
    
    CentralMassInfo ctrinfo;
    edm::Wrapper<CentralMassInfo> wctrinfo;
    
    RPMulTrackCandidateCollection mcoll;
    RPTrackCandidateDistinctCollectionsSet coll_set;
    edm::Wrapper<RPMulTrackCandidateCollection> RPMulTrackCandidateCollectionWrapper;
    edm::Wrapper<RPTrackCandidateDistinctCollectionsSet> RPTrackCandidateDistinctCollectionsSetWrapper;
    
    edm::Wrapper<RPMulFittedTrackSetsCollection> wraprpmulttracksetcol; 
    
    RPMulTrackCandidateSetsCollection rpmultracsetcandcol;
    edm::Wrapper<RPMulTrackCandidateSetsCollection> wraprpmultracsetcandcol;
    
    std::less<unsigned int> lui;
    std::binary_function<unsigned int,unsigned int,bool> bf;

    RPMulFittedTrackCollection the_mtrack_cand_col;
    edm::Wrapper<RPMulFittedTrackCollection> the_w_mrpftc;
    
    RPReconstructedProton rprecprot;
    edm::Wrapper<RPReconstructedProton> wraprprecprot;
    std::vector<RPReconstructedProton> prprecprotvec;
    
    RPReconstructedProtonCollection rprpcol;
    edm::Wrapper<RPReconstructedProtonCollection> wrrprpcol;
    
    RP2DHit rp2dhit;
    edm::Wrapper<RP2DHit> wrp2dhit;
    RP2DHitDebug rp2debugdhit;
    edm::Wrapper<RP2DHitDebug> rwp2debugdhit;
    
    std::map<unsigned int, RP2DHitDebug> vrp2debugdhit;
    edm::Wrapper<std::map<unsigned int, RP2DHitDebug> > wvrp2debugdhit;
    std::pair<unsigned int, RP2DHitDebug> vrp2debughitpair;
    
    RPReconstructedProtonPair p;
    edm::Wrapper<RPReconstructedProtonPair> wp;
    std::vector<RPReconstructedProtonPair> vp;
    
    RPReconstructedProtonPairCollection pc;
    edm::Wrapper<RPReconstructedProtonPairCollection> wpc;

	RPStationTrackFit rpstf;
	std::vector<RPStationTrackFit> vrpstf;
	std::map<unsigned int, std::vector<RPStationTrackFit> > m_ui_vrpstf;
	RPStationTrackFitCollection rpstfc;
	edm::Wrapper<RPStationTrackFitCollection> w_rpstfc;

  }
}
