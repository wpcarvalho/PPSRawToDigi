/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
 */

// user include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Provenance/interface/BranchDescription.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"

#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"

//
// class decleration
//

class T1TrackAnalyzer : public edm::EDAnalyzer 
{
	enum {maxHits = 10000};

public:
	explicit T1TrackAnalyzer(const edm::ParameterSet&);
	~T1TrackAnalyzer();
	float Eta(float x,float y,float z);
	float Phi(float, float);

private:
	edm::InputTag simVertexContainerLabel;
	edm::InputTag simTrackContainerLabel;
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;


	double _SeeTracks;
	double _SeeHits;
	double _ChiOverNCut;
	double _ZRange;
	int eeevvv;
	int _CUTS;

	TFile* theFile;

	std::auto_ptr<TH1D> hDeltaEta;
	std::auto_ptr<TH1D> hDeltaPhi;
	std::auto_ptr<TH1D> hDeltaR;
	std::auto_ptr<TH1D> hDeltaEtaZRange;
	std::auto_ptr<TH1D> hDeltaPhiZRange;
	std::auto_ptr<TH1D> hDeltaRZRange;

	std::auto_ptr<TH1D> hDeltaEtaChiCut;
	std::auto_ptr<TH1D> hDeltaPhiChiCut;
	std::auto_ptr<TH1D> hDeltaRChiCut;

	std::auto_ptr<TH1D> hEventType;
	std::auto_ptr<TH1D> hEventTypeNoT1;
	std::auto_ptr<TH1D> hEventTypeHalfT1;
	std::auto_ptr<TH1D> hTriggerType;
	std::auto_ptr<TH1D> hTriggerTypeNoT1;
	std::auto_ptr<TH1D> hTriggerTypeHalfT1;
	std::auto_ptr<TH1D> hNTracks;
	std::auto_ptr<TH1D> hNTracksNoT1;
	std::auto_ptr<TH1D> hNTracksHalfT1;
	std::auto_ptr<TH1D> hNTracksLeft;
	std::auto_ptr<TH1D> hNTracksLeftNoT1;
	std::auto_ptr<TH1D> hNTracksLeftHalfT1;
	std::auto_ptr<TH1D> hNTracksRight;
	std::auto_ptr<TH1D> hNTracksRightNoT1;
	std::auto_ptr<TH1D> hNTracksRightHalfT1;

	std::auto_ptr<TH1D> hSimTracks;
	std::auto_ptr<TH1D> hRecTracksZRange;
	std::auto_ptr<TH1D> hRecTracksChiCut;
	std::auto_ptr<TH1D> hRecTracks;
	std::auto_ptr<TH1D> hRecSimTrackRatio;
	std::auto_ptr<TH1D> hRecZRangeSimTrackRatio;
	std::auto_ptr<TH1D> hGoodTracks;
	std::auto_ptr<TH1D> hGoodTracksOverAllReco;
	std::auto_ptr<TH1D> hNotRecoTracks;
	std::auto_ptr<TH1D> hNotRecoTracksRatio;
	std::auto_ptr<TH1D> hNumOfTracksInLostEvents;
	std::auto_ptr<TH1D> hNumOfTracksInEventsWnoGoodTracks;
	std::auto_ptr<TH1D> hEtaOfTracksInLostEvents;
	std::auto_ptr<TH1D> hAllEtaSim;
	std::auto_ptr<TH1D> hAllEtaSimT1Range;
	std::auto_ptr<TH1D> hAllEtaRec;
	std::auto_ptr<TH1D> hAllEtaRecZRange;
	std::auto_ptr<TH1D> hAllEtaRecChiCut;

	std::auto_ptr<TH2D> hDEvsCHIrid;
	std::auto_ptr<TH1D> hChiSquaredOverN;

	std::auto_ptr<TH1D> hTrackNumberSingleArmCondition2;
	std::auto_ptr<TH1D> hTrackNumberDoubleArmCondition2;
	std::auto_ptr<TH1D> hTrackNumberSingleArmCondition12;
	std::auto_ptr<TH1D> hTrackNumberDoubleArmCondition12;
	std::auto_ptr<TH1D> hTrackNumberSingleArmConditionHalf12;
	std::auto_ptr<TH1D> hTrackNumberDoubleArmConditionHalf12;

	TTree* tree;

	unsigned int nTr;
	double eeTr[maxHits],  ffTr[maxHits],  chTr[maxHits], pTr[maxHits];
	int evto[maxHits], tyTr[maxHits];
};
