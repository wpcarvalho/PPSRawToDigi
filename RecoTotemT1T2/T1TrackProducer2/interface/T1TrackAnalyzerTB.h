/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
  Modified by Marcin Borratynski
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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"

#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"

//
// class decleration
//
class T1TrackAnalyzerTB : public edm::EDAnalyzer 
{


public:
	explicit T1TrackAnalyzerTB(const edm::ParameterSet&);
	~T1TrackAnalyzerTB();
	float Eta(float x,float y,float z);
	float Phi(float, float);

private:
	virtual void beginJob(const edm::EventSetup&) ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;

	double _SeeTracks;
	double _SeeHits;
	double _ChiOverNCut;
	double _ZRange;
	double _Zmin;
	double _Zmax;
	int _Verbosity;
	int eeevvv;
	double _realBeamPosX;
	double _realBeamPosY;
	double _realBeamAngle;

	TFile* theFile;


	TH1D * hRecTracksZRange;
	TH1D * hRecTracksChiCut;
	TH1D * hRecTracks;

	TH1D * hGoodTracks;
	TH1D * hGoodTracksOverAllReco;
	TH1D * hNotRecoTracks;

	TH1D * hNumOfTracksInLostEvents;
	TH1D * hNumOfTracksInEventsWnoGoodTracks;
	TH1D * hEtaOfTracksInLostEvents;

	TH1D * hAllEtaRec;
	TH1D * hAllEtaRecZRange;
	TH1D * hAllEtaRecChiCut;

	TH2D * hDEvsCHIrid;
	TH1D * hChiSquaredOverN;

	TH2D *hXYatZ5000;
	TH2D *hXYatZ0;
	TH2D *hXYatZ5000cut;
	TH2D *hXYatZm500cut;
	TH2D *hXYatZ0cut;
	TH2D *hXYatZ5700cut;
	TH2D *hXYatZ6750cut;
	TH2D *hRZ;


	T1T2TrackCollection *allTracks;
	void TwoTracksMinD(T1T2Track, T1T2Track, pair<GlobalPoint, float> & );
};
