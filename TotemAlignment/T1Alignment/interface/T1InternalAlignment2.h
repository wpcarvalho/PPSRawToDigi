#include <memory>
/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
*/

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DataFormats/T1Road/interface/T1Road.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Provenance/interface/BranchDescription.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"

#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "TotemAlignment/T1Alignment/interface/func2.h"
//
// class decleration
//


#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"


class T1InternalAlignment2 : public edm::EDAnalyzer
{


 public:
  explicit T1InternalAlignment2(const edm::ParameterSet&);
  ~T1InternalAlignment2();
  float Eta(float ,float ,float );
   float Phi(float,float);
   int sextant(float ,float ,float );
   int plane(float);
 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  TFile* theFile;
  int _Verbosity;
/*
  double _SeeTracks;
  double _SeeHits;
  double _ChiOverNCut;
  double _ZRange;
  double _Zmin;
  double _Zmax;
  double _realBeamPosY;
  double _realBeamPosX;
  double _realBeamAngle;




  TH1D * hRecTracksZRange;
  TH1D * hRecTracksChiCut;
  TH1D * hRecTracks;

  TH1D * hZ_at_Rmin;
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
*/
  TH1D * hDeltaX_CSC0_Pl0 ;
  TH1D *  hDeltaY_CSC0_Pl0 ;

  TH1D * hDeltaX_CSC0_Pl1 ;
  TH1D * hDeltaY_CSC0_Pl1 ;

  TH1D * hDeltaX_CSC0_Pl2 ;
  TH1D *  hDeltaY_CSC0_Pl2 ;

  TH1D *  hDeltaX_CSC0_Pl3 ;
  TH1D *  hDeltaY_CSC0_Pl3 ;

  TH1D *  hDeltaX_CSC0_Pl4 ;
  TH1D *  hDeltaY_CSC0_Pl4 ;



  TH1D * hDeltaX_CSC1_Pl0 ;
  TH1D *  hDeltaY_CSC1_Pl0 ;

  TH1D * hDeltaX_CSC1_Pl1 ;
  TH1D * hDeltaY_CSC1_Pl1 ;

  TH1D * hDeltaX_CSC1_Pl2 ;
  TH1D *  hDeltaY_CSC1_Pl2 ;

  TH1D *  hDeltaX_CSC1_Pl3 ;
  TH1D *  hDeltaY_CSC1_Pl3 ;

  TH1D *  hDeltaX_CSC1_Pl4 ;
  TH1D *  hDeltaY_CSC1_Pl4 ;



  TH1D * hDeltaX_CSC2_Pl0 ;
  TH1D *  hDeltaY_CSC2_Pl0 ;

  TH1D * hDeltaX_CSC2_Pl1 ;
  TH1D * hDeltaY_CSC2_Pl1 ;

  TH1D * hDeltaX_CSC2_Pl2 ;
  TH1D *  hDeltaY_CSC2_Pl2 ;

  TH1D *  hDeltaX_CSC2_Pl3 ;
  TH1D *  hDeltaY_CSC2_Pl3 ;

  TH1D *  hDeltaX_CSC2_Pl4 ;
  TH1D *  hDeltaY_CSC2_Pl4 ;



  TH1D * hDeltaX_CSC3_Pl0 ;
  TH1D *  hDeltaY_CSC3_Pl0 ;

  TH1D * hDeltaX_CSC3_Pl1 ;
  TH1D * hDeltaY_CSC3_Pl1 ;

  TH1D * hDeltaX_CSC3_Pl2 ;
  TH1D *  hDeltaY_CSC3_Pl2 ;

  TH1D *  hDeltaX_CSC3_Pl3 ;
  TH1D *  hDeltaY_CSC3_Pl3 ;

  TH1D *  hDeltaX_CSC3_Pl4 ;
  TH1D *  hDeltaY_CSC3_Pl4 ;



  TH1D * hDeltaX_CSC4_Pl0 ;
  TH1D *  hDeltaY_CSC4_Pl0 ;

  TH1D * hDeltaX_CSC4_Pl1 ;
  TH1D * hDeltaY_CSC4_Pl1 ;

  TH1D * hDeltaX_CSC4_Pl2 ;
  TH1D *  hDeltaY_CSC4_Pl2 ;

  TH1D *  hDeltaX_CSC4_Pl3 ;
  TH1D *  hDeltaY_CSC4_Pl3 ;

  TH1D *  hDeltaX_CSC4_Pl4 ;
  TH1D *  hDeltaY_CSC4_Pl4 ;



  TH1D * hDeltaX_CSC5_Pl0 ;
  TH1D *  hDeltaY_CSC5_Pl0 ;

  TH1D * hDeltaX_CSC5_Pl1 ;
  TH1D * hDeltaY_CSC5_Pl1 ;

  TH1D * hDeltaX_CSC5_Pl2 ;
  TH1D *  hDeltaY_CSC5_Pl2 ;

  TH1D *  hDeltaX_CSC5_Pl3 ;
  TH1D *  hDeltaY_CSC5_Pl3 ;

  TH1D *  hDeltaX_CSC5_Pl4 ;
  TH1D *  hDeltaY_CSC5_Pl4 ;






  TH1D * hDeltaMENOX_CSC0_Pl0 ;
  TH1D *  hDeltaMENOY_CSC0_Pl0 ;

  TH1D * hDeltaMENOX_CSC0_Pl1 ;
  TH1D * hDeltaMENOY_CSC0_Pl1 ;

  TH1D * hDeltaMENOX_CSC0_Pl2 ;
  TH1D *  hDeltaMENOY_CSC0_Pl2 ;

  TH1D *  hDeltaMENOX_CSC0_Pl3 ;
  TH1D *  hDeltaMENOY_CSC0_Pl3 ;

  TH1D *  hDeltaMENOX_CSC0_Pl4 ;
  TH1D *  hDeltaMENOY_CSC0_Pl4 ;



  TH1D * hDeltaMENOX_CSC1_Pl0 ;
  TH1D *  hDeltaMENOY_CSC1_Pl0 ;

  TH1D * hDeltaMENOX_CSC1_Pl1 ;
  TH1D * hDeltaMENOY_CSC1_Pl1 ;

  TH1D * hDeltaMENOX_CSC1_Pl2 ;
  TH1D *  hDeltaMENOY_CSC1_Pl2 ;

  TH1D *  hDeltaMENOX_CSC1_Pl3 ;
  TH1D *  hDeltaMENOY_CSC1_Pl3 ;

  TH1D *  hDeltaMENOX_CSC1_Pl4 ;
  TH1D *  hDeltaMENOY_CSC1_Pl4 ;



  TH1D * hDeltaMENOX_CSC2_Pl0 ;
  TH1D *  hDeltaMENOY_CSC2_Pl0 ;

  TH1D * hDeltaMENOX_CSC2_Pl1 ;
  TH1D * hDeltaMENOY_CSC2_Pl1 ;

  TH1D * hDeltaMENOX_CSC2_Pl2 ;
  TH1D *  hDeltaMENOY_CSC2_Pl2 ;

  TH1D *  hDeltaMENOX_CSC2_Pl3 ;
  TH1D *  hDeltaMENOY_CSC2_Pl3 ;

  TH1D *  hDeltaMENOX_CSC2_Pl4 ;
  TH1D *  hDeltaMENOY_CSC2_Pl4 ;



  TH1D * hDeltaMENOX_CSC3_Pl0 ;
  TH1D *  hDeltaMENOY_CSC3_Pl0 ;

  TH1D * hDeltaMENOX_CSC3_Pl1 ;
  TH1D * hDeltaMENOY_CSC3_Pl1 ;

  TH1D * hDeltaMENOX_CSC3_Pl2 ;
  TH1D *  hDeltaMENOY_CSC3_Pl2 ;

  TH1D *  hDeltaMENOX_CSC3_Pl3 ;
  TH1D *  hDeltaMENOY_CSC3_Pl3 ;

  TH1D *  hDeltaMENOX_CSC3_Pl4 ;
  TH1D *  hDeltaMENOY_CSC3_Pl4 ;



  TH1D * hDeltaMENOX_CSC4_Pl0 ;
  TH1D *  hDeltaMENOY_CSC4_Pl0 ;

  TH1D * hDeltaMENOX_CSC4_Pl1 ;
  TH1D * hDeltaMENOY_CSC4_Pl1 ;

  TH1D * hDeltaMENOX_CSC4_Pl2 ;
  TH1D *  hDeltaMENOY_CSC4_Pl2 ;

  TH1D *  hDeltaMENOX_CSC4_Pl3 ;
  TH1D *  hDeltaMENOY_CSC4_Pl3 ;

  TH1D *  hDeltaMENOX_CSC4_Pl4 ;
  TH1D *  hDeltaMENOY_CSC4_Pl4 ;



  TH1D * hDeltaMENOX_CSC5_Pl0 ;
  TH1D *  hDeltaMENOY_CSC5_Pl0 ;

  TH1D * hDeltaMENOX_CSC5_Pl1 ;
  TH1D * hDeltaMENOY_CSC5_Pl1 ;

  TH1D * hDeltaMENOX_CSC5_Pl2 ;
  TH1D *  hDeltaMENOY_CSC5_Pl2 ;

  TH1D *  hDeltaMENOX_CSC5_Pl3 ;
  TH1D *  hDeltaMENOY_CSC5_Pl3 ;

  TH1D *  hDeltaMENOX_CSC5_Pl4 ;
  TH1D *  hDeltaMENOY_CSC5_Pl4 ;


  T1T2TrackCollection *allTracks05;
  T1T2TrackCollection *allTracks01;
  T1T2TrackCollection *allTracks02;
  T1T2TrackCollection *allTracks03;
  T1T2TrackCollection *allTracks04;
  T1T2TrackCollection *allTracks00;

  T1T2TrackCollection *allTracks15;
  T1T2TrackCollection *allTracks11;
  T1T2TrackCollection *allTracks12;
  T1T2TrackCollection *allTracks13;
  T1T2TrackCollection *allTracks14;
  T1T2TrackCollection *allTracks10;


  int fitParam(int, int, T1T2TrackCollection *, double[2][5][6], double[2][5][6], double[2][5][6]);

};
