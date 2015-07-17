// -*- C++ -*-
//
// Package:    TotemRecotrackClassifier
// Class:      TotemRecotrackClassifier
// 
/**\class TotemRecotrackClassifier TotemRecotrackClassifier.h TotemT1T2Validation/TotemRecotrackClassifier/interface/TotemRecotrackClassifier.h

 Description: Classifies reconstructed tracks into primary and secondary tracks based on simulated information.

 Implementation:

*/
//
// Original Author:  Matti Leinonen T2BackgroundAnalysis.h
//         Created:  Thu Aug 19 13:39:39 CEST 2010
// $Id$
//
//

#ifndef _T2BackgroundAnalysis_h
#define _T2BackgroundAnalysis_h
// system include files
#include <memory>
#include <string>

// user include files


#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "FWCore/Utilities/interface/InputTag.h"



//#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Framework/interface/EventSetup.h"



//#include "FWCore/Framework/interface/ESHandle.h"
#include "TotemAnalysis/T2Cuts/interface/T2SelectionCutUtils.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/T1T2Track/interface/T1T2Track.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"

#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2StripClusterCollection.h"


//#include "Geometry/TotemGeometry/interface/T1Geometry.h"

//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
//Root classes
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"




//#include "HepPDT/defs.h"

//#include "HepPDT/TableBuilder.hh"

//#include "HepPDT/ParticleDataTable.hh"

//#include "HepMC/GenEvent.h"



//
// class decleration
//

class T2BackgroundAn : public edm::EDAnalyzer {
public:
  explicit T2BackgroundAn(const edm::ParameterSet&);
  ~T2BackgroundAn();
  virtual void beginJob();

  virtual void analyze(const edm::Event&,
		       const edm::EventSetup&
		       );

  virtual void endJob(); 

private:
 


  void Trigger(const T2PadDigiCollection* PadDigiptr, bool & H0trigger, bool & H1trigger, bool & H2trigger, bool & H3trigger);

  unsigned int  etaToPos(double centre);
  void ZImpactRangeFromFit(int quarter, double eta2, double & ZImpLeft, double & ZImpRight);

  double PartCharge(int pdgcode);
  bool VertexFromIonPump(double Zvtxposition);

  bool IsAcceptedForDNDeta(T1T2Track Trk);

  double GetEta2FromGeantHits(std::vector<std::pair<double,double> > SimHitRZ_ofGeantPrimaryTrk);

  bool CUTfordNdeta(T1T2Track Trk);
 
  bool CUTfordNdetaVtx(T1T2Track trk, bool fitxy);

  void PrintParents(HepMC::GenEvent::particle_const_iterator p);

  bool ThisParticleHasStableDaughters(HepMC::GenEvent::particle_const_iterator p);

  bool AreTrkFromIP_PyChParticle(int Hit_TrkId,const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, const HepMC::GenEvent* evt);
  
  bool AreTrkFrom_K0(int Hit_TrkId,const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, const HepMC::GenEvent* evt);

  unsigned int RawtoSymb(uint32_t thedet);

  bool PrimaryTrackCondition(SimTrack atrk,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
			     const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,int &barcodeMother);

  
  int GetTrkOlderMotherPid(SimTrack atrk,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
				  const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,int &barcodeMother,int &DirectTrkPID);

  int getPythiaPID(const HepMC::GenEvent* evt,int &PyPid,int OldestmotherGeantId);

  int  GeantTrkReconstructed(const edm::Event& event,
			     const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
			     const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
			     const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
			     const std::map<unsigned int, std::list<Local3DPoint> >& trackHitList,
			     const std::auto_ptr<T1T2TrackCollection> &theT2Tracks2,
			     int GeantTrkIndex, 
			     bool printActive,
			     unsigned int &NumGeantTrk_PassingCuts,
			     unsigned int &NumGeantTrk_PassingCuts_56Division,
			     unsigned int &NumGeantTrk_PassingCutsDiscrepancyZCutCond,
			     unsigned int &NumGeantTrk_PassingCutsDiscrepancy,
			     unsigned int &NumGeantTrk_PassingCuts_2, 
			     /*std::vector<Local3DPoint> &PrimaryGeantHitPos,*/
			     T1T2Track &TheGeantTrkReconstr,bool &GeantTrkFound,
			     double &effiasRecRec,
			     int &BestTrkIdMatchedPosition,
			     double &lastpassingZmin, double &lastpassingZImpact,
			     bool discrepancyStudies=false
			     );

  
  
  
  
  // T1T2Track MyLinearfitCorr(std::vector<T2Hit> hitvec,TMatrixD &par_covariance_matrix,double &chi2_,bool Usestrip,std::vector<double> trkparam, int RoadID);
  T1T2Track MyLinearfitCorrDEV(std::vector<T2Hit> hitvec2,TMatrixD &par_covariance_matrix,double &chi2_,bool Usestrip, int RoadID);

  T1T2Track CalculateRecoTrkFromGeant(std::vector<Local3DPoint> PrimaryGeantHitPos, std::vector<double> &theZ,  std::vector<unsigned int> &thecmsswId);

  bool NeutralDecay(int pid);

  void TrkDeltaTheta(T1T2Track &atrk,double &DeltaThetaX,double &DeltaThetaY);

  bool RecoTrk_UnfoldingCut(T1T2Track &atrk);

  bool RecoTrk_UnfoldingCut_56Division(T1T2Track &atrk);

  bool RecoTrk_UnfoldingCut_2(T1T2Track &atrk);
  /*
    Compares reconstructed tracks against simualted tracks and
    classifies reconstructed tracks into primary and simulated tracks.
   */
  void
  MatchRecotracksToSimulatedTracks(const edm::Event& event,
				   const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
				   const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
				   const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
				   const std::map<unsigned int, std::list<Local3DPoint> >& trackHitList,
				   const std::auto_ptr<T1T2TrackCollection> &theT2Tracks2
				   );

  std::vector<int> RecotracksInfoFromSimu(const edm::Event& event,
					  const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
					  const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
					  const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
					  const std::map<unsigned int, std::list<Local3DPoint> >& trackHitList,
					  const HepMC::GenEvent* evt,
					  T1T2Track recTrack,
					  std::vector<double> &PrimaryTrkenergy_ForTheRecoHits
					  );


  std::vector<int>  RecotracksInfoFromSimu_BIS(const edm::Event& event,
					  const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
					  const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
					  const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
					  const std::map<unsigned int, std::list<Local3DPoint> >& trackHitList,
					  const HepMC::GenEvent* evt,
					  T1T2Track recTrack,
					  std::vector<double> &PrimaryTrkenergy_ForTheRecoHits
					  );

  T1T2Track RPhiFit(std::vector<T2Hit> &hitvec2);

  double MostProbableTrkEnergy(std::vector<double> &PrimaryTrkenergy_ForTheRecoHits);

  /*
    Loads corresponding trackproducer tracks and returns number of tracks read
  */
  unsigned int loadTrackproducerRaw(const edm::Event&,
					std::auto_ptr<T1T2TrackCollection>&
					);
  
  //unsigned int loadTrackproducer1Tracks(const edm::Event&, std::auto_ptr<T1T2TrackCollection>&);
  unsigned int loadTrackproducer2Tracks(const edm::Event&,
					std::auto_ptr<T1T2TrackCollection>&
					);

  /*
    Loads simulated tracks, vertices and hits and returns number of items read
  */

  unsigned int loadSimulatedTracks(const edm::Event&,
				   std::auto_ptr<edm::SimTrackContainer>&
				   );
   
  unsigned int loadSimulatedVertices(const edm::Event&,
				     std::auto_ptr<edm::SimVertexContainer>&
				     );

  unsigned int loadSimulatedHits(const edm::Event&,
                                 const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
                                 const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
                                 std::auto_ptr<edm::PSimHitContainer>& theSimHits,
                                 std::map<unsigned int, std::list<Local3DPoint> >& trackHitList,
                                 std::map<unsigned int, int >& unknownTrackList
                                 );

  /*
    Calculates correct XYZ pos from simhit pos
  */
  Local3DPoint CalculateCorrectXYZPos(const Local3DPoint &pos, const unsigned int detUnitId);

  /*
    Returns track corresponding to trackId
  */
  bool getTrack(const unsigned int trackId,
                const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
                SimTrack &track);

  /*
    Returns vertex corresponding to vertexId
  */
  bool getVertex(const  int vertexId,  const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, SimVertex &vertex);

  /*
  bool getVertex(const unsigned int vertexId,
                 const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
                 SimVertex &thevertex);
  */
  /*
    Checks if given TRACKID corrensponds to the primary vertex.
    (Not a secondary track)
    reurns true if TRACKID corrensponds to the primary vertex otherwise false
    also returns false if no vertex found
  */

  bool isTrackPrimary(const int trackId,
                      const std::auto_ptr<edm::SimTrackContainer>&,
                      const std::auto_ptr<edm::SimVertexContainer>&);
  
  double ZimpactFromRecTrack(T1T2Track &mytrk);
  
  double  CalculateEfficiencyAsRecRec(T1T2Track &TheGeantTrkReconstr, T1T2Track &atrk);
  double GeneratorEtaFromTrkBarcode(int barcode,const HepMC::GenEvent* evt);					
  
  edm::InputTag T2PadDigiCollectionLabel;
  edm::InputTag T2StripDigiCollectionLabel;
edm::InputTag SimTrackContainerLabel;
edm::InputTag SimVertexContainerLabel;

  std::vector<T1T2Track> t2trackVectorBothSide; 
  T2SelectionCutUtils T2CutsUtil;
      
  std::vector<T1T2Track> t2trackVectorBothSideRAW;    
  // Histograms
  std::auto_ptr<TH1D> H0H1_PhiGeantTrkWhenEffiFail;
  std::auto_ptr<TH2D> H0H1_EtaPhiGeantTrkWhenEffiFail;
  std::auto_ptr<TH1D> PdgSecondaryInT2EtaAcceptance;
  std::auto_ptr<TH2D> RecoVtx_ZvsR;
  std::auto_ptr<TH2D> SecondaryVtx_ZvsR_ContributingInT2;

  std::auto_ptr<TH1D> FractionOfSecondaryTrksRecoAsPrimary;
  std::auto_ptr<TH1D> FractionOfSecondaryTrksInPrimaryVtx;
  
  std::auto_ptr<TH1D> FractionOfPrimaryTrksInSecVtx;
  std::auto_ptr<TH1D> FractionOfPrimaryTrksRecoAsSecondary;
  
  std::auto_ptr<TH1D> NumRecoPrimVtx;

  std::auto_ptr<TH2D> AllRecoVtxYZ;
  std::auto_ptr<TH1D> AllRecoVtxZ;

  std::auto_ptr<TH2F> ParticlePdgvsEinT2;
  std::auto_ptr<TH2F> ParticlePdgvsENotinT2;
  
  std::auto_ptr<TH2F> PrimaryTrksEtavsZatRmin;
  std::auto_ptr<TH1D> PrimaryTrksEta;
  std::auto_ptr<TH1D> AllTrksEta;
  std::auto_ptr<TH1D> SecondaryTrksEta;

  std::auto_ptr<TH1D>CumulativeZimpH1;  
  std::auto_ptr<TH1D>CumulativeZimpH0;
  std::auto_ptr<TH1D>CumulativeZimpH2;  
  std::auto_ptr<TH1D>CumulativeZimpH3;
  std::auto_ptr<TH1D>rechitUnitID;
  std::auto_ptr<TH1D>simhitUnitID;

  std::auto_ptr<TProfile>  DigiPadOccupancy;
  std::auto_ptr<TH1D> RecoTrksEtaALL;
  std::auto_ptr<TH1D> RecoTrksZ0ALL;
  std::auto_ptr<TH1D> RecoTrksR0ALL;

  std::auto_ptr<TH1D> RecoTrksZ0ALL_R0Cut10;
  std::auto_ptr<TH1D> RecoTrksZ0ALL_R0Cut20;
  
  std::auto_ptr<TH1D> RecoTrksZ0ALL_EtaCut4;
  std::auto_ptr<TH2D> SecondarySmallZ0_EtavsR0;

std::auto_ptr<TH1D>   eta2geant_forEffHisto1;
std::auto_ptr<TH1D>   eta2geant_forEffHisto2;
std::auto_ptr<TH1D>   eta2geant_forEffHisto3;



  std::auto_ptr<TH1D> SecondaryRecoTrksZ0_2Arms_R0Cut20;
  std::auto_ptr<TH1D> SecondaryRecoTrksZ0_2Arms_EtaCut4_R0Cut20; 
  std::auto_ptr<TH1D> RecoTrksZ0ALL_EtaCut4_R0Cut20; 
  
  std::auto_ptr<TH2D> SecondaryBigZ0_EtavsR0; 


  std::auto_ptr<TH2D> Reco_TrkZ0_vs_R0_AllH0_FromPrimary; 
  std::auto_ptr<TH2D> Reco_TrkZ0_vs_Eta_AllH0_FromPrimary;
  std::auto_ptr<TH2D> Reco_TrkZ0_vs_Rmin_AllH0_FromPrimary;
  

  std::auto_ptr<TH2D> Reco_TrkZ0_vs_R0_AllH0_FromSecondary; 
  std::auto_ptr<TH2D> Reco_TrkZ0_vs_Eta_AllH0_FromSecondary; 
  std::auto_ptr<TH2D> Reco_TrkZ0_vs_Rmin_AllH0_FromSecondary;
  
  std::auto_ptr<TH2D> Reco_Trk_DeltaThetaX_vs_DeltaThetaY_FromSecondary; 
  std::auto_ptr<TH2D> Reco_Trk_DeltaThetaX_vs_DeltaThetaY_FromPrimary; 

   
  
  std::auto_ptr<TH2D> SecondarySmallRminZVtxPosition;
  std::auto_ptr<TH1D> NumSecPerEvent;
  std::auto_ptr<TH2D> DrDphiSecondary;



  std::auto_ptr<TH2D>   SecondaryTrksZvsEta; 
  std::auto_ptr<TH2D>   SecondaryTrksZvsPId;
  
  std::auto_ptr<TH1D> SecondaryTrksEtaFromIonPump; 

  std::auto_ptr<TH1D> PID_ofSecRecoTrack;
  std::auto_ptr<TH1D> PID_ofSecRecoTrackAsaFakePrimary;
  std::auto_ptr<TH1D> PID_DircetMotherForAFakeRecoPrimaryTrk;
  std::auto_ptr<TH1D> PID_OldestMotherForAFakeRecoPrimaryTrk;
  std::auto_ptr<TH1D> PID_GeantTrk_atleast4HitInT2;
  std::auto_ptr<TH1D> PID_GeantTrk_atleast4HitInT2_StableAtIP;
  std::auto_ptr<TH1D> PID_GeantTrk_atleast4HitInT2_ConsideredAsPrimary;
  std::auto_ptr<TH1D> Geant_Prim_Or_Sec_Unfold;
  


  std::auto_ptr<TH2D> NumRecoSecondary_VsNumRecoPrimary; 
  std::auto_ptr<TH2D> NumRecoSecondary_VsNumRecoPrimary_EtaCut4_R0Cut20;
  std::auto_ptr<TH2D> SimuVtxPositionForAFakeRecoPrimaryTrk;

  std::auto_ptr<TH1D> SecondaryRecoTrksEta;
  std::auto_ptr<TH1D> PrimaryRecoTrksEta;

  std::auto_ptr<TH1D> SecondaryRecoTrksEtaPointingToVtx;
  std::auto_ptr<TH1D> SecondaryRecoTrksR0PointingToVtx;
  std::auto_ptr<TH1D> SecondaryRecoTrksZ0PointingToVtx;
  
  std::auto_ptr<TH1D> SecondaryRecoTrksR02PointingToVtx;
  std::auto_ptr<TH1D> PrimaryRecoTrksZ0PointingToVtx;
  std::auto_ptr<TH1D> PrimaryRecoTrksrR02PointingToVtx;
  std::auto_ptr<TH1D> PrimaryRecoTrksrR0PointingToVtx;
   

  std::auto_ptr<TH1D> RecoTrksEtaALL_EtaCut4_R0Cut20;  
  std::auto_ptr<TH1D> SecondaryRecoTrksEta_2Arms_EtaCut4_R0Cut20;  
  std::auto_ptr<TH1D> PrimaryRecoTrksEta_2Arms_EtaCut4_R0Cut20; 
  std::auto_ptr<TH1D> PrimaryRecoTrksZ0_2Arms_EtaCut4_R0Cut20; 
  std::auto_ptr<TH1D> RecoTrksR0ALL_EtaCut4_R0Cut20;  
  std::auto_ptr<TH1D> SecondaryRecoTrksR0_2Arms_EtaCut4_R0Cut20;  
  std::auto_ptr<TH1D> PrimaryRecoTrksR0_2Arms_EtaCut4_R0Cut20; 

  std::auto_ptr<TProfile>  TrkRecoEfficiencyProfilePlus; 
  std::auto_ptr<TH1D>  GeantTrkMultiplicityfrequencyPlus; 
  std::auto_ptr<TProfile>  TrkRecoEfficiencyProfileMinus; 
  


  std::auto_ptr<TH1D> SecondaryRecoTrksR02_2Arms;
  std::auto_ptr<TH1D> SecondaryRecoTrksR0_2Arms;
  std::auto_ptr<TH1D> SecondaryRecoTrksZ0_2Arms;
  std::auto_ptr<TH1D> PrimaryRecoTrksR02_2Arms;
  std::auto_ptr<TH1D> PrimaryRecoTrksR0_2Arms;
  std::auto_ptr<TH1D> PrimaryRecoTrksZ0_2Arms;
  std::auto_ptr<TH1D> PrimaryRecoTrksZ0_2ArmsEsmall5;

  std::auto_ptr<TH1D> RecoTrkPrimaryNumHitEnergy;
  std::auto_ptr<TH1D> RecoTrkPrimaryMostProbEnergy;

  std::auto_ptr<TH1D> H0_PrimaryRecoTrksZ0_2ArmsEbig5;
  std::auto_ptr<TH1D> H0_PrimaryRecoTrksZ0_2Arms;

  std::auto_ptr<TH1D> H0_PrimaryRecoTrksZ0_2Arms_BisMatch;
  std::auto_ptr<TH1D> H0_PrimaryZ0_OrtogImpact_BisMatch;
  std::auto_ptr<TH2D> H0_PrimaryXY_OrtogImpact_BisMatch;

  std::auto_ptr<TH2D> H0_SecondaryRecoTrksZ0_Vs_Z0OrtogImpact_BisMatch;
  std::auto_ptr<TH2D> H0_PrimaryRecoTrksZ0_Vs_Z0OrtogImpact_BisMatch;

  std::auto_ptr<TH2D> H0_SecondaryRecoTrksZ0Ortog_Vs_Y0Ortog_BisMatch;
  std::auto_ptr<TH2D> H0_PrimaryRecoTrksZ0Ortog_Vs_Y0Ortog_BisMatch;
  std::auto_ptr<TH2D> H0_SecondaryRecoTrksZ0Ortog_Vs_X0Ortog_BisMatch;
  std::auto_ptr<TH2D> H0_PrimaryRecoTrksZ0Ortog_Vs_X0Ortog_BisMatch;
  std::auto_ptr<TH2D> H0_SecondaryRecoTrksX0Ortog_Vs_Y0Ortog_BisMatch;
  std::auto_ptr<TH2D> H0_PrimaryRecoTrksX0Ortog_Vs_Y0Ortog_BisMatch;
  std::auto_ptr<TH2D> H0_SecondaryRecoTrksZ0Ortog_Vs_R0Ortog_BisMatch;
  std::auto_ptr<TH2D> H0_PrimaryRecoTrksZ0Ortog_Vs_R0Ortog_BisMatch;
  std::auto_ptr<TH2D> H0_EtaPhiGeantTrkWhenEffiFailMultLess20;
  
  std::auto_ptr<TH1D> H0_TrkPrimIDWhenEffiFail;
  std::auto_ptr<TProfile> H0_HowManySimHitWhenEffiFail;
  std::auto_ptr<TH1D> H0_EtaGeantTrkWhenEffiFailMultLess20;

  std::auto_ptr<TH1D> H0_RecoTrksZ0;
   


  std::auto_ptr<TH1D> H0_PrimaryRecoTrksBX_2Arms;
  std::auto_ptr<TH1D> H0_PrimaryRecoTrksBY_2Arms;
  std::auto_ptr<TH1D> H0_SecondaryRecoTrksBX_2Arms;
  std::auto_ptr<TH1D> H0_SecondaryRecoTrksBY_2Arms;
  std::auto_ptr<TH1D> H0_RecoTrksEta;
  std::auto_ptr<TH1D> H0_SecondaryRecoTrksEta;
  std::auto_ptr<TH1D> H0_RecoTrksEtaR0Small20;
  std::auto_ptr<TH1D> H0_SecondaryRecoTrksEtaR0Small20;


  std::auto_ptr<TH1D> H0_PrimaryRecoTrksR0_2ArmsEbig5; 
  std::auto_ptr<TH1D> H0_PrimaryRecoTrksEta_2ArmsEbig5; 
  std::auto_ptr<TH1D> H0_PrimaryRecoTrksZ0_2ArmsEbig5R0Less20;
  std::auto_ptr<TH1D> H0_PrimaryRecoTrksEta_2ArmsEbig5R0Less20;
  
  std::auto_ptr<TH2D> H0_PrimaryXY_OrtogImpactEbig5;
  std::auto_ptr<TH2D> H0_PrimaryXY_OrtogImpact; 
  std::auto_ptr<TH2D> H0_XY_OrtogImpact;

  std::auto_ptr<TH1D> H0_PrimaryZ0_OrtogImpact;
  std::auto_ptr<TH1D> H0_Z0_OrtogImpact;
  std::auto_ptr<TH1D> H0_PrimaryZ0_OrtogImpactEbig5;
  std::auto_ptr<TH1D> PrimaryRecoTrksZ0_2ArmsEbig5;
  std::auto_ptr<TH1D> PrimaryRecoTrksZ0_2ArmsMixedEnergy;
  
 std::auto_ptr<TProfile> CutZimpRightH0;
  std::auto_ptr<TProfile> CutZimpLeftH0;


  std::auto_ptr<TProfile> CutZimpRightH1;
  std::auto_ptr<TProfile> CutZimpLeftH1;

  std::auto_ptr<TProfile> CutZimpRightH2;
  std::auto_ptr<TProfile> CutZimpLeftH2;


  std::auto_ptr<TProfile> CutZimpRightH3;
  std::auto_ptr<TProfile> CutZimpLeftH3;



  std::auto_ptr<TH2D> PrimaryVertexPositionYZ;
  std::auto_ptr<TH2D> SecondaryVertexPositionYZ;
  std::auto_ptr<TH2D> ALLVertexPositionYZ;
  std::auto_ptr<TH2D> SecondaryVtxPositionContributingInT2;
  std::auto_ptr<TH2D> SecondaryVtxPositionContributingInT2MultiCount;
  std::auto_ptr<TH2D> PrimaryVtxPositionContributingInT2;


  //std::auto_ptr<TProfile> ProbabilityToRecoRealPrimary_vs_Eta; //Trk Efficiency
  //std::auto_ptr<TProfile> ProbabilityToRecoFakePrimary_vs_Eta; //BackGroung Contamination
  
  std::auto_ptr<TProfile> FakeSecondaryRecoDnDeta2FromVtx; //Trk Efficiency
  std::auto_ptr<TProfile> PrimaryRecoDnDeta2FromVtx; //Trk Efficiency
  std::auto_ptr<TProfile> PrimaryRecoDnDeta2NoCut;
  std::auto_ptr<TProfile> RecoDnDeta2FromVtx; //BackGroung Contamination
  std::auto_ptr<TProfile> PrimaryGeantDNDeta2; //BackGroung Contamination
  std::auto_ptr<TProfile> PrimaryGeantDNDeta2_IfOneTrkInT2; 
  std::auto_ptr<TProfile> PrimaryGeantDNDeta2_IfOneRecoZImpInH0;
  std::auto_ptr<TProfile> PrimaryGeantDNDeta2_IfOneTrkInT2_H1;
  std::auto_ptr<TProfile> PrimaryGeantDNDeta2_IfOneTrkInT2_H1sometrkexcluding64;
  std::auto_ptr<TProfile> PrimaryGeantDNDeta2_IfOneTrkInT2_H0;
  std::auto_ptr<TProfile>   PrimaryGeantDNDeta2H0_ZCut;
  std::auto_ptr<TProfile> DNDetaMBALLMCGenerator_GeneratorTriggered;
  std::auto_ptr<TProfile> DNDetaMBALLMCGenerator_GeneratorTriggeredAtLeastOne5364;
  
  std::auto_ptr<TProfile> RecoDnDeta2NoCut; 
  std::auto_ptr<TProfile> RecoDnDeta2UnfoldCut2;
  std::auto_ptr<TProfile> RecoDnDeta2UnfoldCut2_Primary;
  std::auto_ptr<TProfile> RecoDnDeta2UnfoldCut;
  std::auto_ptr<TProfile> RecoDnDeta2UnfoldCut_Primary;

  std::auto_ptr<TH1D> PrimaryGeanteta2;
  std::auto_ptr<TH1D>  VtxTrkAssociation;

  std::auto_ptr<TH1D> GeantThetaXofPrimary;std::auto_ptr<TH1D> GeantThetaYofPrimary;
  std::auto_ptr<TH1D>   NumHitInRecTrk;
  std::auto_ptr<TH1D> Pythia_UnstablePDG; std::auto_ptr<TH1D> Pythia_PDG_HavingStableDesc;
  std::auto_ptr<TH1D> Pythia_PDG_Status2_Theta003; std::auto_ptr<TH1D> Pythia_PDG_HavingStableDesc_Theta003;
  std::auto_ptr<TH1D> Pythia_PDG_HavingStableDesc_Theta001;
  std::auto_ptr<TH1D> PlaneClusterActiveCounter;
  std::auto_ptr<TH1D> PrimTrkTestedQuarter; 
  std::auto_ptr<TH1D> PrimTrkTestedPlaneQuarter;
  std::auto_ptr<TH1D> PrimTrkTestedQuarterUnamb;
  std::auto_ptr<TH1D> H0_sameTrk_StripNumb;
  std::auto_ptr<TH1D> H0_sameTrk_PadRow;  
  std::auto_ptr<TH1D> H0_sameTrk_PadCol;  
  std::auto_ptr<TH1D> H0_BinMultStat;   std::auto_ptr<TH1D> H1_BinMultStat;  
  std::auto_ptr<TH1D> H2_BinMultStat;  std::auto_ptr<TH1D> H3_BinMultStat;   
  std::auto_ptr<TH1D> H0_BinMultStatFail;
  std::auto_ptr<TH1D> H1_BinMultStatFail;
  std::auto_ptr<TH1D> H2_BinMultStatFail;
  std::auto_ptr<TH1D> H3_BinMultStatFail;

  std::auto_ptr<TH1D> H0_sameTrk_PadColMaxDist;
  std::auto_ptr<TH1D> H0_sameTrk_PadRowMaxDist;
  std::auto_ptr<TH1D> H0_sameTrk_StripColMaxDist; 
  std::auto_ptr<TH1D> H0_sameTrk_StripRowMaxDist;

  std::auto_ptr<TH1D> eta2geant_forEffQ0;
  std::auto_ptr<TH1D> eta2geant_forEffQ1;
  std::auto_ptr<TH1D> eta2geant_forEffQ2;
  std::auto_ptr<TH1D> eta2geant_forEffQ3;
   std::auto_ptr<TH1D> eta2geant_forEffQ0eff;
  std::auto_ptr<TH1D> eta2geant_forEffQ1eff;
  std::auto_ptr<TH1D> eta2geant_forEffQ2eff;
  std::auto_ptr<TH1D> eta2geant_forEffQ3eff;
  

  std::auto_ptr<TH2F> SelectedEventGeantPad;
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3Hit;
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3HitEntry_Weight; 
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3HitExit_Weight; 
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3Tracing;
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3TracingMultiCount;
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3TracingGEO;
  std::auto_ptr<TH2F> SelectedEventGeantEntryPad; 
  std::auto_ptr<TH2F> SelectedEventGeantExitPad;
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY;
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_Noelectron;
  
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3_PythiaChinIP;
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3_PythiaChinIPGEOXY;

  
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXY[50];
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEORPhi[50];
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOXZ[50];
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3_PythiaK0GEOYZ[50];
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsPadRow[50];
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadColvsZ[50]; 
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3_PythiaK0PadRowvsZ[50];
  

  std::auto_ptr<TH2F> SelectedEventDIGIPad_GTrkWithMoreThan3_GEOYZ[50]; 
  std::auto_ptr<TH2F> SelectedEventDIGIPad_GTrkWithMoreThan3_GEOXZ[50]; 
  std::auto_ptr<TH2F> SelectedEventDIGIPad_GTrkWithMoreThan3_GEOXY[50]; 
  std::auto_ptr<TH2F> SelectedEventDIGIPad_GTrkWithMoreThan3_PadColvsPadRow[50];
  std::auto_ptr<TH2F> SelectedEventDIGIPad_GTrkWithMoreThan3_PadColvsZ[50]; 
  std::auto_ptr<TH2F> SelectedEventDIGIPad_GTrkWithMoreThan3_PadRowvsZ[50]; 

  std::auto_ptr<TH2F> SelectedEventRoadPadFinderClu_PadRowvsZ[50]; 
  std::auto_ptr<TH2F> SelectedEventRoadPadFinderClu_PadColvsZ[50]; 
  std::auto_ptr<TH2F> SelectedEventRoadPadFinderClu_GEOXZ[50]; 
  std::auto_ptr<TH2F> SelectedEventRoadPadFinderClu_GEOYZ[50]; 
  

 

  std::auto_ptr<TH2F> SelectedEventRoadPadFinder_FullR_Clu_PadRowvsZ[50]; 
  std::auto_ptr<TH2F> SelectedEventRoadPadFinder_FullR_Clu_PadColvsZ[50]; 
  std::auto_ptr<TH2F> SelectedEventRoadPadFinder_FullR_Clu_GEOXZ[50]; 
  std::auto_ptr<TH2F> SelectedEventRoadPadFinder_FullR_Clu_GEOYZ[50];
  

  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3TracingGEOXZ_NoelectronPrimary[50];
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3TracingGEOYZ_NoelectronPrimary[50];
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3TracingPadRowvsZ_NoelectronPrimary[50];
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3TracingPadColvsZ_NoelectronPrimary[50];
  std::auto_ptr<TH2F> SelectedEventGeantPadTrkWithMoreThan3TracingGEOXY_NoelectronPrimary[50];


  std::auto_ptr<TH1D> ArmPlus_TrkRecoCutMult_FixedGeantMult[50];
  std::auto_ptr<TH1D> ArmPlus_TrkRecoCutMult_FixedGeantMult2[50];
  std::auto_ptr<TProfile> ArmPlus_TrkEtaEfficiency_FixedGeantMult[50];
 
  std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division[10];
  std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant[10];
  std::auto_ptr<TProfile> H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  std::auto_ptr<TProfile> H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  std::auto_ptr<TProfile> H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
   std::auto_ptr<TProfile> H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_PhiCut[10];
 
  std::auto_ptr<TProfile> H0_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult[10];
  std::auto_ptr<TProfile> H1_TrkPhiEfficiencyCumulative_UnfoldCutsBinMult[10];

  std::auto_ptr<TH2F> PlaneGeantPadColXY[40];
  std::auto_ptr<TH2F> PlaneGeantToRecXY[40];

  std::auto_ptr<TH1D> Associated_PrimaryTrkZImpact_H0[22];
  std::auto_ptr<TH1D> Associated_PrimaryTrkZImpact_H1[22];
  std::auto_ptr<TH1D> Associated_PrimaryTrkZImpact_H2[22];
  std::auto_ptr<TH1D> Associated_PrimaryTrkZImpact_H3[22];
  


  std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative;
  std::auto_ptr<TProfile> H1_TrkEtaEfficiencyCumulative;
  std::auto_ptr<TProfile> H2_TrkEtaEfficiencyCumulative;
  std::auto_ptr<TProfile> H3_TrkEtaEfficiencyCumulative;

  std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative_UnfoldCuts;
  std::auto_ptr<TProfile> H1_TrkEtaEfficiencyCumulative_UnfoldCuts;
  std::auto_ptr<TProfile> H2_TrkEtaEfficiencyCumulative_UnfoldCuts;
  std::auto_ptr<TProfile> H3_TrkEtaEfficiencyCumulative_UnfoldCuts;
  std::auto_ptr<TProfile> H1_TrkPhiEfficiencyCumulative_UnfoldCutsCumul;
  std::auto_ptr<TProfile> H0_TrkPhiEfficiencyCumulative_UnfoldCutsCumul;

  std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative_UnfoldCuts_2;
  std::auto_ptr<TProfile> H1_TrkEtaEfficiencyCumulative_UnfoldCuts_2;
  std::auto_ptr<TProfile> H2_TrkEtaEfficiencyCumulative_UnfoldCuts_2;
  std::auto_ptr<TProfile> H3_TrkEtaEfficiencyCumulative_UnfoldCuts_2;
  std::auto_ptr<TH1D>  NumTrkInT2PerEvt7mcut;
  std::auto_ptr<TH1D>  NumTrkInT2PerEvt;
  std::auto_ptr<TH1D>  NumChPartInT2PerEvt;


  std::auto_ptr<TProfile> ArmPlus_TrkEtaEfficiencyCumulative;
  std::auto_ptr<TProfile> ArmMinus_TrkEtaEfficiencyCumulative;
  std::auto_ptr<TProfile> ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts;
  std::auto_ptr<TProfile> ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts;
  std::auto_ptr<TProfile> ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2;
  std::auto_ptr<TProfile> ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts_2;
  

  std::auto_ptr<TProfile> ArmPlus_TrkEtaEfficiencyCumulativeSecondary;
  std::auto_ptr<TProfile> ArmMinus_TrkEtaEfficiencyCumulativeSecondary;
  std::auto_ptr<TProfile> ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCutsSecondary;
  std::auto_ptr<TProfile> ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCutsSecondary;
  std::auto_ptr<TProfile> ArmPlus_TrkEtaEfficiencyCumulative_UnfoldCuts_2Secondary;
  std::auto_ptr<TProfile> ArmMinus_TrkEtaEfficiencyCumulative_UnfoldCuts_2Secondary;


  

  std::auto_ptr<TProfile> SingleTrkEffiTrigvsE;
  std::auto_ptr<TProfile> DR_vsEta;
  //  T1Geometry * theT1Geometry;
  std::auto_ptr<TH1D> ZDistanceRespToVtz_GEMUpgrade_TB;
  std::auto_ptr<TH1D> YDistanceRespToVtz_GEMUpgrade_TB;
  std::auto_ptr<TH1D> XDistanceRespToVtz_GEMUpgrade_TB;
  std::auto_ptr<TH2D> T1Chi2VsE;
  // std::auto_ptr<TProfile>  Fraction_OfPrimaryTrk_vsEta_UnfoldCuts_2; 
  

  
  //Trk Efficiency as a function of multiplicity, like Halftrk efficiency but using Geant Information.
  //We want to calculate the efficiency in finding trk with zimpact cut satisfied, matching at least 
  //an hit with the primary Geant trk. The Geant Trk should release at least 4 hits:
  //2 in the firsts two plane, two in the last. If you require 4 GeantHit istead of recHit the discrepancy
  //found when comparing with the Simu TrkEffi can be used directly as an error on the nominal primary efficiency
  //that you need to use in the unfolding correction, since this case is the denominator has the Truth count. 
  //The fraction of Geant trk failing
  //to release 4 hits in the ref-quarter is corrected basing on the acceptance correction. You need also to be sure that
  //the ref Geant trk satisfy the Zimpact cut. The fraction of primary  with 4 hit but without ZimpCut should be
  //corrected from acceptance factor as well. Therefore the acceptance factor is the number of geantTrk with satisfied4HitandZimp/#ChPart
  //vs eta. This effi4HitZImpGeant*AccCorrB should be almost equal to effi4Hit*AccCorrA since in the effi4Hit is already cutted for 
  //the Zimp at the recotrk level.

  std::auto_ptr<TProfile> H0TrkEffi_AxisAvgMult_ZImpactCut;
  std::auto_ptr<TProfile> H1TrkEffi_AxisAvgMult_ZImpactCut;
  std::auto_ptr<TProfile> H2TrkEffi_AxisAvgMult_ZImpactCut;
  std::auto_ptr<TProfile> H3TrkEffi_AxisAvgMult_ZImpactCut;

  std::auto_ptr<TProfile> H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut;
  std::auto_ptr<TProfile> H1TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut;
  std::auto_ptr<TProfile> H2TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut;
  std::auto_ptr<TProfile> H3TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut;

  std::auto_ptr<TProfile> H0TrkEffi_AxisAvgMult_NORecCut_GeantZCut;
  std::auto_ptr<TProfile> H1TrkEffi_AxisAvgMult_NORecCut_GeantZCut;
  std::auto_ptr<TProfile> H2TrkEffi_AxisAvgMult_NORecCut_GeantZCut;
  std::auto_ptr<TProfile> H3TrkEffi_AxisAvgMult_NORecCut_GeantZCut;
  
  std::auto_ptr<TProfile> H0TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut_MatchProj;
  std::auto_ptr<TProfile> H1TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut_MatchProj;
  std::auto_ptr<TProfile> H2TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut_MatchProj;
  std::auto_ptr<TProfile> H3TrkEffi_AxisAvgMult_ZImpactCut_GeantZCut_MatchProj;

  std::auto_ptr<TProfile> H0TrkEffi_AxisAvgMult_NORecCut_GeantZCut_MatchProj;
  std::auto_ptr<TProfile> H1TrkEffi_AxisAvgMult_NORecCut_GeantZCut_MatchProj;
  std::auto_ptr<TProfile> H2TrkEffi_AxisAvgMult_NORecCut_GeantZCut_MatchProj;
  std::auto_ptr<TProfile> H3TrkEffi_AxisAvgMult_NORecCut_GeantZCut_MatchProj;
  
  std::auto_ptr<TH1D> NumHitInRecTrkH0;
  std::auto_ptr<TH1D> NumHitInRecTrkH1;
  std::auto_ptr<TH1D> NumHitInRecTrkH2;
  std::auto_ptr<TH1D> NumHitInRecTrkH3;
  std::auto_ptr<TH1D>  AveragePadCLSDistrH1;std::auto_ptr<TH1D>  AveragePadCLSDistrH0;
  std::auto_ptr<TH2D>  corposprimaryXYQ0;
  std::auto_ptr<TH2D>  corposprimaryXYQ1;
  std::auto_ptr<TH2D> Q0PadRCPrimary;
  std::auto_ptr<TH2D> Q1PadRCPrimary;

  std::auto_ptr<TH1D>  Cumulative_GeantTrk_ZImpactRefQ;
  std::auto_ptr<TH1D>    Cumulative_GeantTrk_ZImpactRefQ_EAboveZ5;
  std::auto_ptr<TH1D>    Cumulative_GeantTrk_ZImpactRefQ_PtAboveZ5;
  std::auto_ptr<TH1D>   Cumulative_GeantTrk_ZImpactRefQ_E;
  //  std::auto_ptr<TH1D>  Cumulative_GeantTrkEfficiency;


  std::auto_ptr<TH1D> H1PrimaryTrkPdgID;
  std::auto_ptr<TH1D> H0PrimaryTrkPdgID;
  std::auto_ptr<TH1D> H0PrimaryTrkE;
  std::auto_ptr<TH1D> H1PrimaryTrkE;
  

  std::auto_ptr<TH1D> NumH0TrkWhenGeantEffiFail; 
  std::auto_ptr<TH1D> NumH0PlaneActiveWhenGeantEffiFail; 
  std::auto_ptr<TH1D> H0ZimpDistrWhenGeantEffiFail; 
  std::auto_ptr<TH1D>   minDphiWhenGeantEffiFail;
  std::auto_ptr<TH1D>   minDrWhenGeantEffiFail;
  std::auto_ptr<TH2D>  minDxDyWhenGeantEffiFail;
  std::auto_ptr<TH2D>  minDxDyWhenGeantEffiFailH1;
  std::auto_ptr<TH2D>  Numtrk_vs_NumtrkpasCutsWhenEffiFail;
  std::auto_ptr<TH2D>  NumH0PlaneActive_vsMultipl_WhenGeantEffiFail; 
  std::auto_ptr<TH1D>  ZImpOfARecTrkWhenEffiFail;

  std::auto_ptr<TH1D>  ZImpOfAGeantTrkWhenEffiFailH1;
  std::auto_ptr<TH1D>  ZImpOfAGeantTrkWhenEffiFailH0;
  std::auto_ptr<TH1D>  ZImpOfAPrimaryGeantTrk_Cumulative;
  std::auto_ptr<TH1D>  ZImpOfAPrimaryGeantTrk_CumulativeZMinCut;
  std::auto_ptr<TH1D>  ZImpOfAPrimaryGeantTrk_Cumulative_match;
  std::auto_ptr<TH1D>  ZImpOfAPrimaryGeantTrk_CumulativeZMinCut_match;




  ////////////////////////////////////////////
  std::auto_ptr<TH1D>  Histonumtoeventgen;
  std::auto_ptr<TH1D>  Histonumshowevent;  
  std::auto_ptr<TH1D>  Histonumcleanevent; 
  std::auto_ptr<TH1D> DNDetaT2Clean;
  std::auto_ptr<TH1D> DNDetaT2ALL;
  std::auto_ptr<TH1D> DNDetaT2Show;
  std::auto_ptr<TH1D> NumEventGeneratorTrigger;
  std::auto_ptr<TH1D> NumEventTrackTrigger; 
  
  unsigned int numshowevent;
  unsigned int numtotevent;
  unsigned int numcleanevent;
  //////////////////////////////////////////////
  



  std::auto_ptr<TH1D>  PSimHitT2_EnergyLoss;
  std::auto_ptr<TH1D>  CumulativeGeantTrkHadronR;  
  std::auto_ptr<TH1D>  CumulativeGeantTrkElectronR;
  std::auto_ptr<TH2D>  CumulativeGeantTrkElectronXY;
 
  std::auto_ptr<TH1D> SecondaryVtx_R_ContributingInT2;
  std::auto_ptr<TH1D> SecondaryVtx_Z_ContributingInT2;
  std::auto_ptr<TH1D> Purity_andEfficiency; 
  std::auto_ptr<TH1D>  Purity_andEfficiencyNoGhost;
  std::auto_ptr<TH1D>  RecoTrkEta_MatchingPrimary;
  std::auto_ptr<TH2D> RPhiGeantHitinH0Col32;
  std::auto_ptr<TH2D> XYGeantHitinH0Col32;
  std::auto_ptr<TProfile>   PadRowGeantRecoH1_bis; 
  std::auto_ptr<TProfile>   PadRowGeantRecoH1; 
  std::auto_ptr<TProfile>   PadRowGeantRecoH0_bis; 
  std::auto_ptr<TProfile>   PadRowGeantRecoH0;


  std::auto_ptr<TProfile>   EtaResolGeant;
  std::auto_ptr<TProfile>   EtaResolGenerator;
  std::auto_ptr<TProfile>   EtaResolFreeGeant;
  std::auto_ptr<TProfile>   EtaResolFreeGenerator;



  // ----------member data ---------------------------
  int verbosity;
  std::string TrackLabel2;

  unsigned int numberOfPrimaryTracks;
  unsigned int numberOfSecondaryTracks;

  unsigned int K0EventCounter; unsigned int HadrEventCounter;
  bool fastAnalysis;

  bool PadRoadFinderAnalysis;

  double VtxClusterDistance;
  std::string outputFileName;
  std::string CluLabel;
  std::string HitLabel;
  std::string RoadLabel;
  std::string RoadInstanceLabel;
  std::string TrackLabel;
  int selected_event;
  unsigned int numevent;
  bool UseselectedplanesforAPM;
  int MaxPadCluOfFittingCurves;
  std::string NameOfGenerator;

  unsigned int totnumberbinfordndeta;
  double maxetafordndetahisto; //range assumed simmetric around eta=0
  double etabinsize;
  double etamax;
  double etamin;
  int  numhitRequiredFormMatching;
  bool wantedconditionfordndeta;
  bool fastSimulation;
  double PtCutinPrimaryEfficiency;
  double EnergyCutinPrimaryEfficiency;
  std::vector<int> T2_QuarterUsed;
  std::vector<double> ZEffiCutImpact;

  // std::vector<Local3DPoint> ThePrimaryGeantHitPos;
  
  std::vector<double> etacentreBinV;

  std::vector<double> ZimpCutat95DataVH0; 
  std::vector<double> ZimpCutat95DataLVH0;
  std::vector<double> ZimpCutat95DataRVH0;
 
  std::vector<double> ZimpCutat95DataVH1; 
  std::vector<double> ZimpCutat95DataLVH1;
  std::vector<double> ZimpCutat95DataRVH1;

  std::vector<double> ZimpCutat95DataVH2; 
  std::vector<double> ZimpCutat95DataLVH2;
  std::vector<double> ZimpCutat95DataRVH2;

  std::vector<double> ZimpCutat95DataVH3; 
  std::vector<double> ZimpCutat95DataLVH3;
  std::vector<double> ZimpCutat95DataRVH3;




  std::vector<int> knowPart; std::vector<int> AllowedDecayToCount;
  std::vector<unsigned int> vectorEta2CounterFakePrimary;
  std::vector<unsigned int> vectorEta2CounterDnDetaNoCut;
  std::vector<unsigned int> vectorEta2CounterDnDeta;
  std::vector<unsigned int> vectorEta2SimulatedPrimaryTrkInT2_H0;
  std::vector<unsigned int> vectorEta2SimulatedPrimaryTrkInT2_H1;
  std::vector<unsigned int> vectorEta2CounterDnDetaPrimary;
  std::vector<unsigned int> vectorEta2CounterDnDetaPrimaryNoCut;
  std::vector<unsigned int> vectorEta2CounterRecoDnDeta2UnfoldCut2; 
  std::vector<unsigned int> vectorEta2CounterRecoDnDeta2UnfoldCut2_Primary;
  std::vector<unsigned int> vectorEta2CounterRecoDnDeta2UnfoldCut; 
  std::vector<unsigned int> vectorEta2CounterRecoDnDeta2UnfoldCut_Primary;
  std::vector<unsigned int> vectorMBALLMCGenerator;
  std::vector<unsigned int> vectorEta2GeantPrimaryTrkInH0_ZCut;
  std::vector<unsigned int> vectorEta2PrimaryGeantDNDeta2_IfOneRecoZImpInH0;
};

//DEFINE_FWK_MODULE(T2BackgroundAn);
#endif
