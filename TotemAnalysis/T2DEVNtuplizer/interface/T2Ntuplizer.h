/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Mirko Berretti
*
****************************************************************************/


// system include files
/*
#include <memory>
#include <math.h>
#include <TMath.h>
*/
// user include files
/*
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/T2Cluster/interface/T2Cluster.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2StripClusterCollection.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/T2Road/interface/T2Road.h"
#include "DataFormats/T2Road/interface/T2RoadCollection.h"
#include "DataFormats/T1T2Track/interface/T1T2Track.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
#include "TotemAnalysis/T2Cuts/interface/T2SelectionCutUtils.h"
*/
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "TotemAnalysis/T2DEVNtuplizer/interface/Ntuplizer.h"
#include "TotemAnalysis/T2DEVNtuplizer/interface/T2EventDEV.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "Geometry/TotemGeometry/interface/T2Geometry.h"
#include <boost/shared_ptr.hpp>
//#include "TotemAnalysisT2MakeNtplesDEVELT2NtplizerDEVLinkDef.h"
//#include "TotemAnalysis/T2MakeNtplesDEVEL/src/LinkDef.h"
//#include "TotemAnalysis/T2MakeNtplesDEVEL/interface/T2EventDEV.h"
#include "TTree.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/IO_GenEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"




#include "TotemAnalysis/T2Cuts/interface/T2SelectionCutUtils.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"

#include <memory>

/**
 *\brief Saves T2 reconstruction data into the common ntuple.
 **/
class T2Ntuplizer : public Ntuplizer
{
public:
  T2Ntuplizer(const edm::ParameterSet&);
  virtual ~T2Ntuplizer() {}
  
  virtual void CreateBranches(const edm::EventSetup&, TTree *);
  virtual void FillEvent(const edm::Event&, const edm::EventSetup&);

  T2EventDEV t2obj;

private:

  edm::InputTag simTrackContainerLabel;
  edm::InputTag simVertexContainerLabel;
  edm::InputTag t2PadDigiCollectionLabel;
  edm::InputTag t2StripDigiCollectionLabel;
  double PartCharge(int pdgcode);

unsigned int loadSimulatedTracks(const edm::Event&,
				   std::auto_ptr<edm::SimTrackContainer>&
				   );


 double GetEta2FromGeantHits(std::vector<std::pair<double,double> > SimHitRZ_ofGeantPrimaryTrk,double vtx_z);







   
  unsigned int loadSimulatedVertices(const edm::Event&,
				     std::auto_ptr<edm::SimVertexContainer>&
				     );

  unsigned int loadSimulatedHits(const edm::Event&,
                                 const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
                                 const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
                                 std::auto_ptr<edm::PSimHitContainer>& theSimHits
                                 );

  Local3DPoint CalculateCorrectXYZPos(const Local3DPoint &pos, const unsigned int detUnitId);
   
  double EtaFromAveragesVtxCorr(T1T2Track thetrk,double vtx_x,double vtx_y,double vtx_z);

  

  void GeantTrkAnalyzer(const std::auto_ptr<edm::PSimHitContainer>& theSimHits,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt, std::vector<double> &PrimGeantTrkEta2InT2, std::vector<int> &PrimGeantTrkHemiInT2,std::vector<double> &PrimGeantTrkPTInT2,std::vector<int> &PrimGeantTrkBarcodeInT2 ,double vtxZ);



  void TracksInfoAssociationDEBUG(const edm::Event& event,
				  const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
				  const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
				  const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,				 
				  const HepMC::GenEvent* evt,
				  T1T2Track recTrack,
				  std::vector<double> &PrimaryTrkenergy_ForTheRecoHits,
				  std::vector<int> &ThePrimTrkIdFromHits,std::vector<int> &ThePrimTrkIdFromHitsReal,std::vector<int> &TheSecondTrkIdFromHits,std::vector<int> &TheSecondTrkIdFromHitsReal
				  );

  std::vector<int>  RecotracksInfoFromSimu_BIS(const edm::Event& event,
					       const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
					       const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
					       const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,				 
					       const HepMC::GenEvent* evt,
					       T1T2Track recTrack,
					       std::vector<double> &PrimaryTrkenergy_ForTheRecoHits,
					       std::vector<double> &PrimaryTrkpt_ForTheRecoHits,
					       std::vector<int> &ThePrimTrkIdFromHits,std::vector<int> &ThePrimTrkIdFromHitsReal,std::vector<int> &TheSecondTrkIdFromHits,std::vector<int> &TheSecondTrkIdFromHitsReal
					       );
  
  int CaluclateCommonPDGSECMoth(std::vector<int> &TheSecondaryTrkIdFromRecHits,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,bool printing,int& motherbarc);
  int CaluclateCommonPDGSEC(std::vector<int> &TheSecondaryTrkIdPDGID);
  int GetTrkOlderMotherPid(SimTrack atrk,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
			   const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,int &barcodeMother,int &DirectTrkPID, bool printing);
  

  unsigned int RawtoSymb(uint32_t thedet);
  /*
    Returns track corresponding to trackId
  */
  bool getTrack(const /*unsigned*/ int trackId,
                const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
                SimTrack &track);

  /*
    Returns vertex corresponding to vertexId
  */
  bool getVertex(const /*unsigned*/ int vertexId,
                 const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
                 SimVertex &thevertex);
  
  /*
    Checks if given TRACKID corrensponds to the primary vertex.
    (Not a secondary track)
    reurns true if TRACKID corrensponds to the primary vertex otherwise false
    also returns false if no vertex found
  */
  bool isTrackPrimary(const /*unsigned*/ int trackId,
                      const std::auto_ptr<edm::SimTrackContainer>&,
                      const std::auto_ptr<edm::SimVertexContainer>&);
  
  bool PrimaryTrackCondition(SimTrack atrk,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
			     const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,int &barcodeMother);

  double CalculateGeneratorEtaFromTrkId(std::vector<int> &ThePrimTrkIdFromHits,const HepMC::GenEvent* evt);
  

  T1T2Track CalculateRecoTrkFromGeant(std::vector<Local3DPoint> PrimaryGeantHitPos, std::vector<double> &theZ,  std::vector<uint32_t> &thecmsswId);
  bool CalculateGeantEtaFromTrkId(std::vector<int> &ThePrimTrkIdFromHitsReal,const std::auto_ptr<edm::PSimHitContainer>& theSimHits,const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,T1T2Track &TheGeantTrkReconstr);
 

  double MostProbableTrkEnergy(std::vector<double> &PrimaryTrkenergy_ForTheRecoHits);  
  //double MostProbableTrkPt(std::vector<double> &PrimaryTrkpt_ForTheRecoHits);  
  T1T2Track MyLinearfitCorrDEV(std::vector<T2Hit> hitvec2,TMatrixD &par_covariance_matrix,double &chi2_,bool Usestrip, int RoadID);

  T1T2Track RPhiFit(std::vector<T2Hit> &hitvec2);
  
  // virtual void beginJob(const edm::EventSetup&) ;
  //virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void endJob() ;

private:
  std::auto_ptr<TH1F> NumhitinTrack;
  std::auto_ptr<TH1F> NumhitinTrackGood;

  std::auto_ptr<TH1F> clusterstripentries;
  std::auto_ptr<TH1F> clusterpadentries;
  std::auto_ptr<TH1F> diffphiCluGun;
  std::auto_ptr<TH1F> RLocHRecoH;
  std::auto_ptr<TH1F> Trketa;
  std::auto_ptr<TH1F> Trkphi;
  std::auto_ptr<TH1F> DPhiGoodTrk;
  std::auto_ptr<TH1F> DEtaGoodTrk;
  std::auto_ptr<TH1F> diffphiGUNHIT;
  std::auto_ptr<TH1F> diffphiCluHit;
  std::auto_ptr<TH1F> Trkphigood;
  std::auto_ptr<TH1F> Trketagood;
  std::auto_ptr<TH1F> diffRCluHit;
  std::auto_ptr<TH1F> Chi2RProb;
  std::auto_ptr<TH1F> Chi2PhiProb;
  std::auto_ptr<TH1F> TrKPhi0degstudyg;
  std::auto_ptr<TH1F> SourceofSecondary;
  std::auto_ptr<TH1F> SourceofReco;
  std::auto_ptr<TH1F> Z0distro;
  std::auto_ptr<TH1F> R0distro;
  std::auto_ptr<TH1F> tantheta;
  std::auto_ptr<TH1F> tanthetam0;
  std::auto_ptr<TH1F> RecoZHit;
  std::auto_ptr<TH1F> SimuZHit;


  std::auto_ptr<TProfile>  MeanT2PadvsEta; 
  std::auto_ptr<TProfile>  MeanT2StripvsEta;
  std::auto_ptr<TProfile>  MeanT2RoadEntriesvsEta;
  std::auto_ptr<TProfile>  MeanT2HitCl1vsEta;
  std::auto_ptr<TProfile>  SingleParticleEfficiency;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyAllCuts;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyCutEta;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyCutChi;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyCutDZ;
  std::auto_ptr<TProfile>  MultiParticleEfficiencyCut;
  std::auto_ptr<TProfile>  MultiParticleEfficiencyCutNorm;
  std::auto_ptr<TProfile>  MeanT2HitvsEta;
  std::auto_ptr<TProfile>  MeanT2GeantHitvsEta;
  std::auto_ptr<TProfile>  P0NumEtacorr;
  std::auto_ptr<TProfile>  P0NumEnergycorr;

   std::auto_ptr<TH1F> TrkEtaGoodH0;
   std::auto_ptr<TH1F> TrkEtaGoodH1;
     std::auto_ptr<TH1F> TrkEtaH0;
   std::auto_ptr<TH1F> TrkEtaH1;

  // std::auto_ptr<TProfile>   SingleParticleEfficiencyAllCutsH0; 
  // std::auto_ptr<TProfile>   SingleParticleEfficiencyAllCutsH1;
  std::auto_ptr<TH1F> EtaParticle;  
  std::auto_ptr<TH1F> EtaParticleAll;  
  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem2; 
  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem1; 
  std::auto_ptr<TProfile>  MeanEtaResvsPhi;

  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem2Part; 
  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem1Part; 
  std::auto_ptr<TProfile>  MeanEtaResvsPhiPart;
  std::auto_ptr<TProfile>  P0EtaEnergycorr;
  std::auto_ptr<TProfile>  MeanEtaResvsPhiCut;
  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem2Cut; 
  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem1Cut;
  std::auto_ptr<TProfile>   MeanT2StripDigivsEta;
  std::auto_ptr<TProfile>   MeanT2PadDigivsEta;

  std::auto_ptr<TProfile>  NumChPartVsNumP0; 
  std::auto_ptr<TProfile>  NumP0VsNumTracks;
  std::auto_ptr<TProfile>  NumChPartInVsNumChPartOut;

  std::auto_ptr<TH1F> HNumrecotrackcutmag0;
  std::auto_ptr<TH1F> HTrackcounter;
  std::auto_ptr<TH1F> chPartinT2;

  std::auto_ptr<TH2F>  P0EtaEnergy;
  std::auto_ptr<TH2F>  muEtaEnergy;
  std::auto_ptr<TH2F>  ChOutEtaEnergy;
  std::auto_ptr<TH1F>  ChEnergyinT2;
  std::auto_ptr<TH1F> HTrackcounterNonCrit;
  std::auto_ptr<TH1F> HTrackcounterCrit;
  std::auto_ptr<TH1F> stablepdg;
  std::auto_ptr<TCanvas> Chi2PhiProbLogy;
  std::auto_ptr<TCanvas> Chi2RProbLogy;

  std::auto_ptr<TH2F>  PimenoEtaEnergy[15];

  std::auto_ptr<TH1F> chiPhidistro;
  std::auto_ptr<TH1F> chiRdistro;
 
  std::auto_ptr<TH1F> SymbIdHT0;
  std::auto_ptr<TH1F> Trketachicut;
  std::auto_ptr<TH1F> Trketaetachicut; 
  std::auto_ptr<TH1F> Trketaetacut;
  std::auto_ptr<TH1F> Trketadzcut;
  std::auto_ptr<TH1F> Trketadzetacut;
  std::auto_ptr<TH1F> Trketadzetachicut;
  std::auto_ptr<TH1F> DPhiChiCutOneTrk;
  std::auto_ptr<TH1F> DEtaChiCutOneTrk;
  std::auto_ptr<TH1F> MinDz2Plane;
  std::auto_ptr<TH1F>  energyP04575;
  std::auto_ptr<TH1F>  etaP04575;
  bool singleparticle;

  unsigned int numevent;
  double Chicut;
  double energy;
  double DZScale;
  double tracketamin;
  double tracketamax;
  double PhiChiProbCut;
  double RChiProbCut;
  //bool onlyMultiParticle;
   T2SelectionCutUtils T2CutsUtil;
  
  std::vector<int> knowPart; std::vector<int> AllowedDecayToCount;

  bool IncludeT1;
  std::string RawDataName;
  std::string outputFileName, HepMCProductLabel, CluLabel, HitLabel, RoadLabel,TrackLabel;
  bool LoadPrimaryGeantInformation;
  bool EasyGeneratorPrimaryDef;






};

