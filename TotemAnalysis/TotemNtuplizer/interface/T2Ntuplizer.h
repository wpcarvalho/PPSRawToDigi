/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Mirko Berretti
*
****************************************************************************/

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "TotemAnalysis/TotemNtuplizer/interface/Ntuplizer.h"
//#include "TotemAnalysis/T2MakeNtples/interface/T2Event.h"
#include "TotemAnalysis/TotemNtuplizer/interface/T2Event.h"
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
  double EtaFromAveragesVtxCorr(T1T2Track thetrk,double vtx_x,double vtx_y,double vtx_z);
  unsigned int RawtoSymb(uint32_t thedet);
  double PartCharge(int pdgcode);

  T2Event t2obj;

private:
  edm::InputTag t2PadDigiCollectionLabel;
  edm::InputTag t2StripDigiCollectionLabel;
  edm::InputTag HepMCProductLabel;
  edm::InputTag t2PadClusterCollectionLabel;
  edm::InputTag t2StripClusterCollectionLabel;
  edm::InputTag HitLabel;
  edm::InputTag RoadLabel;
  edm::InputTag TrackLabel;

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
  double vtxpos;
  double energy;
  double DZScale;
  double tracketamin;
  double tracketamax;
  double PhiChiProbCut;
  double RChiProbCut;
  //bool onlyMultiParticle;
   T2SelectionCutUtils T2CutsUtil;
    
  std::string RawDataName;
  std::string outputFileName;

};
