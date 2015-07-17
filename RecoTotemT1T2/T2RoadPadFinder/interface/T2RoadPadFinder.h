/**
 * Class T2RoadPadFinder.
 *
 * Author: Mirko Berretti / University of Siena 
 * Email:  mirko.berretti@gmail.com
 * Date:   2010-12-08
 */

#ifndef _T2RoadPadFinder_h
#define _T2RoadPadFinder_h

#include <map>
#include <cstdlib>

#include <boost/shared_ptr.hpp>

#include "TMath.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Road/interface/T2Road.h"
#include "RecoTotemT1T2/T2RoadPadFinder/interface/Tube.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include "DataFormats/T2Cluster/interface/T2Cluster.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2StripClusterCollection.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2Hit/interface/T2HitCollectionMapping.h"
#include "DataFormats/T2Hit/interface/T2PadStripAssociator.h"
#include "DataFormats/T2Cluster/interface/cluster_entry.h"
#include "DataFormats/T1T2Track/interface/T1T2Track.h"



class T2RoadPadFinder : public edm::EDProducer {
public:
  /// Constructor
  T2RoadPadFinder(const edm::ParameterSet &config);


  /// Destructor
  virtual ~T2RoadPadFinder();

  /// The method which produces the Roads
  virtual void produce(edm::Event& event, const edm::EventSetup& setup);

  

  //Load a maps containing <ClustuniqueId->CluRoadId> and <ClustuniqueId->CluStatus>  
  //The map will be filled for the quarters eventually included in QuartedSelected
  //In the meantime also a single map of <ClustuniqueId->index position inside the t2padclcoll will filled>
  void InitializeMapClusterAssociators(const edm::Event&);

  void PropagateRoads(unsigned int seedplanelevel,unsigned int matchplanelevel,std::vector<Tube> &existingRoad); 

  Tube  MatchConditonAndTubePropagation(Tube atube,int cluuniqueid,const T2PadClusterCollection* T2PadClusterCollectionPtr,double &distance,bool &ismatching);
  
  void MakethreeWithThisPlaneSeed(unsigned int seedplane,std::vector<Tube> &AllInitialTubes, std::vector<Tube> &AllFinalTubes,const T2PadClusterCollection *PtrToPadColl);
  
  Tube MergeSmallTubes(Tube &Tube1,Tube &Tube2);

  
  bool SmallAngleTube(int candidateCluIdStart, int candidateCluIdStop, const T2PadClusterCollection* T2PadClusterCollectionPtr);

  bool QJumpCheckForTubeGen(int candidateCluIdStart, int candidateCluIdStop,const T2PadClusterCollection *PtrToPadColl);

  Tube New_TwoPoints_Tube(int candidateCluIdStart, int candidateCluIdStop,bool isbackw,const T2PadClusterCollection *PtrToPadColl);

  bool CheckMatchingDistance(std::vector<double> &propagatedP,int cluuniqueid, T2Cluster &clusterchecked, std::vector<int> &tubeCluIds,const T2PadClusterCollection *T2PadClusterCollectionPtr,double &distance);

  std::vector<double> TubePredictedPoint(Tube &atube, double planeZ, int &countgoodhit);
  
  T2Cluster PadExtractionFromBlob(T2Cluster &blobclu, std::vector<double> &propagatedP, bool &allOK);

  bool QuantumJumpCheck(T2Cluster &clusterchecked, int cluuniqueid, std::vector<int> &tubeCluIds,const T2PadClusterCollection *T2PadClusterCollectionPtr);


  bool ZombieForwCompatibility(Tube &TubeForw,Tube &TubeZombie,const T2PadClusterCollection *T2PadClusterCollectionPtr);


  //std::vector<std::vector <long> > RoadPadCL1HitsCombination(std::map<double, std::vector<long> > &ZvsUniqueHitId_inroad);
     
  std::vector<std::vector <unsigned long int> > RoadPadCL1HitsCombination(std::map<double, std::vector<unsigned long int> > &ZvsUniqueHitId_inroad);


  int GiveMeBestCandidate(std::vector<std::vector<T2Hit> > &allcandidateRoads, double &chi2X, double &chi2Y,std::vector<std::vector<double> > &BxByAxAy_AllRoads);


  void RemoveOutliers(std::vector<std::vector<T2Hit> > &analyzedRoadS,std::vector<std::vector<double> > &BxByAxAy_AllRoads);
    
  int findWorstPoint(std::vector<T2Hit> &trk,std::vector<double>  &BxByAxAy);

  T1T2Track fitTrackspecialXY(std::vector<T2Hit>  &analyzedRoad);//;bool UseAnyQuarter=false);
 
  T1T2Track fitTrackspecialXYDEV(std::vector<T2Hit>  &analyzedRoad);

  T1T2Track TrackerFitter(std::vector<T2Hit> &hitvec2);

  void RemoveHitsFromLoserCombinations(std::vector<std::vector<T2Hit> > &allcandidateRoadsCL1,std::vector<T2Hit> &winnigRoad); 

  std::vector<std::vector<T2Hit> > ConcurrentBranchResolver(std::vector<std::vector<T2Hit> > &ParallelBranches);

  bool TrkInOverlapRegion(T2Road &theRoad);
  
  bool SameTrksInOverlap(std::vector<T2Hit> &Road1,std::vector<T2Hit> &Road2,double &Deta, double &dphi, double &DR);
  
  
  Tube MakeStraightTubeForThisclusters(std::vector<int> &cluIds,const T2PadClusterCollection*  T2PadClusterCollectionPtr);

  void StraightTubeFinder(int quarter, std::vector<Tube> &WinningTubes,const T2PadClusterCollection*  T2PadClusterCollectionPtr);

  double GetClusterRadialError(T2Cluster &clu);
  
  template <typename T> void FreeVectorMemory(T & t);

 private:

  T2GeometryUtil convGeo;
  T2GeometryUtil::T2DetInfo planeinfo;
  
  std::vector<int> QuartedSelected;
  std::vector<int> VplaneToExclude;
  int MinCluSize_considered_asBlobs;
  double TolleranceThetaX;
  double TolleranceThetaY;
  double TolleranceDX;
  double TolleranceDY;
  int InefficiencyMaxJump;
  bool AllowsPadReAssociation, AllowsConcurrentBranches; 
  int Nmin_padsFinal;            
  string HitLabel;
  string CluLabel;
  int MinimumNumCl1Hit;
  double chi2XYProb_Thr;

  string T2RoadCollProdName;
  string HitToRoadAssociatorName;
  unsigned int TubeGeneratorId;
  double TwoPointsTubesAngularCollinearity;
  int verbosity;

  bool ResolveOverlapDoubleCount;
  double OverlapDoubleCountDR;
  double OverlapDoubleCountDPhi;
  double OverlapDoubleCountDTheta;   
  double BiggestTubeAngleConsideredINPUT;
  
  double BiggestTubeAngleConsidered;
  double NumSigma;

  double BiggestTubeAngleConsideredH0;
  double BiggestTubeAngleConsideredH1;
  double BiggestTubeAngleConsideredH2;
  double BiggestTubeAngleConsideredH3;
  
  bool useStraightPadTowers; 
  double NumPadCluOccupancyAlert;

  double occupadcluH0;
  double occupadcluH1;
  double occupadcluH2;
  double occupadcluH3;

  std::vector<Tube> WinningTubes;
  // std::auto_ptr<T2PadClusterCollection> T2PadClusterCollectionPtr;
  const T2PadClusterCollection *T2PadClusterCollectionPtr;


  edm::Handle< T2PadClusterCollection > t2padclcoll; //I want to see this collection from all the procedures

  edm::Handle< T2StripClusterCollection > t2strclcoll;
  
  edm::Handle<T2HitCollection> t2hitcoll;

  edm::Handle<T2HitCollectionMapping> t2hitcollmap;

  edm::Handle<T2PadStripAssociator> t2hitpadonstripMap;

  //std::vector<int> H0CluPadIdsMatrix[24][65];
  //std::vector<int> H1CluPadIdsMatrix[24][65];
  //std::vector<int> H2CluPadIdsMatrix[24][65];
  //std::vector<int> H3CluPadIdsMatrix[24][65];
  //unsigned int LookUpCluPadIdsMatrix[25][64]; 

  
  std::vector< std::vector< std::vector<int> > > H0CluPadIdsMatrix;
  std::vector< std::vector< std::vector<int> > > H1CluPadIdsMatrix;
  std::vector< std::vector< std::vector<int> > > H2CluPadIdsMatrix;
  std::vector< std::vector< std::vector<int> > > H3CluPadIdsMatrix;
  


  std::map<int,std::vector<int> > Map_ClustuniqueId_VectorCluRoadId; 
  std::map<int, int> Map_ClustuniqueId_CluRoadId; 
  std::map<int, int> Map_ClustuniqueId_CluStatus; 
  std::map<int, std::pair <uint32_t,unsigned int> > Map_ClustuniqueId_CluIndexPos; 
  std::map<int, double> Map_ClustuniqueId_ZPos; 
  std::map<unsigned int, std::vector<int> > Map_PlaneSymb_ClustersuniqueId;
  std::map<int,unsigned int> Map_ClustuniqueId_PlaneSymb;
  std::map<int,bool> Map_ClustuniqueId_ISBLOB;
  std::map<int,int> Map_RealClusterID_ToSymbUsed;
  std::map<int,unsigned int> Map_ConcurrentTubeId; 
  std::map<int,std::vector<unsigned int> > Map_cloneRoadID_ToRoadPosIndex;
};
#endif

//For seedplanelevel=0..10 //Piano da dove parte la costr del seed   
// For matchplanelevel=seedplanelevel..10
  //CreateRoads(seedplanelevel,matchplanelevel,&existingRoad)-> creo la road da seedplane in avanti 
    
