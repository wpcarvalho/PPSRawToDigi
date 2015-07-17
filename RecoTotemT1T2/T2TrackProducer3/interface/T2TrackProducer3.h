/*
Class T2TrackProducer3
Author: Mirko Berretti
mirko.berretti@gmail.com

Fit best track, also in XY plane

*/


#ifndef _T2TrackProducer3_h
#define _T2TrackProducer3_h

// system include files
#include <memory>
#include <math.h>
#include <TMath.h>
// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Road/interface/T2Road.h"
#include "DataFormats/T1T2Track/interface/T1T2Track.h"
#include "TotemAnalysis/T2Cuts/interface/T2SelectionCutUtils.h"
 




class T2TrackProducer3 : public edm::EDProducer {
public:
  /// Constructor
  T2TrackProducer3(const edm::ParameterSet&);
  virtual ~T2TrackProducer3();

  virtual void produce(edm::Event& event, const edm::EventSetup& setup);

  //virtual T1T2Track fitTrackspecialRZ(T1T2Track, bool  CheckWorstHit=false);

  //T1T2Track fitTrackspecialXY(std::vector<T2Hit> &hitvec);

   T1T2Track MyLinearfitCorr(std::vector<T2Hit> hitvec2,TMatrixD &par_covariance_matrix,double &chi2_,bool Usestrip,std::vector<double> trkparam, int RoadID);
   
   T1T2Track MyLinearfitCorrDEV(std::vector<T2Hit> hitvec2,TMatrixD &par_covariance_matrix,double &chi2_,bool Usestrip,std::vector<double> trkparam, int RoadID,T2Hit &worsthit);

   T1T2Track RPhiFit(std::vector<T2Hit> &hitvec2);

   std::vector<T2Hit> OutliersRemoving(T1T2Track &firstTrkRaw, T2Hit &worsthit);
   
   std::vector<T2Hit> VectOutliersRemoving(std::vector<T2Hit> &hitv, T2Hit &worsthit);

  private:
   
   T2SelectionCutUtils T2CutsUtil;

   bool RemoveOutliers;   
   bool PickUpDisplacedHit;
   double PickUpRadius;
   bool StripFitting;
   bool verbosity;
   bool forceRZfit;
   bool forceXYfit;
   double theChiRidThr;
   unsigned int MinHitFinalTrk;
   unsigned int MinHitFinalTrkCl1Hit;
   unsigned int   MaxHitsInRoad;
   std::string RoadLabel;
   std::vector<double> RecoHitRErrorPerStripcount;
   bool UseRefittingProcedure;
   double DropWorstRecHitChisquareRThreshold;
   double VtxPositionEX,VtxPositionEY,VtxPositionEZ,VtxPositionX,VtxPositionY,VtxPositionZ;
   bool FitVertexPosition;
   string RoadModuleLabel; string RoadInstanceLabel;
   bool GhostSuppression;
   std::string theT2Hit_to_Track_Map0_Module;
   std::string theT2Hit_to_Track_Map0_Instance;
   std::string TrackInstanceLabel;
   std::string HitLabel;

};
#endif
