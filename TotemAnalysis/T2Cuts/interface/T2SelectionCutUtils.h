#ifndef T2SelectionCutUtils_h
#define T2SelectionCutUtils_h


#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"

 


class T2SelectionCutUtils {

 private:
  int i;

 public:
  
   
  typedef struct{
    long double x;
    long double y;
    long double z;
  } point3d ;
  
  typedef struct{
    std::vector<T1T2Track> VtxTrks;
    point3d VtxPos;
    point3d VtxPosRZ;
    point3d VtxPosXY;
  } CandidateVtxWithTrks;

  typedef struct{
    long double distance;
    point3d midPoint;
  } TracksSeparation;


  T2SelectionCutUtils();
  virtual ~T2SelectionCutUtils();
  
  double FirstPlaneForThishit(T2Hit hit);

  void SetCuts(double T2_TrkEtamin_,double T2_TrkEtaMAX_,int T2_trkMultimin_,int T2_trkMultiMAX_,double T2_DZMultiplicator_,double T2_PhiChiProbCut_,double T2_RChiProbCut_,std::vector<int> T2_QuarterUsed_, bool XYFitUsed_);

  void SetCuts(double T2_TrkEtamin_,double T2_TrkEtaMAX_,int T2_trkMultimin_,int T2_trkMultiMAX_,double T2_DZMultiplicator_,double T2_PhiChiProbCut_,double T2_RChiProbCut_,std::vector<int> T2_QuarterUsed_,double IgnoredSmallAnglePrimarySlopeCut_,bool XYFitUsed_);

  bool DZCutPassed(T1T2Track trk);
  bool ThetaCutPassed(T1T2Track trk);
  bool AcceptThisT2Track(T1T2Track trk);
  bool AcceptThisT2TrackWhatEverQuarter(T1T2Track trk);
  bool T2TrkMultCond(unsigned int trkmultiplicity_);
  bool IsinT2(double eta);
  double sigmaZfunctE10(double eta);
  bool ChiCutCond(T1T2Track trk);
  bool EtaTrkCompatibleWithImpact(T1T2Track trk,double cut);
  double PhiFitAverage(std::vector<T2Hit> hitvec);

  bool ChiCutCond(T1T2Track trk, bool xyfitused,double prob1, double prob2);
  bool TrkAlsoInQuarter(T1T2Track trk);
  bool TrkAlsoInQuarter(T1T2Track trk, std::vector<int>  QuarterWanted);

  bool TrkInQuarter(T1T2Track trk);
  bool TrkInQuarter(T1T2Track trk,unsigned int selectedquarter);

  bool TrkInQuarterTight(T1T2Track trk);
  bool TrkInQuarterTight(T1T2Track trk,unsigned int selectedquarter);

  bool HitVectorInAQuarter(std::vector<T2Hit> hitvect, unsigned int selectedquarter);


  bool ThetaAsPrimary(T1T2Track trk, bool usexytrack); 
  std::vector<T2Hit> HitsFromTrk(T1T2Track trk);
  std::vector<T2Hit> HitsFromTrkInOneQuarter(T1T2Track trk, unsigned int selectedquarter);
  std::vector<T2Hit> HitsFromTrkInOneQuarter(std::vector<T2Hit> vv, unsigned int selectedquarter);
  TracksSeparation TwoTracksSeparation(T1T2Track t1, T1T2Track t2, bool usesXYtracking);
  TracksSeparation TrkZ0(T1T2Track t1, bool usesXYtracking);

  //void FindPrimaryVtx(std::vector<T1T2Track> stratingtracks, bool usexytrack, bool UseAnyQuarter=false);


  void FindPrimaryVtx(std::vector<T1T2Track> startingtracks2,bool usexytrack, bool UseAnyQuarter=false, double Maxbxtrk=100.,double Maxbytrk=100.);


  std::vector<T2SelectionCutUtils::CandidateVtxWithTrks> FindAllVtxs(std::vector<T1T2Track> startingtracks2, bool usexytrack, std::vector<int> Quartertouse, double ClusteringDistance_mm);



  //void FindPrimaryVtx(std::vector<T1T2Track> stratingtracks, bool usexytrack, bool UseAnyQuarter, double trkmaxbx, double trkmaxby);
  //void T2SelectionCutUtils::FindPrimaryVtx(std::vector<T1T2Track> startingtracks2,bool usexytrack, bool UseAnyQuarter=false, double Maxbxtrk=100.,double Maxbytrk=100.){ 

  //Vector size=2: Return XY and RZ Vtx position for an event
  std::vector<point3d> GivePrimaryVtx(std::vector<T1T2Track> startingtracks2, bool usesXYtracking, bool UseAnyQuarter=false);


  //An XY Fit is assumed, so you need precise tracks
  double FindAverageTrackZ(std::vector<std::vector<T2Hit> >*theTracks, char XorY);

  
  //Vector size=3: Return the average X-Y-Z of a vertex collection
  std::vector<double> vtxPosToRet(std::vector<point3d> allvtxs); 


  long double CartesianDistance(point3d p1,point3d p2);
  /*
  double GetInternalMisErrorDX();
  double GetInternalMisErrorDY();
  double SetInternalMisErrorDX(double dxmisalError_);
  double SetInternalMisErrorDY(double dxmisalError_);
  */
  T1T2Track fitTrackspecialXY(std::vector<T2Hit> vect, bool worsthitcheck);
  int findWorstPoint(const T1T2Track &trk);

  
  T1T2Track RPhiFit(std::vector<T2Hit> hitvec);
  double EtaFromAverages(T1T2Track thetrk);
  double EtaFromAveragesVtxCorr(T1T2Track thetrk,double vtx_x,double vtx_y,double vtx_z);
  double Chi2R_Prob_Calculator(T1T2Track thetrk); //Take a track (whatever fit), and return the chi2-RZ fit probability 

  T1T2Track TrackFromHits(bool DoXYTrkFit, std::vector<T2Hit> hitvec);
  std::vector<double> MyLinearfitCorr(std::vector<T2Hit> hitvec,TMatrixD &par_covariance_matrix,double &chi2_);

  std::vector<double> MyLinearfitXYSimple(std::vector<T2Hit> hitvec);
  std::vector<T2Hit> GiveMeTheHitForResidual(unsigned int planeIdToExclude, std::vector<T2Hit>  hitv);
  T1T2Track TrackXYFromHitsExcludingOneplane(unsigned int planeIdToExclude,  T1T2Track origtrk);

  std::vector<T2Hit> GetHitsInQuarterFrame(std::vector<T2Hit> hitvec);

  CandidateVtxWithTrks storeVtx;

  bool XYFitUsed;
  
  double dxMisalError;
  double dyMisalError;

  double T2_TrkEtamin;
  double T2_TrkEtaMAX;
  int T2_trkMultimin;
  int T2_trkMultiMAX;
  double T2_DZMultiplicator;
  double T2_PhiChiProbCut;
  double T2_RChiProbCut;
  double IgnoredSmallAnglePrimarySlopeCut;
  std::vector<int> T2_QuarterUsed;
  

  //    static double Zposition(uint32_t cmsswid);
  
};

#endif




