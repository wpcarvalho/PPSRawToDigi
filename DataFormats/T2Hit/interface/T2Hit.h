/*
  Class T2Hit
  Author: Mirko Berretti
  mirko.berretti@gmail.com
*/

#ifndef T2HIT__
#define T2HIT__


#include <stdint.h>
#include <vector>
#include <iostream>
#include "Rtypes.h"

#include "DataFormats/T2Hit/interface/hit_entry.h"    //a class that contain data of the clusters forming hit
#include "DataFormats/T2Cluster/interface/cluster_entry.h"


using namespace std;


class T2Hit
{
 public:

  T2Hit();
  T2Hit(double expX, double expY, double Z, uint32_t cmsswid);
  T2Hit(double X, double Y, double Z, double EX, double EY, double EZ, uint32_t cmsswid); //Used for trk imposing vtx.
  
  virtual ~T2Hit();

  typedef vector<hit_entry> vecthit_entries;
  //typedef map<int, map<int,int> > CluChMap;
   
  inline float GetHitR() const {return hitR_;}
  inline float GetHitPhi() const {return hitPhi_;}
  inline float GetHitZ() const {return hitZ_;}
  inline float GetHitDR() const {return hitDR_;}
  inline float GetHitDPhi() const {return hitDPhi_;}
  inline float GetHitDZ() const {return hitDZ_;}
    
  inline float GetHitX() const {return hitX_;}
  inline float GetHitDX() const {return hitDX_;}
  inline float GetHitY() const {return hitY_;}
  inline float GetHitDY() const {return hitDY_;}
	
  //CluChMap hitClusterPadVect;
  //CluChMap hitClusterStripVect;
  // inline CluChMap GetHitClusterPadVect(){return hitClusterPadVect;}
  // inline CluChMap GetHitClusterStripVect(){return hitClusterStripVect;}

  inline unsigned int GetHitClass() const {return hitClass_;}
  inline unsigned int GetHitNumStrip() const {return hit_numstrip_;}    
  inline unsigned int GetHitNumPad() const {return hit_numpad_;}
  inline unsigned int GetHitPadCol() const {return hit_numpadcol_;}


  inline unsigned int GetHitArm() const {return hitarm_;}
  inline unsigned int GetHitHalftele() const {return hithalftele_;}
  inline unsigned int GetHitPlane() const {return hitplane_;}
  inline unsigned int GetHitPlaneSide() const {return hitplaneside_;}
  inline uint32_t GetHitDetRawId() const {return detrawid_;}

  inline void SetHitArm(unsigned int harm){hitarm_=harm;}
  inline void SetHitHalftele(unsigned int hht){hithalftele_=hht;}
  inline void SetHitPlane(unsigned int hpl){hitplane_=hpl;}
  inline void SetHitPlaneSide(unsigned int hpls){hitplaneside_=hpls;}
  inline void SetHitDetRawId(uint32_t  drid){detrawid_=drid;}

  inline void SetHitR(float hitR) {hitR_=hitR;}
  inline void SetHitPhi(float hitPhi) {hitPhi_=hitPhi;}
  inline void SetHitZ(float hitZ) {hitZ_=hitZ;}
  inline void SetHitDR(float hitDR) {hitDR_=hitDR;}
  inline void SetHitDPhi(float hitDPhi) {hitDPhi_=hitDPhi;}
  inline void SetHitDZ(float hitDZ) {hitDZ_=hitDZ;}
    
  inline void SetHitX(float hitX) {hitX_=hitX;}
  inline void SetHitDX(float hitDX) {hitDX_=hitDX;}
  inline void SetHitY(float hitY) {hitY_=hitY;}
  inline void SetHitDY(float hitDY) {hitDY_=hitDY;}

  inline void SetHitClass(unsigned int hitClass) {hitClass_=hitClass;}                      
  inline void SetHitNumStrip(unsigned int hit_numstrip) {hit_numstrip_=hit_numstrip;}             // 1: Strip Cluster + Pad Cluster 2: Other configuration
  inline void SetHitNumPad(unsigned int hit_numpad) {hit_numpad_=hit_numpad;}
  inline void SetHitNumPadCol(unsigned int hit_numpadcol) {hit_numpadcol_=hit_numpadcol;}
  inline void AddFormingCluster(hit_entry ahit_entry) {thehitentries.push_back(ahit_entry);}      // add a cluster in the group of the clusters belong to 1 hit. 
  /*
   bool IsVtxHit() const {return VtxHits;} 
   bool IsInTrk() const {return HitInTrk;} 
   void SetVtxHit(bool theVtxHits) {VtxHits=theVtxHits;}
   void SetInTrk(bool theHitInTrk) {HitInTrk=theHitInTrk;}
  */

  void ComputeHit();

  //typedef vector<cluster_entry> cluster_entries_vector; 

  //Added for noise studies
  inline const std::vector<cluster_entry> & GetCluStripEntries() const {return ClusterStrip_entries;}
  inline const std::vector<cluster_entry> & GetCluPadEntries() const {return ClusterPad_entries;}
  inline void AddCluPadEntries(std::vector<cluster_entry> tocpp)
    {
      
      ClusterPad_entries.clear();
      for(unsigned i=0; i<tocpp.size();i++)
	{
	  ClusterPad_entries.push_back(tocpp.at(i));
	}
    };

  inline void AddCluStripEntries(std::vector<cluster_entry> tocps)
    {
      ClusterStrip_entries.clear();
      for(unsigned i=0; i<tocps.size();i++)
	{
	  ClusterStrip_entries.push_back(tocps.at(i));
	}
    };   


  friend bool operator<(const T2Hit &h1,const T2Hit &h2);
  friend bool operator>(const T2Hit &h1,const T2Hit &h2);
  friend bool operator==(const T2Hit &h1,const T2Hit &h2);
  friend bool operator!=(const T2Hit &h1,const T2Hit &h2);

  std::vector<cluster_entry>  ClusterStrip_entries; 
  std::vector<cluster_entry>  ClusterPad_entries;

//  bool UtilizedForTracking;
  //std::pair<int,std::vector<int> > PadCluIdToMatchingStripIds;

  std::pair<long int,long int> HitUniqueId; //PadCLUId,StripCLUId, -1 if some cluster is missing.
  std::vector<long int>  StripVectorBelowThePadId; //useful in track studies in order to know what strip cluster are associatet to a pad.

  unsigned int Hit_PosInCollection; //Associate at the Hit  the position in the T2HitCollection


  // const int hitRoadID;

 private:

  float hitR_;
  float hitPhi_;
  float hitZ_;
  float  hitDR_;
  float  hitDPhi_;
  float  hitDZ_;
  
  float hitX_;
  float hitY_;
  float  hitDX_;
  float  hitDY_;
 
  bool VtxHits;
  bool HitInTrk;


  unsigned int hithalftele_;
  unsigned int hitplane_;
  unsigned int hitplaneside_;
  unsigned int hitarm_;
  uint32_t detrawid_;

  unsigned int hitClass_;
   
  unsigned int hit_numstrip_;              // number of strips forming the hit 
  unsigned int  hit_numpad_; 
  unsigned int  hit_numpadcol_;  
                
  vecthit_entries thehitentries;           // vector of clusters data forming hit

};




#endif
