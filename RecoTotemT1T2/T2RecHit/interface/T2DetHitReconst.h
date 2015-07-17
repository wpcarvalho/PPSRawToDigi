/****************************************************************************
*
* 
* Authors: Mirko Berretti   
* mirko.berretti@gmail.com   
* University of Siena
* $Revision: 1.2 $
* $Date: 2009/03/20 10:49:24 $
*
****************************************************************************/


#ifndef T2DETHITRECONS__
#define T2DETHITRECONS__

#include <vector>
#include <map>
#include <utility>
 

#include "Rtypes.h"
#include "DataFormats/T2Cluster/interface/T2Cluster.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"


using namespace std;

class T2DetHitReconst                            //Find the matching cluster in the plane to give as input for T2Hit class; After Reconstruct the plane-hits 
{
  public:
  //T2DetHitReconst();
  T2DetHitReconst(uint32_t det_id,unsigned int  Cl1MaxPad_, unsigned int  Cl1MaxStrip_);     //Raw Id of the detector is needed to calculate z of the Hit
  virtual ~T2DetHitReconst();

  //bool MatchClust(T2Cluster c1, T2Cluster c2);
    //   void FindMatchingClusters();                 //Loading of matchingplanecl (type VMathcingCluster) 
    void FindHits();

    typedef vector<T2Cluster> VCluster;
    //typedef vector<VCluster> VMatchingCluster;   //1 Plane-> more possible hits; 1 Hits-> more possible forming-clusters
    //typedef pair<T2Cluster,bool> ClustBool;      //necessary type to avoid double counting in the FindMatchingCluster() procedure
    //typedef vector<ClustBool> ClustBoolV; 
    
    void AddStripClusters(VCluster StripClv)     //Take all the strip ClusterS of a plane
      {strip_clusters=StripClv;}
    
    void AddPadClusters(VCluster PadClv)        
      {pad_clusters = PadClv;}
   
    unsigned int GetDetId() 
      {return det_id_;}
 
    float GetDetZ()
      {return Zdet;}
    
    
    void CalculateZDet();
   
    
   
        
    // VMatchingCluster TakeMatchedClusters()
    //{return matchingPlaneCl;}
    
    T2HitCollection TakePlaneHits()            // Return all the hits of a plane
      {return theHitV;}
    

    




    inline void ClearClusterContainers() 
      {pad_clusters.clear(); strip_clusters.clear(); /*matchingPlaneCl.clear();*/}

   
  protected:

    VCluster pad_clusters;
    VCluster strip_clusters;
    // VMatchingCluster matchingPlaneCl;
    
    uint32_t det_id_;
    unsigned int  Cl1MaxPad, Cl1MaxStrip;

    T2HitCollection theHitV;                    // The final output: a vector of hits

   
    
    T2DetId currentDet; 
    float Zdet;
    float DZdet; 
};

#endif
