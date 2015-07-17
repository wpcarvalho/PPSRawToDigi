/****************************************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors: 
*	 Giuseppe Latino
*    
* $Id: T2DetClustReconst.h,v 1.3.2.1 2009/11/07 20:05:51 berretti Exp $
* $Revision: 1.3.2.1 $
* $Date: 2009/11/07 20:05:51 $
*
****************************************************************************/
//Created by Giuseppe Latino: reconstruction of clusters (September 06).
//Hubert Niewiadomski: cluster storage added (November 06).
//Giuseppe Latino: geometrical reconstruction of clusters (May 07).


#ifndef T2DETCLUSTRECONST__
#define T2DETCLUSTRECONST__

#include <vector>
#include <map>
#include "DataFormats/T2Cluster/interface/T2Cluster.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"

using namespace std;

class T2DetClustReconst
{

public:

  T2DetClustReconst ();
  virtual ~ T2DetClustReconst ();
  void SetStripHits (vector<int> strip_ch);
  void SetPadHits (vector<int> pad_ch);
  void BuildClusters ();
  void SetDetId (uint32_t det_id)
  {
    det_id_ = det_id;
  }
  
  uint32_t  GetDetId ()
  {
    return det_id_;
  }

  typedef vector<T2Cluster> clusters_vector;

  inline const clusters_vector & GetStripClusters () const
  {
    return strip_clusters;
  }

  inline const clusters_vector & GetPadClusters () const
  {
    return pad_clusters;
  }

  typedef map<int, map<int,bool> > HitsMap;
  typedef map<int, map<int,int> > ChMap;

  

  HitsMap StrHits, PadHits;
  ChMap StrCluStrRawID, StrCluStrColID;
  ChMap PadCluPadRawID, PadCluPadColID;

  inline const ChMap GetMapDet_PadColRowID()
    {
      return PadCluPadColID;
    }


  inline const ChMap GetMapDet_StrColRowID()
    {
      return StrCluStrColID;
    }


  double Projectionthreshold;
  int BlobMinSize;
 
protected:

  clusters_vector pad_clusters;
  clusters_vector strip_clusters;
  uint32_t det_id_;

  inline void ClearClusterContainers ()
  {
    pad_clusters.clear();
    strip_clusters.clear();
  }

};

#endif
