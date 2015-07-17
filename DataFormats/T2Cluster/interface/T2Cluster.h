/****************************************************************************
 *
 * This is a part of TOTEM testbeam/monitoring software.
 * Authors: 
 *	 Hubert Niewiadomski
 *    
 * $Id: T2Cluster.h,v 1.1 2008/09/12 16:19:56 lgrzanka Exp $
 * $Revision: 1.1 $
 * $Date: 2008/09/12 16:19:56 $
 *
 ****************************************************************************/
//Hubert Niewiadomski: the GEM cluster element
//Giuseppe Latino: geometrical reconstruction of GEM clusters added (May 07).
//Mirko Berretti: ComputeClusterParamsForTracking (Dec 2010)

#ifndef T2DETCLUSTERS__
#define T2DETCLUSTERS__

#include <stdint.h>
#include <vector>
#include <iostream>
#include "Rtypes.h"
#include "cluster_entry.h"
#include <map>

using namespace std;


class T2Cluster
{
 public:
  virtual ~T2Cluster() {}
  typedef vector<cluster_entry> cluster_entries_vector;
  enum cluster_type {pad, strip};
    
  inline unsigned short GetNoOfEntries() const {return entry_numb_;}
  inline short GetRadialSpread() const {return radial_spread_;}
  inline short GetAngularSpread() const {return angular_spread_;}
  inline float GetRadialCentrePos() const {return radial_centre_;}
  inline float GetAngularCentrePos() const {return angular_centre_;}
  inline float GetClusterR() const {return R_centre_;}
  inline float GetClusterDR() const {return DR_centre_;}
  inline float GetClusterPhi() const {return Phi_centre_;}
  inline float GetClusterDPhi() const {return DPhi_centre_;}

  inline cluster_type GetClusterType() const {return type_;}
  // inline unsigned int GetDetID() const {return det_id_;}
  inline uint32_t GetDetID() const {return det_id_;}  //Bug found int 12/1/2011
  inline const cluster_entries_vector & GetEntries() const {return cluster_entries;}
  //    
  inline void SetNoOfEntries(unsigned short entry_numb) {entry_numb_ = entry_numb;}
  inline void SetRadialSpread(unsigned short radial_spread) {radial_spread_ = radial_spread;}
  inline void SetAngularSpread(unsigned short angular_spread) {angular_spread_ = angular_spread;}
  inline void SetRadialCentrePos(float radial_centre_pos) {radial_centre_ = radial_centre_pos;}
  inline void SetAngularCentrePos(float angular_centre_pos) {angular_centre_ = angular_centre_pos;}
  inline void SetClusterR(float R_centre_pos) {R_centre_ = R_centre_pos;}
  inline void SetClusterDR(float DR_centre_pos) {DR_centre_ = DR_centre_pos;}
  inline void SetClusterPhi(float Phi_centre_pos) {Phi_centre_ = Phi_centre_pos;}
  inline void SetClusterDPhi(float DPhi_centre_pos) {DPhi_centre_ = DPhi_centre_pos;}
  inline void SetClusterType(cluster_type clu_type) {type_=clu_type;}
  inline void SetDetID(uint32_t det_id) {det_id_ = det_id;}

  inline void AddEntry(const cluster_entry &point) {cluster_entries.push_back(point);}
  void ComputeClusterParams();

  void ComputeClusterParamsForTracking(double Projectionthreshold,int BlobMinSize);

  void Blobangular_centre_spread(std::vector<double> &AllColscenters,float &angular_centre_,unsigned short &angular_spread_);
  void Blobradial_centre_spread(std::vector<double> &AllRowscenters,float &radial_centre_,unsigned short &radial_spread_);
  void BlobPhi_Dphi(std::vector<double> &AllAzimuthcenters,std::vector<double> &AllAzimDphi,float &Phi_centre_,float &DPhi_centre); 
  void BlobR_DR(std::vector<double> &AllRadialcenters,std::vector<double> &AllRadialDR,float &R_centre_,float &DR_centre_);
  
  //Used in T2RoadPadProducer
  //std::pair<unsigned int, unsigned int> PadStatus_RoadId;
  long int ClustId_unique;

 private:
  uint32_t det_id_;  //id of the detector
  cluster_type type_;   //type of the cluster
  unsigned short entry_numb_;  //number of pads/strips in the cluster
  unsigned short radial_spread_;  //number of rows covered
  unsigned short angular_spread_;  //number of columns covered
  float radial_centre_;  //centre of gravity in radial coordinate, in rows
  float angular_centre_;  //centre of gravity in phi clockwise coordinate, in columns
  float R_centre_; //centre of gravity in radial coordinate, in mm
  float DR_centre_; //spread in radial coordinate, in mm //NOT EXACT!!!
  float Phi_centre_; //centre of gravity in phi anti-clockwise coordinate, in degrees
  float DPhi_centre_; //spread in phi coordinate, in degrees
  cluster_entries_vector cluster_entries;  //vector of pads/strips forming a cluster
  //ClassDef(T2Cluster,1)
};

ostream & operator<<(ostream & out, const T2Cluster& clust);


#endif
