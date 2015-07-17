/**
 * Class 
 *
 * Author: Mirko Berretti / University of Siena 
 * Email:  mirko.berretti@gmail.com
 * Date:   2010-12-05
 */

#ifndef _Tube_h
#define _Tube_h
#include <vector>
#include <map>
#include <utility>
#include "DataFormats/T2Cluster/interface/T2Cluster.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Cluster/interface/cluster_entry.h"
//#include "DataFormats/T2Road/interface/T2Road.h"


class Tube {
public:
  /// Constructor
  Tube();


  /// Destructor
  virtual ~Tube();

  unsigned int uniquetubeId;
  int seedId;
  std::vector<T2Cluster> ClustV;
  std::vector<int> uniqueCluIdV;

  std::vector<double> VectR;
  std::vector<double> VectPhi;
  std::vector<double> VectX;
  std::vector<double> VectY;
  std::vector<double> VectZ;
  std::vector<double> VectEX;
  std::vector<double> VectEY;
  
  double Ax,Ay,Xi,Yi,Xf,Yf,Zi,Zf,Ri,Rf,Phii,Phif,planei,planef;
  
  
  // Update Ax,Ay,X0,Y0,Xi,Yi,Xf,Yf;
  void updateTubeparam(T2Cluster addNewCluster);
  
  void ResetTube();

  // void Tube_extrapolation(double wantedZ, double &retX, double &retY, double &retEX, double &retEY);

};
#endif
