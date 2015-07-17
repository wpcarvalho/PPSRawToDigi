/*
  Class T2Road
  Author: Mirko Berretti
  mirko.berretti@gmail.com
*/
#ifndef DataFormats_T2Road_H
#define DataFormats_T2Road_H

#include <vector>
#include <iterator>

#include <DataFormats/T2Hit/interface/T2Hit.h>
#include <DataFormats/T2Cluster/interface/T2Cluster.h>

class T2Road
{
 public:
  T2Road() {}
  std::vector<T2Hit>  thisRoad;
  std::vector<T2Cluster>  thisPadRoad;
  
  void SetRoadDR(float dR) {roadDR_=dR;}
  void SetRoadDPhi(float dPhi){roadDPhi_=dPhi;}
  double GetRoadDPhi() const {return roadDPhi_;}
  double GetRoadDR() const {return roadDR_;}
  unsigned int GetRoadSize() const {return thisRoad.size();}
    
  double GetRoadPhimin() const {return roadPhimin;}
  double GetRoadPhimax() const {return roadPhimax;}
  double GetRoadRmin() const {return roadRmin;}
  double GetRoadRmax() const {return roadRmax;}
    
  void SetRoadPhimin(double rm) {roadPhimin=rm;}
  void SetRoadPhimax(double rM) {roadPhimax=rM;}
  void SetRoadRmin(double pm) {roadRmin=pm;}
  void SetRoadRmax(double pM)  {roadRmax=pM;}
  
  void CalculateRoadExtreme();
 

 unsigned int GetRoadID() const {return RoadID;}
 void SetRoadID(unsigned int roadid_){RoadID=roadid_;}
 
unsigned int RoadID;

 private:
  double roadDR_;
  double roadDPhi_;
   
  double roadPhimin;
  double roadPhimax;
  double roadRmin;
  double roadRmax;
  
  //unsigned int RoadType; //0=Old-standard; 1=PadRoad.

};




#endif
