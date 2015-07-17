#include "RecoTotemT1T2/T2RoadPadFinder/interface/Tube.h"

Tube::Tube()
{}

 


Tube::~Tube()
{}

void Tube::ResetTube(){

  //uniquetubeId;
  //seedId;

  Ax=0;Ay=0;Xi=0;Yi=0;Xf=0;Yf=0;Zi=0;Zf=0;Ri=0;Rf=0;Phii=0;Phif=0;
  planei=0;
  planef = 0.0;

  ClustV.clear();
  uniqueCluIdV.clear();

  VectR.clear();
  VectPhi.clear();
  VectX.clear();
  VectY.clear();
  VectZ.clear();
  VectEX.clear();
  VectEY.clear();
   
}

void Tube::updateTubeparam(T2Cluster addNewCluster)
{
  std::cout<<"Update-Tube"<<std::endl;
  ClustV.push_back(addNewCluster);
}
  
/*
void Tube::Tube_extrapolation(double wantedZ, double &retX, double &retY, double &retEX, double &retEY)
{
  retX=0.;
  retY=0.; 
  retEX=0.; 
  retEY=0.;
  retX=retX*wantedZ;
}
*/
