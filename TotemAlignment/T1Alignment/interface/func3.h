#ifndef Reco_FUNC2_h
#define Reco_FUNC2_h


#include <vector>
#include "TF1.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TStyle.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TVector2.h"


#define RADICE3 1.732050808
#define PITCH 5.0
#define _30DEG 0.5236
#define _25DEG 0.4363
#define _60DEG 1.0472
#define _90DEG 1.5708
#define _120DEG 2.0944
#define _150DEG 2.618
#define _210DEG 3.665
#define _240DEG 4.1888
#define _270DEG 4.7124
#define _300DEG 5.2360
#define _330DEG 5.7596
#define _3DEG 0.05239878


class MyFCN3 : public ROOT::Minuit2::FCNBase { 

public: 

  MyFCN3(vector<TVector3> xyz, vector<TVector2> xy_prev, vector<TVector2> xy_error) {


    _xyz = xyz;
    _xy_prev = xy_prev;
    _xy_error = xy_error;
}

  MyFCN3() {}

  void setMyFCN3(vector<TVector3> xyz, vector<TVector2> xy_prev, vector<TVector2> xy_error) {

    _xyz = xyz;
    _xy_prev = xy_prev;
    _xy_error = xy_error;

    return;
  }

  double operator() (const std::vector<double> & par) const {
  
    assert(_xyz.size() == _xy_prev.size());
    assert(_xyz.size() != 0);
    assert(_xy_prev.size() != 0);
    double CHI=0;
    double XX=0;
    double YY=0;
    for(unsigned int i=0; i< _xyz.size();i++){

      int sestante = sextant((double)_xyz[i].X(),(double)_xyz[i].Y(),(double)_xyz[i].Z());

      switch(sestante){
      case 0:
	
	  XX = _xyz[i].X()+par[0];
	  YY = _xyz[i].Y()+par[1];
	  break;

      case 1:
	
	  XX = _xyz[i].X()+par[2];
	  YY = _xyz[i].Y()+par[3];
	  break;

      case 2:
	
	  XX = _xyz[i].X()+par[4];
	  YY = _xyz[i].Y()+par[5];
	  break;

      case 3:
	
	  XX = _xyz[i].X()+par[6];
	  YY = _xyz[i].Y()+par[7];
	  break;

      case 4:
	
	  XX = _xyz[i].X()+par[8];
	  YY = _xyz[i].Y()+par[9];
	  break;

      case 5:
	
	  XX = _xyz[i].X()+par[10];
	  YY = _xyz[i].Y()+par[11];
	  break;



      }


	CHI += (XX-_xy_prev[i].X())*(XX-_xy_prev[i].X())/_xy_error[i].X()/_xy_error[i].X()+(YY-_xy_prev[i].Y())*(YY-_xy_prev[i].Y())/_xy_error[i].Y()/_xy_error[i].Y();
      




    }


    return CHI;

  

  } 
  
  double Up() const { return 1.; }


double Eta(double x,double y,double z){
  double xyt;
  double c=0;
  double eta2;
  xyt = sqrt(x*x + y*y);
  //theta
  if(z>0) c = atan(xyt/z);
  if(z<0) c = atan(xyt/z)+3.14159;
  if(z==0) {c = 3.14159;}
  //pseudorapidity
  eta2 = -log(tan(c/2.));
  return eta2;
}
double Phi(double x,double y) const{
  double c=0;
  if(x>0 && y>0) c = atan(y/x);
  if(x<0) c = atan(y/x)+3.14159;
  if(x>0 && y<0) c = atan(y/x)+6.28318;
  return c;
}

int sextant(double x, double y, double z) const{
    int sestante = -1;
    double fi = Phi(x,y);
//  cout << " FI " << fi << endl;
    if(z>0){
      if(plane(z)==4){
	if( fi < _30DEG || fi > (2*PI - _30DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG) && fi < (_60DEG +_30DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG) && fi < (_120DEG +_30DEG) ) sestante = 1;
	if(fi > (PI -_30DEG) && fi < (PI +_30DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG) && fi < (_240DEG +_30DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG) && fi < (_300DEG +_30DEG) ) sestante = 4;
      } 
      if(plane(z)==3){
	if( fi < (_30DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG-_3DEG) && fi < (PI +_30DEG-_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG) ) sestante = 4;
      } 
      if(plane(z)==2){
	if( fi < (_30DEG-_3DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG-_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG-_3DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG-_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG-_3DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG-_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG-_3DEG-_3DEG) && fi < (PI +_30DEG-_3DEG-_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG-_3DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG-_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG-_3DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG-_3DEG) ) sestante = 4;
      } 
      if(plane(z)==1){
	if( fi < (_30DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG+_3DEG) && fi < (PI +_30DEG+_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG) ) sestante = 4;
      } 
      if(plane(z)==0){
	if( fi < (_30DEG+_3DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG+_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG+_3DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG+_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG+_3DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG+_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG+_3DEG+_3DEG) && fi < (PI +_30DEG+_3DEG+_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG+_3DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG+_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG+_3DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG+_3DEG) ) sestante = 4;
      } 
    
    }


    if(z<0){
      if(plane(z)==4){
	if( fi < _30DEG || fi > (2*PI - _30DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG) && fi < (_60DEG +_30DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG) && fi < (_120DEG +_30DEG) ) sestante = 0;
	if(fi > (PI -_30DEG) && fi < (PI +_30DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG) && fi < (_240DEG +_30DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG) && fi < (_300DEG +_30DEG) ) sestante = 3;
      } 
      if(plane(z)==3){
	if( fi < (_30DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG+_3DEG) && fi < (PI +_30DEG+_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG) ) sestante = 3;
      } 
      if(plane(z)==2){
	if( fi < (_30DEG+_3DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG+_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG+_3DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG+_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG+_3DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG+_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG+_3DEG+_3DEG) && fi < (PI +_30DEG+_3DEG+_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG+_3DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG+_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG+_3DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG+_3DEG) ) sestante = 3;
      } 
      if(plane(z)==1){
	if( fi < (_30DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG-_3DEG) && fi < (PI +_30DEG-_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG) ) sestante = 3;
      } 
      if(plane(z)==0){
	if( fi < (_30DEG-_3DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG-_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG-_3DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG-_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG-_3DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG-_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG-_3DEG-_3DEG) && fi < (PI +_30DEG-_3DEG-_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG-_3DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG-_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG-_3DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG-_3DEG) ) sestante = 3;
      } 
    
    }

  
    return sestante;
  }
  int plane(double z) const {

    int piano = -1;
    if(fabs(z)<8000){
      piano = 0;
 
    }       
    if(fabs(z)<8400 && fabs(z)>8000){
      piano = 1;
    }
    if(fabs(z)<9000 && fabs(z)>8400){
      piano = 2;
    }
    if(fabs(z)<9500 && fabs(z)>9000){
      piano = 3;
    }
    if(fabs(z)>9500){
      piano = 4;
    }
    return piano;
 	
  }





private: 

  vector<TVector3> _xyz;
  vector<TVector2> _xy_prev;
  vector<TVector2> _xy_error;
};
#endif
