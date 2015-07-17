#ifndef Reco_FUNC22h
#define Reco_FUNC22h


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

class MyFCN22: public ROOT::Minuit2::FCNBase { 

public: 

  MyFCN22(vector<TVector3> xyz, vector<TVector2> xy_prev, vector<TVector2> xy_error) {


    _xyz = xyz;
    _xy_prev = xy_prev;
    _xy_error = xy_error;
}

  MyFCN22() {}

  void setMyFCN22(vector<TVector3> xyz, vector<TVector2> xy_prev, vector<TVector2> xy_error) {

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
	if(fabs(_xyz[i].Z())<8000){
	  XX = _xyz[i].X()*cos(par[10])-_xyz[i].Y()*sin(par[10])+par[0];
	  YY = _xyz[i].X()*sin(par[10])+_xyz[i].Y()*cos(par[10])+par[1];
	}       
	if(fabs(_xyz[i].Z())<8400 && fabs(_xyz[i].Z())>8000){
	  XX = _xyz[i].X()*cos(par[11])-_xyz[i].Y()*sin(par[11])+par[2];
	  YY = _xyz[i].X()*sin(par[11])+_xyz[i].Y()*cos(par[11])+par[3];
	}
	if(fabs(_xyz[i].Z())<9000 && fabs(_xyz[i].Z())>8400){
	  XX = _xyz[i].X()*cos(par[12])-_xyz[i].Y()*sin(par[12])+par[4];
	  YY = _xyz[i].X()*sin(par[12])+_xyz[i].Y()*cos(par[12])+par[5];

 	}
	if(fabs(_xyz[i].Z())<9500 && fabs(_xyz[i].Z())>9000){

 	  XX = _xyz[i].X()*cos(par[13])-_xyz[i].Y()*sin(par[13])+par[6];
	  YY = _xyz[i].X()*sin(par[13])+_xyz[i].Y()*cos(par[13])+par[7];


 	}
	if(fabs(_xyz[i].Z())>9500){
 	  XX = _xyz[i].X()*cos(par[14])-_xyz[i].Y()*sin(par[14])+par[8];
	  YY = _xyz[i].X()*sin(par[14])+_xyz[i].Y()*cos(par[14])+par[9];

 	}


/*
      double XX = _xyz[i].X()*cos(par[2])-_xyz[i].Y()*sin(par[2])+par[0];
      double YY = _xyz[i].X()*sin(par[2])+_xyz[i].Y()*cos(par[2])+par[1];
*/
	CHI += (XX-_xy_prev[i].X())*(XX-_xy_prev[i].X())/_xy_error[i].X()/_xy_error[i].X()+(YY-_xy_prev[i].Y())*(YY-_xy_prev[i].Y())/_xy_error[i].Y()/_xy_error[i].Y();
      




    }

/*
    cout << _xy_prev[0].X() << " " << _xy_prev[0].Y() << endl;
  CHI = (par[0]-_xy_prev[0].X())+(par[1]-_xy_prev[0].Y());

*/
    return CHI;

    // Rosebrock function
//    return  fA*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]) + fB*(1 - x[0])*(1 - x[0]);

  } 
  
  double Up() const { return 1.; }

private: 

  vector<TVector3> _xyz;
  vector<TVector2> _xy_prev;
  vector<TVector2> _xy_error;
};
#endif
