/*
Hit class for analysis in TotemTest
*/

#ifndef TOTEMHIT_H
#define TOTEMHIT_H

#include<iostream>

class TotemHit{

 public:

//  TotemHit(int det,float eloss, float pabs, int pid, int trackid,
//               int parent,float x, float y, float z,float vx, float vy, float vz);

  TotemHit(int, float, float, int, int,
               int,float, float, float,float, float, float,float, float, float,float, float, float);
TotemHit(){}
virtual ~TotemHit();
int det() {return _det;}
 int pid() {return _pid;}
 int trackid() {return _trackid;}
 int parent() {return _parent;}
  float eloss() {return _eloss;}
 float pabs() {return _pabs;}
 float x() {return _x;}
float y() {return _y;}
float z() {return _z;}
float vx() {return _vx;}
float vy() {return _vy;}
float vz() {return _vz;}
 float px() {return _px;}
 float py() {return _py;}
 float pz() {return _pz;}

 float vpx() {return _vpx;}
 float vpy() {return _vpy;}
 float vpz() {return _vpz;}


 void Print();

 private:

  int _det, _pid, _trackid, _parent;
  float _eloss, _pabs, _x, _y, _z, _vx, _vy, _vz, _px,_py,_pz, _vpx,_vpy,_vpz;



};






#endif
