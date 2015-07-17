#include "MyBeamSpot.h"
#include <iostream>

using namespace std;

ClassImp(MyBeamSpot)

MyBeamSpot::MyBeamSpot(){
  this->Reset();
}

MyBeamSpot::~MyBeamSpot(){}

void MyBeamSpot::Print() {

  cout<<"beam spot information: "<<endl;
  cout<<"x: "<<this->x<<endl;
  cout<<"y: "<<this->y<<endl;
  cout<<"z: "<<this->z<<endl;
 
  cout<<"error x: "<<this->ex<<endl;
  cout<<"error y: "<<this->ey<<endl;
  cout<<"error z: "<<this->ez<<endl;

  cout<<"sigma z: "<<this->sigmaZ<<endl;
  cout<<"dxdz: "<<this->dxdz<<endl;
  cout<<"dydz: "<<this->dydz<<endl;

  cout<<"error sigma z: "<<this->esigmaZ<<endl;
  cout<<"error dxdz: "<<this->edxdz<<endl;
  cout<<"error dydz: "<<this->edydz<<endl;

  cout<<"BeamWidthX: "<<this->BeamWidthX<<endl;
  cout<<"BeamWidthY: "<<this->BeamWidthY<<endl;

  cout<<"error BeamWidthX: "<<this->eBeamWidthX<<endl;
  cout<<"error BeamWidthY: "<<this->eBeamWidthY<<endl;
}
