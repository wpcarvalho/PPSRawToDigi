#include <DataFormats/T1Cluster/interface/T1Cluster.h>
using namespace std;

T1Cluster::T1Cluster (int first, int last, float center, float sigma,float width, int bx ) :
  _fstrip(first),
  _lstrip(last),
  _bunchx(bx), 
  _center(center),
  _sigma(sigma),
  _width(width)
{}

T1Cluster::T1Cluster ():
  _fstrip(0),
  _lstrip(0),
  _bunchx(0), 
  _center(0),
  _sigma(0),
  _width(0)
{}


int T1Cluster::firstStrip() const {

  return _fstrip;
}
int T1Cluster::lastStrip() const {

  return _lstrip;
}
float T1Cluster::Center() const {

  return _center;
}
float T1Cluster::Sigma() const {

  return _sigma;
}
float T1Cluster::Width() const {

  return _width;
}
int T1Cluster::bx() const {

  return _bunchx;
}

bool T1Cluster::operator<(const T1Cluster& cl) const{
  
  if(cl.bx() == this->bx())
    return cl.Center()<this->Center();
  else 
    return cl.bx()<this->bx();
}

bool 
T1Cluster::operator==(const T1Cluster& cl) const {
  return ( (this->Center() == cl.Center()) &&
	   (this->bx()          == cl.bx())          && 
	   (this->firstStrip()  == cl.firstStrip()) );
}
/*
  void
  T1Cluster::print() const {
  std::cout << "First Strip " << firstStrip() << std::endl;
  std::cout << "Last Strip " << lastStrip() << std::endl;
  std::cout << "Center Strip " << Center()<< std::endl;
  std::cout << "Sigma " << Sigma()<< std::endl;
  std::cout << " bx " << bx() <<std::endl;
  }
*/
