#ifndef T1Cluster_T1Cluster_h
#define T1Cluster_T1Cluster_h

#include <boost/cstdint.hpp>
using namespace std;
class T1Cluster{

 public:
  explicit T1Cluster (int first, int last, float center, float sigma, float width, int bx );
  T1Cluster ();
 
  ~T1Cluster(){};
 
  bool operator==(const T1Cluster& cl) const;
  bool operator<(const T1Cluster& cl) const;

  int firstStrip() const;
  int lastStrip() const;
  float Center() const;
  float Sigma() const;
  float Width() const;
  int bx() const;
  
  //  void print() const ;
 private:
  int _fstrip;
  int _lstrip;
  int _bunchx;
  float _center;
  float _sigma;
  float _width;


};
#endif
