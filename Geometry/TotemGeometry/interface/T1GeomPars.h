#ifndef T1_GEOM_PARS
#define T1_GEOM_PARS

#include <string>

using namespace std;

class T1GeomPars {

 public:
  T1GeomPars(string parFileName = "");
  ~T1GeomPars() {}

 private:
  void _defaultValues();
  int _readFromFile(string parFileName);

  // Geometry parameters
  float _pitchW;
  float _pitchS;
  int _nWires[5][2];
  int _nStrips[5][2];
  float _wireOffset[5][6];  // in n. of wires
  float _rW0[5];  // radial position of wire 0
  float _drS0[5][6];  // radial distance between strip-0 intersection and wire 0
  float _dx[6];  // tangent offset of sensitive volume axis wrt. radius (seen from IP)
  float _zLayer[5];  // z position of layers
  float _dz[6];  // CSC z offset wrt. layer position
  float _tilt[5];  // layer tilt angle

 public:

  // Getters for geometry parameters
  const float getPitchW() {return _pitchW;}
  const float getPitchS() {return _pitchS;}
  const float getNWires(int layer, int type) {return _nWires[layer][type];}
  const float getNStrips(int layer, int type) {return _nStrips[layer][type];}
  const float getWireOffset(int layer, int sextant) {return _wireOffset[layer][sextant];}
  const float getRW0(int layer) {return _rW0[layer];}
  const float getDrS0(int layer, int sextant) {return _drS0[layer][sextant];}
  const float getDx(int sextant) {return _dx[sextant];}
  const float getZLayer(int layer) {return _zLayer[layer];}
  const float getDz(int sextant) {return _dz[sextant];}
  const float getTilt(int layer) {return _tilt[layer];}

  const float xzParity(int arm) {return (arm == 0) ? 1. : -1.;}
};

#endif
