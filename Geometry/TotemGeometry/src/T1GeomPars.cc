#include "Geometry/TotemGeometry/interface/T1GeomPars.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace edm;

T1GeomPars::T1GeomPars(string parFileName) {
  _defaultValues();
  if (parFileName.compare("") != 0) {
    int status = _readFromFile(parFileName);
    if (status != 0)
      LogInfo("T1GeomPars") << "Could not read parameters from file " << parFileName
			    << ". Using default values.";
  }
}

void T1GeomPars::_defaultValues() {
  _pitchW = 3.;
  _pitchS = 5.;

  _nWires[0][0] = 165;
  _nWires[0][1] = 127;
  _nWires[1][0] = 181;
  _nWires[1][1] = 150;
  _nWires[2][0] = 196;
  _nWires[2][1] = 158;
  _nWires[3][0] = 213;
  _nWires[3][1] = 180;
  _nWires[4][0] = 226;
  _nWires[4][1] = 204;

  _nStrips[0][0] = 118;
  _nStrips[0][1] = 97;
  _nStrips[1][0] = 129;
  _nStrips[1][1] = 108;
  _nStrips[2][0] = 138;
  _nStrips[2][1] = 114;
  _nStrips[3][0] = 151;
  _nStrips[3][1] = 127;
  _nStrips[4][0] = 164;
  _nStrips[4][1] = 150;

  for (int iSext = 0; iSext < 6; iSext++) {
    _wireOffset[0][iSext] = -54.66;
    _wireOffset[1][iSext] = -57.80;
    _wireOffset[2][iSext] = -60.82;
    _wireOffset[3][iSext] = -64.50;
    if (iSext == 2 || iSext == 5)
      _wireOffset[4][iSext] = -71.26;
    else
      _wireOffset[4][iSext] = -62.50;
  }

  _rW0[0] = 144.5;
  _rW0[1] = 155.5;
  _rW0[2] = 164.5; 
  _rW0[3] = 176.5;
  _rW0[4] = 185.5;

  for (int iSext = 0; iSext < 6; iSext++) {
    _drS0[0][iSext] = 164.46;
    _drS0[1][iSext] = 174.42;
    _drS0[2][iSext] = 183.42;
    _drS0[3][iSext] = 194.42;
    if (iSext == 2 || iSext == 5)
      _drS0[4][iSext] = 213.91;
    else
      _drS0[4][iSext] = 187.93;
  }

  _dx[0] = 3.035;
  _dx[1] = -3.035;
  _dx[2] = -3.035;
  _dx[3] = 3.035;
  _dx[4] = -3.035;
  _dx[5] = 3.035;

  _zLayer[0] = 7543.;
  _zLayer[1] = 8213.;
  _zLayer[2] = 8820.;
  _zLayer[3] = 9427.;
  _zLayer[4] = 10190.;

  for (int iSext = 0; iSext < 6; iSext++)
    _dz[iSext] = (iSext%2 == 0) ? 25. : -25.;

  _tilt[0] = -6.;
  _tilt[1] = -3.;
  _tilt[2] = 6.;
  _tilt[3] = 3.;
  _tilt[4] = 0.;
}

int T1GeomPars::_readFromFile(string parFileName)
{
  return 0;
}
