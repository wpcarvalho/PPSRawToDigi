#ifndef _CFECS_
#define _CFECS_

class cfecConn
{
public:
 
cfecConn(int cfecSide, unsigned int view, unsigned int firstStrip, unsigned int lastStrip) :
  _cfecSide(cfecSide), _view(view), _firstStrip(firstStrip), _lastStrip(lastStrip)
{
  _dir = _cfecSide*2*(_view - 1.5);
  for (int iHalf = 0; iHalf < 2; iHalf++) {
    int nStrips = std::min(16, _lastStrip - (_firstStrip + 16*iHalf) + 1);
    if (nStrips <= 0)
      _firstPin[iHalf] = 16;
    else if (nStrips <= 6)
      _firstPin[iHalf] = (_dir == 1) ? 10 : 6 - nStrips;
    else
      _firstPin[iHalf] = (_dir == 1) ? 0 : 16 - nStrips;
 }
}

  ~cfecConn() {};
  bool isConnected() {if (_cfecSide != 0 && _view != 0) return true; else return false;};
  unsigned int getView() {return _view;};
  int getStrip(unsigned int channel){
  int half = channel/16;
  int startStrip = (_dir == 1) ? _firstStrip + 16*half : std::min(_lastStrip, _firstStrip + 15 + 16*half);
  int strip = startStrip + _dir*(channel - (_firstPin[half] + 16*half));
  if (strip < _firstStrip + 16*half || strip > _lastStrip)
    strip = -1;
  return strip;
};

private:


  int _cfecSide;
  int _view;
  int _firstStrip;
  int _lastStrip;
  int _dir;  // direction of increasing strips
  int _firstPin[2];
};



#endif
