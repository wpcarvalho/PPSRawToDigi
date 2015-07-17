#ifndef T2_PAD_DIGI_H
#define T2_PAD_DIGI_H

/**
 * Class that describes the characteristica of a T2 pad
 * including position in row and columns ccordinates
 * and adc.
 *
 * author: Erik Brücken / University of Helsinki
 * email: serres@gmx.li
 * updated: 2008-05-27 
 */

#include "boost/cstdint.hpp"

class T2PadDigi {

 private:

  int hitCounter;
  int adc;
  int padNr;
  int row; 
  int col;

 public:

  T2PadDigi(int padNr, int row, int col, int adc);
  T2PadDigi();
  ~T2PadDigi(){};
  void reset();

  int getHitCounter() const { return hitCounter; };
  int getAdc() const { return adc; };
  int getPadNr() const { return padNr; };
  int getRow() const { return row; };
  int getCol() const { return col; };

}; // T2DigiPad

#endif // T2_PAD_DIGI_H
