#ifndef T2_STRIP_DIGI_H
#define T2_STRIP_DIGI_H

/**
 * Class that describes the characteristica of a T2 strip
 * including position in row and columns ccordinates
 * and adc.
 *
 * Author: Erik Brücken / University of Helsinki
 * Email: brucken@cc.helsinki.fi
 * Updated 2008-05-27
 */

#include "boost/cstdint.hpp"

class T2StripDigi {

 private:

  int hitCounter;
  int adc;
  int stripNr;
  int row; 
  int col;

 public:

  T2StripDigi(int stripNr, int row, int col, int adc);
  T2StripDigi();
  ~T2StripDigi(){};
  void reset();

  int getHitCounter() const { return hitCounter; };
  int getAdc() const { return adc; };
  int getStripNr() const { return stripNr; };
  int getRow() const { return row; };
  int getCol() const { return col; };

}; // T2DigiStrip

#endif // T2_STRIP_DIGI_H
