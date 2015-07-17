#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include <iostream>

/**
 * Constructor
 */

T2StripDigi::T2StripDigi(int stripNr, int row, int col, int adc) {
  
  this->stripNr = stripNr;
  this->row = row;
  this->col = col;
  this->adc = adc;
  
} // T2StripDigi

/**
 * Default constructor
 */

T2StripDigi::T2StripDigi() {
  
  stripNr = 0;
  row = 0;
  col = 0;
  adc = 0;
  
} // T2StripDigi

/**
 * Resets the strip class.
 */

void T2StripDigi::reset() {

  hitCounter = 0;   // Set the hit counter to zero.

} // reset
