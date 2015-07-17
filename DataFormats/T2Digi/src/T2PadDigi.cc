#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include <iostream>

/**
 * Constructor
 */

T2PadDigi::T2PadDigi(int padNr, int row, int col, int adc) {
  
  this->padNr = padNr;
  this->row = row;
  this->col = col;
  this->adc = adc;
  
} // T2PadDigi

/**
 * Default constructor
 */

T2PadDigi::T2PadDigi() {
  
  padNr = 0;
  row = 0;
  col = 0;
  adc = 0;
  
} // T2PadDigi

/**
 * Resets the pad class.
 */

void T2PadDigi::reset() {

  hitCounter = 0;   // Set the hit counter to zero.

} // reset
