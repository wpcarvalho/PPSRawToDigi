#ifndef T2_Detector_Hit_h
#define T2_Detector_Hit_h

/**
 * Class T2DetectorHit describes a strip or pad during 
 * the digitization process.
 * An almost straight copy from CSCDetectorHit
 * 
 * Author: Erik Brücken / University of Helsinki
 * email:  brucken@cc.helsinki.fi
 * Date:   2007-11-26
 */

#include <iosfwd>

class PSimHit;

class T2DetectorHit {

 public:
  
  T2DetectorHit(int element, float charge, int row, int col, float time, 
		const PSimHit * hitp) : 
    theElement(element), 
    theCharge(charge), 
    theRow(row), 
    theCol(col), 
    theTime(time), 
    theHitp(hitp) 
    {}
  
  int getElement() const { return theElement; }
  float getCharge() const { return theCharge; }
  int getRow() const { return theRow; }
  int getCol() const { return theCol; }
  float getTime() const { return theTime; }
  const PSimHit * getSimHit() const { return theHitp; }

  friend std::ostream & operator<<(std::ostream &, const T2DetectorHit &);

 private:
  
  int theElement; // strip or pad number
  float theCharge;
  int theRow; // the the row and column position
  int theCol;  
  float theTime; // start counting time at beam crossing
  const PSimHit * theHitp;

};

#endif // T2_Detector_Hit_h
