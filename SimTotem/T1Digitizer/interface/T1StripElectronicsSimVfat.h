/*
 * modified by Marcin Borratynski (mborratynski@gmail.com)
 */

#ifndef T1_STRIP_ELECTRONICS_SIM_VFAT_H
#define T1_STRIP_ELECTRONICS_SIM_VFAT_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "SimMuon/CSCDigitizer/src/CSCBaseElectronicsSim.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfatCollection.h"
#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"
#include "SimTotem/T1Digitizer/interface/T1DeadChannelManager.h"

// declarations
/*
  class CSCLayer;
  class CSCDetectorHit;
  class CSCWireDigi;
  class CSCAnalogSignal;
*/

class T1StripElectronicsSimVfat 
{
 public:
  /// configurable parameters
  T1StripElectronicsSimVfat(const edm::ParameterSet &p);



  //  void setFraction(float newFraction)  {theFraction = newFraction;};

  void fillDigis(int,int,std::vector<T1DetectorHit>,T1DigiVfatCollection & digis, T1DeadChannelManager & channelManager);

 private:
  // helper functions
  //  void init();
  /// initialization for each layer
  //  virtual void initParameters();

  // will return wire group, given wire.
  //  virtual int readoutElement(int element) const;

  //  float calculateAmpResponse(float t) const;
 
  //  virtual float signalDelay(int element, float pos) const;
  //  virtual float timeOfFlightCalibration(int wireGroup) const;

  /// we code strip indices from 1-80, and wire indices start at 100
  //  virtual int channelIndex(int channel) const {return channel+100;}

  // member data
  // the fractional discriminator returns the time when the signal
  // reaches this fraction of its maximum
  //  float theFraction;
  float theStripNoise;
  float theStripThreshold1;
  float theStripThreshold2;

  float _threshold1;
  float _threshold2;
  float _rumore;
  //  float theTimingCalibrationError; // in ns
};

#endif
