#include "L1TriggerTotem/CoincidenceChip/interface/CoincidenceChip.h"
#include <bitset>
#include <iostream>

/**
 * bitset referece : http://www.cplusplus.com/reference/stl/bitset/
 */

/*
void CoincidenceChip::configure(CoincidenceChipConfiguration _config) {
	// do we need this method?
	assert(0);
}
*/

unsigned short CoincidenceChip::NumberOfDetectorPlanes() const{ 
	if (getNP() == 0x1) 
		return  5;
	 else 
		return 10;
}

CoincidenceChip::OutputBits CoincidenceChip::process(InputBits inputBits) const{
  InputBits resOr1;
  std::bitset<16> resOr2;
  std::bitset<16> resAO;
  std::bitset<16> resAO2;
  std::bitset<16> resZ;
  std::bitset<16> resVOutOfNP;
  std::bitset<16> resWOutOfNP;

  if( getVerbose()) edm::LogInfo("CClogic") << "inputBits";
  if( getVerbose()) edm::LogInfo("CClogic") << "0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF";
  if( getVerbose()) edm::LogInfo("CClogic") << inputBits;

#if 0
  !!!! THE LI REGISTER IS RESPONSIBLE FOR SOMETHING ELSE, NOT FOR 
  !!!! CHANGING POLARITY OF INPUTS BITS (see bitset::flip() method)
  // Stage 1 - Inverting inputs, according to LI register
  // LI == 1 => invert inputs

  if (getLI() == 1) {
    inputBits.flip();  dont use this method !!!
    if( getVerbose()) edm::LogInfo("CClogic") << "Inverting inputs";
    if( getVerbose()) edm::LogInfo("CClogic") << inputBits;
  }
#endif

  // Step 2 - Reduction from 80 to 16 bits: "V out of NP" , "Or 1" and "W out of NP" blocks
  // NP == 1 => We have 5 planes ( NP = 5 )
  if( getVerbose()) edm::LogInfo("CClogic") << NumberOfDetectorPlanes() << " planes";

  // "V out of NP" block
  if (getV() <= NumberOfDetectorPlanes()) {
    // NP == 5
    if (NumberOfDetectorPlanes() == 5) {
      // for each coordinate:
      for (unsigned short i = 0; i < 16; i++) {
        // if at least V out of NP bits are on
        if (inputBits[i] + inputBits[i + 16] + inputBits[i + 32] + inputBits[i + 48] + inputBits[i + 64] >= getV()) {
          // then coincidence signal is generated
          resVOutOfNP.set(i);
        }
      }
      // NP == 10
    } else {
      // for each coordinate:
      for (unsigned short i = 0; i < 8; i++) {
        // if at least V out of NP bits are on
        if (inputBits[i] + inputBits[i + 8] + inputBits[i + 16] + inputBits[i + 24] + inputBits[i + 32] + inputBits[i + 40] + inputBits[i + 48] + inputBits[i + 56] + inputBits[i + 64] + inputBits[i + 72] >= getV()) {
          // then coincidence signal is generated
          resVOutOfNP.set(i);
          resVOutOfNP.set(8+i);
        }
      }
    }
  } else {
    if( getVerbose()) edm::LogInfo("CClogic") << "V > " << NumberOfDetectorPlanes() << " , resetting output";
    resVOutOfNP.reset();
  }
  if( getVerbose()) edm::LogInfo("CClogic") << "V = " << getV();
  if( getVerbose()) edm::LogInfo("CClogic") << resVOutOfNP;

  // "OR1" block
  unsigned short np_current = NumberOfDetectorPlanes();
  if( getLogicWithWrongNPFlag() == true ){
     np_current = 15 - NumberOfDetectorPlanes();
     // 15 - 5 = 10  ...  5 -> 10
     // 15 - 10 = 5  ... 10 -> 5
  } 
  resOr1 = inputBits;

  // For every input signal
  for (short i = 0; i < 80; i++) {
    // For every neighbor in distance OV
    for (short j = 0; j <= getOV(); j++) {
      // Do the OR operation of input signal with its neighbors
      // Check boundary conditions
      // Take into account left and right neighbors

      // 80 / 5 = 16 , 80 / 10 = 8
      if ((i % (80 / np_current)) + j < 16) resOr1[i] = resOr1[i] | inputBits[i + j];
      if ((i % (80 / np_current)) - j >= 0) resOr1[i] = resOr1[i] | inputBits[i - j];
    }
  }

  if( getVerbose()) edm::LogInfo("CClogic") << "OV = " << getOV() << " , result:";
  if( getVerbose()) edm::LogInfo("CClogic") << resOr1;

  // "W out of NP" block
  // Similar as "V out of NP" block, but takes as input data coming from "OR1" block
  if (getW() <= NumberOfDetectorPlanes()) {
    // NP == 5
    if (NumberOfDetectorPlanes() == 5) {
      // for each coordinate:
      for (unsigned short i = 0; i < 16; i++) {
        // if at least W out of NP bits are on
        if (resOr1[i] + resOr1[i + 16] + resOr1[i + 32] + resOr1[i + 48] + resOr1[i + 64] >= getW()) {
          // then coincidence signal is generated
          resWOutOfNP.set(i);
        }
      }
      // NP == 10
    } else {
      // for each coordinate:
      for (unsigned short i = 0; i < 8; i++) {
        // if at least W out of NP bits are on
        if (resOr1[i] + resOr1[i + 8] + resOr1[i + 16] + resOr1[i + 24] + resOr1[i + 32] + resOr1[i + 40] + resOr1[i + 48] + resOr1[i + 56] + resOr1[i + 64] + resOr1[i + 72] >= getW()) {
          // then coincidence signal is generated
          resWOutOfNP.set(i);
          resWOutOfNP.set(8+i);
        }
      }
    }
 

  } else {
    if( getVerbose()) edm::LogInfo("CClogic") << "W > " << NumberOfDetectorPlanes() << ", resetting output";
    resWOutOfNP.reset();
  }
  if( getVerbose()) edm::LogInfo("CClogic") << "W = " << getW() << " , result:";
  if( getVerbose()) edm::LogInfo("CClogic") << resWOutOfNP;

  // Step 3 - "And/Or" block
  // Combine result of "V out of NP" and "W out of NP" blocks
  // Configuration stored in AO register
  if (getAO() == 0) {
    resAO = resVOutOfNP & resWOutOfNP;
  } else if (getAO() == 1) {
    resAO = resVOutOfNP | resWOutOfNP;
  } else if (getAO() == 2) {
    resAO = resVOutOfNP;
  } else if (getAO() == 3) {
    resAO = resWOutOfNP;
  } else {
    if( getVerbose())  edm::LogInfo("CClogic") << "AO > 3, resetting output";
    resAO.reset();
  }
  if( getVerbose()) edm::LogInfo("CClogic") << "A0 = " << getAO();
  if( getVerbose()) edm::LogInfo("CClogic") << resAO;

  // Step 4 - "Z out of 16 or 8" block
  // Configuration stored in AO register
  // If more than Z groups are on, then results are masked
  // As result of this block we have mask

  // np/5 : 10/5 = 2 , 5/5 = 1
  np_current = NumberOfDetectorPlanes();
  if( getLogicWithWrongNPFlag() == true ){
     np_current = 15 - NumberOfDetectorPlanes();
     // 15 - 5 = 10  ...  5 -> 10
     // 15 - 10 = 5  ... 10 -> 5
  } 
  // (5/np_current) is equal to 1 or 2; we need to divide resAO.count() by 2 in same cases, 
  // such as in the case of T2 (10 planes and 8 sectors) and correct logic (np_current is also 10). 
  if((unsigned short) resAO.count()*(5/np_current) > (getZ() + 1)) { 
    resZ.set(); // resZ = 0xF
  } else {
    resZ.reset(); // resZ = 0x0
  }
  if( getVerbose()) edm::LogInfo("CClogic") << "Z = " << getZ();
  if( getVerbose()) edm::LogInfo("CClogic") << resZ;

  // Step 5 - "And/Or 2" block
  // Configuration stored in LO register
  // Combine mask coming from "Z out of 16 or 8" block and
  // output of "And/Or" block
  if (getLO() == 0) {
    resAO2 = resAO & (~resZ);
  } else if (getLO() == 1) {
    resAO2 = resAO & resZ;
  } else if (getLO() == 2) {
    resAO2 = ~(resAO & (~resZ));
  } else if (getLO() == 3) {
    resAO2 = ~(resAO & resZ);
  } else {
    if( getVerbose()) edm::LogInfo("CClogic") << "LO > 3, resetting output";
    resAO2.reset();
  }
  if( getVerbose()) edm::LogInfo("CClogic") << "LO = " << getLO();
  if( getVerbose()) edm::LogInfo("CClogic") << resAO2;

  // Step 6 - "Or 2" block
  // Responsible for grouping outputs
  resOr2.reset();
  if (getO2() == 0) {
    resOr2 = resAO2;
  } else if (getO2() == 1) {
    resOr2[0] = resAO2[0] | resAO2[1];
    resOr2[2] = resAO2[2] | resAO2[3];
    resOr2[4] = resAO2[4] | resAO2[5];
    resOr2[6] = resAO2[6] | resAO2[7];
    resOr2[8] = resAO2[8] | resAO2[9];
    resOr2[10] = resAO2[10] | resAO2[11];
    resOr2[12] = resAO2[12] | resAO2[13];
    resOr2[14] = resAO2[14] | resAO2[15];
  } else if (getO2() == 2) {
    resOr2[0] = resAO2[0] | resAO2[1] | resAO2[2] | resAO[3];
    resOr2[4] = resAO2[4] | resAO2[5] | resAO2[6] | resAO[7];
    resOr2[8] = resAO2[8] | resAO2[9] | resAO2[10] | resAO[11];
    resOr2[12] = resAO2[12] | resAO2[13] | resAO2[14] | resAO[15];
  } else if (getO2() == 3) {
    resOr2[0] = resAO2[0] | resAO2[1] | resAO2[2] | resAO[3] | resAO2[4] | resAO2[5] | resAO2[6] | resAO[7];
    resOr2[8] = resAO2[8] | resAO2[9] | resAO2[10] | resAO[11] | resAO2[12] | resAO2[13] | resAO2[14] | resAO[15];
  } else if (getO2() == 4) {
    resOr2[0] = resAO2.any();
  } else {
    if( getVerbose()) edm::LogInfo("CClogic") << "O2 > 4, resetting output";
    resOr2.reset();
  }

  return resOr2;
}

