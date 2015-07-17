/****************************************************************************
 *
 * Author: Fabrizio Ferro - INFN Genova
 *
 ****************************************************************************/

#ifndef DataFormatsTotemL1TriggerT1TriggerBits_h
#define DataFormatsTotemL1TriggerT1TriggerBits_h

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iosfwd>
#include <iostream>
#include <cstdlib>
#include <bitset>

class T1TriggerBits
{
 public:
 
  explicit T1TriggerBits( int arm, std::bitset<32> bs) { _arm=arm; setBS(bs); };

 T1TriggerBits(){ };

  void setArm( int );

  inline int getArm() const {return _arm;};

  inline std::bitset<32> getBS() const {
    std::bitset<32> res;
    for( int i = 0 ; i < 32 ; ++i ){
     res[i] = _bs[i];
   }
    return res;
 };


  bool getBit(int i){
    return _bs[i];
  }

  void setBS( std::bitset<32> bs ){
    for( int i = 0 ; i < 32 ; ++i ){
     _bs[i] = bs[i];
   }
 };

 private:
  int _arm;
  bool _bs[32];

}; 




#endif  //DataFormatsTotemL1TriggerT1TriggerBits_h
