/****************************************************************************
 *
 * Author: Fabrizio Ferro - INFN Genova
 *
 * Copied by Fredrik Oljemark 02.02.2009
 *
 ****************************************************************************/

#ifndef DataFormatsTotemL1TriggerT2TriggerBits_h
#define DataFormatsTotemL1TriggerT2TriggerBits_h

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"

#include <iosfwd>
#include <iostream>
#include <cstdlib>
#include <bitset>

class T2TriggerBits
{
 public:


  /// Construct from a packed id. It is required that the Detector part of
  /// id is Totem and the SubDet part is RP, otherwise an exception is thrown.
  //  explicit RPCCBits(RPCCIdRaw id, std::bitset<8> bs) : id_(id){ setBS(bs); };

  explicit T2TriggerBits(uint32_t id,int row, std::bitset<8> bs) : id_(id){setBS(bs); setSector(row);};

  /// Construct from fully qualified identifier.
  explicit T2TriggerBits( T2DetId id, int row, std::bitset<8> bs) { id_ =(uint32_t) id; setBS(bs); setSector(row);};

  T2TriggerBits() : id_(0){ };

  // void setId( T2DetId id ){id_ = id;};
  void setId( T2DetId id ){id_ = (uint32_t) id;};

  void setSector(int row){row_ = row;};
 

  // T2TriggerBits(){ };

  //  void setArm( int );

  inline uint32_t getId() const {return id_;};

  inline int getSector() const {return row_;};

  inline std::bitset<8> getBS() const {
    std::bitset<8> res;
    for( int i = 0 ; i < 8 ; ++i ){
     res[i] = _bs[i];
   }
    return res;
 };


  bool getBit(int i){
    return _bs[i];
  }

  void setBS( std::bitset<8> bs ){
    for( int i = 0 ; i < 8 ; ++i ){
     _bs[i] = bs[i];
   }
 };

 private:
  //  int _arm;
  uint32_t id_;
  bool _bs[8];
  int row_;
}; 




#endif  //DataFormatsTotemL1TriggerT2TriggerBits_h
