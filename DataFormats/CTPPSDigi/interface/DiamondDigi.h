#ifndef CTPPSDigi_DiamondDigi_h
#define CTPPSDigi_DiamondDigi_h

/** \class DiamondDigi
 *
 * Digi for PPSTiming.
 *  
 *
 * \author Seyed Mohsen Etesami
 * March 2016
 */

#include <boost/cstdint.hpp>

class DiamondDigi{

public:
  
  DiamondDigi (int chid_,int ledgt_, int tedgt_, int threvolt,int mhit_, int hptdcerror_);
  // Default construction.
  
  DiamondDigi ();
  ~DiamondDigi() {};
  

  /// Digis are equal if they are have same chid, ledt and tedt, threshold voltage, multihit flag, hptdcerror flags
  bool operator==(const DiamondDigi& digi) const;


  /// Return digi values number
  int getCHID() const;
  int getLEDT() const;
  int getTEDT() const;
  int getTHREVOLT() const;
  int getMHIT() const;
  int getHPTDCERROR() const;


  /// Set digi values
  void setCHID (int chid_);
  void setLEDT(int ledgt_);  
  void setTEDT(int tedgt_);
  void setTHREVOLT(int threvolt_);
  void setMHIT(int mhit_);
  void setHPTDCERROR(int hptdcerror_);


  /// Print content of digi
 //   void print() const;

private:

  int chid;
  int ledgt;
  int tedgt;
  int threvolt;
  int mhit;
  int hptdcerror;
};


inline bool operator< (const DiamondDigi& one, const DiamondDigi& other)
{
   return one.getCHID() < other.getCHID();

}


  #include<iostream>
  inline std::ostream & operator<<(std::ostream & o, const DiamondDigi& digi) {
return o   << " " << digi.getCHID()
           << " " << digi.getLEDT()
	   << " " << digi.getTEDT()
           << " " << digi.getTHREVOLT()
             << " " << digi.getMHIT()
           << " " << digi.getHPTDCERROR();

}
#endif

