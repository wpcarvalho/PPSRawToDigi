/** \file
 * 
 *
 * \author Seyed Mohsen Etesami
 */


#include <DataFormats/CTPPSDigi/interface/DiamondDigi.h>


using namespace std;




   DiamondDigi::DiamondDigi (int chid_, int ledgt_, int tedgt_, int threvolt_, int mhit_, int hptdcerror_) 
  {
  chid=chid_;   
  ledgt=ledgt_;
  tedgt=tedgt_;
  threvolt=threvolt_;
  mhit=mhit_;
  hptdcerror=hptdcerror_;

  }




    DiamondDigi::DiamondDigi ()
   {
  chid=0;
  ledgt=0;
  tedgt=0;
  threvolt=0;
  mhit=0;
  hptdcerror=0;
   }


  // Comparison
  bool DiamondDigi::operator == (const DiamondDigi& digi) const {if (chid!=digi.getCHID()|| ledgt!=digi.getLEDT()||tedgt != digi.getTEDT()||threvolt != digi.getTHREVOLT()|| mhit!=digi.getMHIT()|| hptdcerror!=digi.getHPTDCERROR() ) return false;
  return true; 
  } 

   // Getters
  int DiamondDigi::getCHID() const{ return chid;}
  int DiamondDigi::getLEDT() const{ return ledgt;}
  int DiamondDigi::getTEDT() const{ return tedgt;}
   int DiamondDigi::getTHREVOLT() const{ return threvolt;}
   int DiamondDigi::getMHIT() const{ return mhit;}
   int DiamondDigi::getHPTDCERROR() const{ return hptdcerror;}


  /// Set digi values
void DiamondDigi::setCHID(int chid_)
   {chid=chid_;}
  void DiamondDigi::setLEDT(int ledgt_)
   {ledgt=ledgt_;}
  void DiamondDigi::setTEDT(int tedgt_)
   {tedgt=tedgt_;}  
void DiamondDigi::setTHREVOLT(int threvolt_)
   {threvolt=threvolt_;}
void DiamondDigi::setMHIT(int mhit_)
   {mhit=mhit_;}
void DiamondDigi::setHPTDCERROR(int hptdcerror_)
   {hptdcerror=hptdcerror_;}


