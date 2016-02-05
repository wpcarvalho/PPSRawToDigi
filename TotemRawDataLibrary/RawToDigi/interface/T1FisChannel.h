#ifndef _T1FisChannel_
#define _T1FisChannel_

//#include "DataFormats/T1DetId/interface/T1DetId.h"


using namespace std;

class T1FisChannel{

 public:
  T1FisChannel(){
    myDetId_.setLayer(0,0,0,0);
    Channel_=0;
    Type_=0;
  }

  T1FisChannel(T1DetId a,int b,int c){
    myDetId_=a;

    Channel_=c;
    Type_=b;
  }    

  ~T1FisChannel(){};

  void set_fisch(T1DetId a, int b,int c){
    myDetId_=a;

    Channel_=c;
    Type_=b;
  }

  void set_type(int type){
    Type_=type;
  }


  T1DetId DetId(){return myDetId_;}
  int Channel(){return Channel_;}
  int Type(){return Type_;}

 private:

  T1DetId myDetId_;
  int Channel_;
  int Type_;

};

#endif
