#ifndef __MyMITEvtSel_H__
#define __MyMITEvtSel_H__

#include "TObject.h"

class MyMITEvtSel : public TObject {

  public :
    MyMITEvtSel();
    virtual ~MyMITEvtSel();

    void Reset();
    void Print();

    double eHcalNeg;       //energy HCAL negative side
    double eHcalPos;       //energy HCAL positive side
    double eHfNeg;         //energy HF negative side
    double eHfPos;         //energy HF positive side
    double eHfNegTime;     //energy weighted HF time on negative side 
    double eHfPosTime;     //energy weighted HF time on positive side 
    double eCaNeg;         //energy CASTOR negative side
    double eCaPos;         //energy CASTOR positive side
    double eCaNegTime;     //energy weighted CASTOR time on negative side 
    double eCaPosTime;     //energy weighted CASTOR time on positive side 
    double eZdcNeg;        //energy ZDC negative side
    double eZdcPos;        //energy ZDC positive side
    double eZdcNegTime;    //energy weighted ZDC time on negative side 
    double eZdcPosTime;    //energy weighted ZDC time on positive side 
    int    ePxbHits;	   //number of pixel rechits in the three barrel layers
    int    ePxHits;	   //number of pixel rechits in all barrel and forward layers
    double eClusVtxQual;   //incompatibility of pixel cluster shapes with vertex (ratio)
    double eClusVtxDiff;   //incompatibility of pixel cluster shapes with vertex (difference)
    int    nHfNegHits;     //hf neg hits above threshold
    int    nHfPosHits;     //hf pos hits above threshold
    int    nHfTowersP;     //hf neg calo towers above threshold
    int    nHfTowersN;     //hf pos calo towers above threshold
    double sumEsubEpPos;   //sum E sub Ep for pos calo towers
    double sumEaddEpPos;   //sum E add Ep for pos calo towers
    double sumEsubEpNeg;   //sum E sub Ep for neg calo towers
    double sumEaddEpNeg;   //sum E add Ep for neg calo towers
    double sumHfEsubEpPos; //sum E sub Ep for pos hf calo towers
    double sumHfEaddEpPos; //sum E add Ep for pos hf calo towers
    double sumHfEsubEpNeg; //sum E sub Ep for neg hf calo towers
    double sumHfEaddEpNeg; //sum E add Ep for neg hf calo towers
    double eHPTrkFrac;     //fraction of high-purity tracks out of all with "loose" cuts


  private:

  ClassDef (MyMITEvtSel,1)
};

#endif

