#ifndef __MyCaloJet_H__
#define __MyCaloJet_H__

#include "MyJet.h"

class MyCaloJet : public MyJet {

  public :
    MyCaloJet();
    ~MyCaloJet();

    virtual void Reset();
    virtual void Print();
    
    //-- estimators (Calojet only)
    Double_t fem;        //-- emEnergyFraction() (for PFJet: fem_ch + fem_ne)
    Double_t eem_EB;     //-- emEnergyInEB()
    Double_t eem_EE;     //-- emEnergyInEE()
    Double_t eem_HF;     //-- emEnergyInHF()

    Double_t fhad;      //-- energyFractionHadronic()  (for PFjet: fhad_ch + fhad_ne)
    Double_t ehad_HB;   //-- hadEnergyInHB()
    Double_t ehad_HE;   //-- hadEnergyInHE()
    Double_t ehad_HF;   //-- hadEnergyInHF()
    Double_t ehad_HO;   //-- hadEnergyInHO()

    Int_t n60;  //-- n60()
    Int_t n90;  //-- n90()

    Double_t emax_ecal; //-- maxEInEmTowers()
    Double_t emax_hcal; //-- maxEInHadTowers()

    Int_t n90hits;
    Double_t HPD;
    Double_t RBX;
    Double_t sigma_eta;
    Double_t sigma_phi;

  private:

  ClassDef (MyCaloJet,1)
};

#endif
