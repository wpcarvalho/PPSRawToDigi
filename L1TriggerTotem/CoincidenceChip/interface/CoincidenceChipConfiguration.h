// -*- C++ -*-
//
// Package:    CoincidenceChip
// Class:      CoincidenceChipConfiguration
//
// Author: Leszek Grzanka

#ifndef _L1TriggerTotemCoincidenceChipCoincidenceChipConfiguration_H_
#define _L1TriggerTotemCoincidenceChipCoincidenceChipConfiguration_H_

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <stdio.h>
#include <string>
#include <bitset>
#include <iostream>

#include "TString.h"

class CoincidenceChipConfiguration {

  public:
    CoincidenceChipConfiguration();
    ~CoincidenceChipConfiguration();

    void configure(const edm::ParameterSet& iConfig);

	/******** control register  **********/
    void setControlReg(std::bitset<24> controlRegisters);
    std::bitset<24> getControlReg() const;

	/******** control register 1 **********/
	void setControlRegister1(unsigned long n);
	unsigned short getControlRegister1() const;

	/******** control register 2 **********/
	void setControlRegister2(unsigned long n);
	unsigned short getControlRegister2() const;

	/******** control register 3 **********/
	void setControlRegister3(unsigned long n);
	unsigned short getControlRegister3() const;


    /******** logic configuration *********/
    void useLogicWithWrongNP();
    void useLogicWithCorrectNP();
    bool getLogicWithWrongNPFlag() const;

    /************  V **************/
    void setV(unsigned short v_);
    unsigned short getV() const;

    /************  NP **************/
    void setNP(unsigned short np_);
    unsigned short getNP() const;

    /************  OV **************/
    void setOV(unsigned short ov_);
    unsigned short getOV() const;

    /************  W **************/
    void setW(unsigned short w_);
    unsigned short getW() const;

    /************  Z **************/
    void setZ(unsigned short z_);
    unsigned short getZ() const;

    /************ O2 **************/
    void setO2(unsigned short o2_);
    unsigned short getO2() const;

    /************ LI **************/
    void setLI(unsigned short li_);
    unsigned short getLI() const;

    /************ LO **************/
    void setLO(unsigned short lo_);
    unsigned short getLO() const;

    /************ AO **************/
    void setAO(unsigned short ao_);
    unsigned short getAO() const;

    //************ summary **************/
    //std::string summary() const;

    /************ summary **************/
    std::string configSummary() const;

	//******** Print configuration of CC **********/
	//void PrintConfig() const;

	void setVerbose(const bool verbose ){ verbose_ = verbose; };
	bool getVerbose() const{ return verbose_; };

  private:
    std::bitset<24> controlRegisters_;
    unsigned int useLogicWithWrongNP_;
	bool verbose_;

	void CheckControlRegister(unsigned short regNumber, unsigned int n) const;

};

#endif
