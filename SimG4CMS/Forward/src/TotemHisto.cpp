#include "SimG4CMS/Forward/interface/TotemHisto.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <cmath>


TotemHisto::TotemHisto(std::string name) :
        file_name(name) {
    edm::LogInfo("TotemHisto") << "===>>>  Start booking user histograms with Root";
    ntuple = new TNtuple("ntuple", "Ntuple",
                         "Event:UnitID:Ptype:TrackID:ParentID:ELoss:PABS:vx:vy:vz:x:y:z:Px:Py:Pz:VPx:VPy:VPz");
    edm::LogInfo("TotemHisto") << "===>>> Done booking user histograms and ntuples ";
}

TotemHisto::~TotemHisto() {
    edm::LogInfo("TotemHisto") << "========================================================";
    edm::LogInfo("TotemHisto") << "=== TotemHisto: Start writing user histograms ===";
    const char *c_file_name;
    c_file_name = file_name.c_str();
    TFile rt_hf(c_file_name, "RECREATE");
    rt_hf.SetCompressionLevel(2);
    edm::LogInfo("TotemHisto") << "I have already created root file";
    TotemHisto::ntuple->Write();
    rt_hf.Close();
    edm::LogInfo("TotemHisto") << "TotemHisto: End writing user histograms ";
}

void TotemHisto::fillNtuple() {
    edm::LogInfo("TotemHisto") << "Saving ntuple to root file";
    rootvec[0] = evt;
    rootvec[1] = UID;
    rootvec[2] = Ptype;
    rootvec[3] = TID;
    rootvec[4] = PID;
    rootvec[5] = ELoss;
    rootvec[6] = PABS;
    rootvec[7] = vx;
    rootvec[8] = vy;
    rootvec[9] = vz;
    rootvec[10] = x;
    rootvec[11] = y;
    rootvec[12] = z;
    rootvec[13] = Px;
    rootvec[14] = Py;
    rootvec[15] = Pz;
    rootvec[16] = VPx;
    rootvec[17] = VPy;
    rootvec[18] = VPz;

    if (ntuple)
        ntuple->Fill(&(rootvec[0]));
}
