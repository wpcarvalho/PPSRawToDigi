{
TFile f("SDValid.root");
f.GetListOfKeys().Print();
TH1F *hSimEta = (TH1F*)f.Get("SimEta");
TH1F *hSimPT = (TH1F*)f.Get("SimPT");
TH1F *hSimMomentum=(TH1F *)f.Get("SimMomentum");
TH1F *hSimMomentumLog=(TH1F *)f.Get("SimMomentumLog");
TH1F *hSimEnergyloss=(TH1F *)f.Get("SimEnergy loss");
TH1F *hSimType=(TH1F *)f.Get("SimType");
TH1F *hSimUnitID=(TH1F *)f.Get("SimUnitID");

TH1F *hClusterWidth=(TH1F *)f.Get("Cluster Width");

TH1F *hWireCh0Arm0Pl0=(TH1F *)f.Get("Wire Ch0 Arm0 Pl0");
TH1F *hWireCh0Arm0Pl1=(TH1F *)f.Get("Wire Ch0 Arm0 Pl1");
TH1F *hWireCh0Arm0Pl2=(TH1F *)f.Get("Wire Ch0 Arm0 Pl2");
TH1F *hWireCh0Arm0Pl3=(TH1F *)f.Get("Wire Ch0 Arm0 Pl3");
TH1F *hWireCh0Arm0Pl4=(TH1F *)f.Get("Wire Ch0 Arm0 Pl4");
TH1F *hWireCh0Arm1Pl0=(TH1F *)f.Get("Wire Ch0 Arm1 Pl0");
TH1F *hWireCh0Arm1Pl1=(TH1F *)f.Get("Wire Ch0 Arm1 Pl1");
TH1F *hWireCh0Arm1Pl2=(TH1F *)f.Get("Wire Ch0 Arm1 Pl2");
TH1F *hWireCh0Arm1Pl3=(TH1F *)f.Get("Wire Ch0 Arm1 Pl3");
TH1F *hWireCh0Arm1Pl4=(TH1F *)f.Get("Wire Ch0 Arm1 Pl4");

TH1F *hStripACh0Arm0Pl0=(TH1F *)f.Get("StripA Ch0 Arm0 Pl0");
TH1F *hStripACh0Arm0Pl1=(TH1F *)f.Get("StripA Ch0 Arm0 Pl1");
TH1F *hStripACh0Arm0Pl2=(TH1F *)f.Get("StripA Ch0 Arm0 Pl2");
TH1F *hStripACh0Arm0Pl3=(TH1F *)f.Get("StripA Ch0 Arm0 Pl3");
TH1F *hStripACh0Arm0Pl4=(TH1F *)f.Get("StripA Ch0 Arm0 Pl4");
TH1F *hStripACh0Arm1Pl0=(TH1F *)f.Get("StripA Ch0 Arm1 Pl0");
TH1F *hStripACh0Arm1Pl1=(TH1F *)f.Get("StripA Ch0 Arm1 Pl1");
TH1F *hStripACh0Arm1Pl2=(TH1F *)f.Get("StripA Ch0 Arm1 Pl2");
TH1F *hStripACh0Arm1Pl3=(TH1F *)f.Get("StripA Ch0 Arm1 Pl3");
TH1F *hStripACh0Arm1Pl4=(TH1F *)f.Get("StripA Ch0 Arm1 Pl4");

TH1F *hStripBCh0Arm0Pl0=(TH1F *)f.Get("StripB Ch0 Arm0 Pl0");
TH1F *hStripBCh0Arm0Pl1=(TH1F *)f.Get("StripB Ch0 Arm0 Pl1");
TH1F *hStripBCh0Arm0Pl2=(TH1F *)f.Get("StripB Ch0 Arm0 Pl2");
TH1F *hStripBCh0Arm0Pl3=(TH1F *)f.Get("StripB Ch0 Arm0 Pl3");
TH1F *hStripBCh0Arm0Pl4=(TH1F *)f.Get("StripB Ch0 Arm0 Pl4");
TH1F *hStripBCh0Arm1Pl0=(TH1F *)f.Get("StripB Ch0 Arm1 Pl0");
TH1F *hStripBCh0Arm1Pl1=(TH1F *)f.Get("StripB Ch0 Arm1 Pl1");
TH1F *hStripBCh0Arm1Pl2=(TH1F *)f.Get("StripB Ch0 Arm1 Pl2");
TH1F *hStripBCh0Arm1Pl3=(TH1F *)f.Get("StripB Ch0 Arm1 Pl3");
TH1F *hStripBCh0Arm1Pl4=(TH1F *)f.Get("StripB Ch0 Arm1 Pl4");

TH1F *hWireCh1Arm0Pl0=(TH1F *)f.Get("Wire Ch1 Arm0 Pl0");
TH1F *hWireCh1Arm0Pl1=(TH1F *)f.Get("Wire Ch1 Arm0 Pl1");
TH1F *hWireCh1Arm0Pl2=(TH1F *)f.Get("Wire Ch1 Arm0 Pl2");
TH1F *hWireCh1Arm0Pl3=(TH1F *)f.Get("Wire Ch1 Arm0 Pl3");
TH1F *hWireCh1Arm0Pl4=(TH1F *)f.Get("Wire Ch1 Arm0 Pl4");
TH1F *hWireCh1Arm1Pl0=(TH1F *)f.Get("Wire Ch1 Arm1 Pl0");
TH1F *hWireCh1Arm1Pl1=(TH1F *)f.Get("Wire Ch1 Arm1 Pl1");
TH1F *hWireCh1Arm1Pl2=(TH1F *)f.Get("Wire Ch1 Arm1 Pl2");
TH1F *hWireCh1Arm1Pl3=(TH1F *)f.Get("Wire Ch1 Arm1 Pl3");
TH1F *hWireCh1Arm1Pl4=(TH1F *)f.Get("Wire Ch1 Arm1 Pl4");

TH1F *hStripACh1Arm0Pl0=(TH1F *)f.Get("StripA Ch1 Arm0 Pl0");
TH1F *hStripACh1Arm0Pl1=(TH1F *)f.Get("StripA Ch1 Arm0 Pl1");
TH1F *hStripACh1Arm0Pl2=(TH1F *)f.Get("StripA Ch1 Arm0 Pl2");
TH1F *hStripACh1Arm0Pl3=(TH1F *)f.Get("StripA Ch1 Arm0 Pl3");
TH1F *hStripACh1Arm0Pl4=(TH1F *)f.Get("StripA Ch1 Arm0 Pl4");
TH1F *hStripACh1Arm1Pl0=(TH1F *)f.Get("StripA Ch1 Arm1 Pl0");
TH1F *hStripACh1Arm1Pl1=(TH1F *)f.Get("StripA Ch1 Arm1 Pl1");
TH1F *hStripACh1Arm1Pl2=(TH1F *)f.Get("StripA Ch1 Arm1 Pl2");
TH1F *hStripACh1Arm1Pl3=(TH1F *)f.Get("StripA Ch1 Arm1 Pl3");
TH1F *hStripACh1Arm1Pl4=(TH1F *)f.Get("StripA Ch1 Arm1 Pl4");

TH1F *hStripBCh1Arm0Pl0=(TH1F *)f.Get("StripB Ch1 Arm0 Pl0");
TH1F *hStripBCh1Arm0Pl1=(TH1F *)f.Get("StripB Ch1 Arm0 Pl1");
TH1F *hStripBCh1Arm0Pl2=(TH1F *)f.Get("StripB Ch1 Arm0 Pl2");
TH1F *hStripBCh1Arm0Pl3=(TH1F *)f.Get("StripB Ch1 Arm0 Pl3");
TH1F *hStripBCh1Arm0Pl4=(TH1F *)f.Get("StripB Ch1 Arm0 Pl4");
TH1F *hStripBCh1Arm1Pl0=(TH1F *)f.Get("StripB Ch1 Arm1 Pl0");
TH1F *hStripBCh1Arm1Pl1=(TH1F *)f.Get("StripB Ch1 Arm1 Pl1");
TH1F *hStripBCh1Arm1Pl2=(TH1F *)f.Get("StripB Ch1 Arm1 Pl2");
TH1F *hStripBCh1Arm1Pl3=(TH1F *)f.Get("StripB Ch1 Arm1 Pl3");
TH1F *hStripBCh1Arm1Pl4=(TH1F *)f.Get("StripB Ch1 Arm1 Pl4");

TH1F *hWireCh2Arm0Pl0=(TH1F *)f.Get("Wire Ch2 Arm0 Pl0");
TH1F *hWireCh2Arm0Pl1=(TH1F *)f.Get("Wire Ch2 Arm0 Pl1");
TH1F *hWireCh2Arm0Pl2=(TH1F *)f.Get("Wire Ch2 Arm0 Pl2");
TH1F *hWireCh2Arm0Pl3=(TH1F *)f.Get("Wire Ch2 Arm0 Pl3");
TH1F *hWireCh2Arm0Pl4=(TH1F *)f.Get("Wire Ch2 Arm0 Pl4");
TH1F *hWireCh2Arm1Pl0=(TH1F *)f.Get("Wire Ch2 Arm1 Pl0");
TH1F *hWireCh2Arm1Pl1=(TH1F *)f.Get("Wire Ch2 Arm1 Pl1");
TH1F *hWireCh2Arm1Pl2=(TH1F *)f.Get("Wire Ch2 Arm1 Pl2");
TH1F *hWireCh2Arm1Pl3=(TH1F *)f.Get("Wire Ch2 Arm1 Pl3");
TH1F *hWireCh2Arm1Pl4=(TH1F *)f.Get("Wire Ch2 Arm1 Pl4");

TH1F *hStripACh2Arm0Pl0=(TH1F *)f.Get("StripA Ch2 Arm0 Pl0");
TH1F *hStripACh2Arm0Pl1=(TH1F *)f.Get("StripA Ch2 Arm0 Pl1");
TH1F *hStripACh2Arm0Pl2=(TH1F *)f.Get("StripA Ch2 Arm0 Pl2");
TH1F *hStripACh2Arm0Pl3=(TH1F *)f.Get("StripA Ch2 Arm0 Pl3");
TH1F *hStripACh2Arm0Pl4=(TH1F *)f.Get("StripA Ch2 Arm0 Pl4");
TH1F *hStripACh2Arm1Pl0=(TH1F *)f.Get("StripA Ch2 Arm1 Pl0");
TH1F *hStripACh2Arm1Pl1=(TH1F *)f.Get("StripA Ch2 Arm1 Pl1");
TH1F *hStripACh2Arm1Pl2=(TH1F *)f.Get("StripA Ch2 Arm1 Pl2");
TH1F *hStripACh2Arm1Pl3=(TH1F *)f.Get("StripA Ch2 Arm1 Pl3");
TH1F *hStripACh2Arm1Pl4=(TH1F *)f.Get("StripA Ch2 Arm1 Pl4");

TH1F *hStripBCh2Arm0Pl0=(TH1F *)f.Get("StripB Ch2 Arm0 Pl0");
TH1F *hStripBCh2Arm0Pl1=(TH1F *)f.Get("StripB Ch2 Arm0 Pl1");
TH1F *hStripBCh2Arm0Pl2=(TH1F *)f.Get("StripB Ch2 Arm0 Pl2");
TH1F *hStripBCh2Arm0Pl3=(TH1F *)f.Get("StripB Ch2 Arm0 Pl3");
TH1F *hStripBCh2Arm0Pl4=(TH1F *)f.Get("StripB Ch2 Arm0 Pl4");
TH1F *hStripBCh2Arm1Pl0=(TH1F *)f.Get("StripB Ch2 Arm1 Pl0");
TH1F *hStripBCh2Arm1Pl1=(TH1F *)f.Get("StripB Ch2 Arm1 Pl1");
TH1F *hStripBCh2Arm1Pl2=(TH1F *)f.Get("StripB Ch2 Arm1 Pl2");
TH1F *hStripBCh2Arm1Pl3=(TH1F *)f.Get("StripB Ch2 Arm1 Pl3");
TH1F *hStripBCh2Arm1Pl4=(TH1F *)f.Get("StripB Ch2 Arm1 Pl4");

TH1F *hWireCh3Arm0Pl0=(TH1F *)f.Get("Wire Ch3 Arm0 Pl0");
TH1F *hWireCh3Arm0Pl1=(TH1F *)f.Get("Wire Ch3 Arm0 Pl1");
TH1F *hWireCh3Arm0Pl2=(TH1F *)f.Get("Wire Ch3 Arm0 Pl2");
TH1F *hWireCh3Arm0Pl3=(TH1F *)f.Get("Wire Ch3 Arm0 Pl3");
TH1F *hWireCh3Arm0Pl4=(TH1F *)f.Get("Wire Ch3 Arm0 Pl4");
TH1F *hWireCh3Arm1Pl0=(TH1F *)f.Get("Wire Ch3 Arm1 Pl0");
TH1F *hWireCh3Arm1Pl1=(TH1F *)f.Get("Wire Ch3 Arm1 Pl1");
TH1F *hWireCh3Arm1Pl2=(TH1F *)f.Get("Wire Ch3 Arm1 Pl2");
TH1F *hWireCh3Arm1Pl3=(TH1F *)f.Get("Wire Ch3 Arm1 Pl3");
TH1F *hWireCh3Arm1Pl4=(TH1F *)f.Get("Wire Ch3 Arm1 Pl4");

TH1F *hStripACh3Arm0Pl0=(TH1F *)f.Get("StripA Ch3 Arm0 Pl0");
TH1F *hStripACh3Arm0Pl1=(TH1F *)f.Get("StripA Ch3 Arm0 Pl1");
TH1F *hStripACh3Arm0Pl2=(TH1F *)f.Get("StripA Ch3 Arm0 Pl2");
TH1F *hStripACh3Arm0Pl3=(TH1F *)f.Get("StripA Ch3 Arm0 Pl3");
TH1F *hStripACh3Arm0Pl4=(TH1F *)f.Get("StripA Ch3 Arm0 Pl4");
TH1F *hStripACh3Arm1Pl0=(TH1F *)f.Get("StripA Ch3 Arm1 Pl0");
TH1F *hStripACh3Arm1Pl1=(TH1F *)f.Get("StripA Ch3 Arm1 Pl1");
TH1F *hStripACh3Arm1Pl2=(TH1F *)f.Get("StripA Ch3 Arm1 Pl2");
TH1F *hStripACh3Arm1Pl3=(TH1F *)f.Get("StripA Ch3 Arm1 Pl3");
TH1F *hStripACh3Arm1Pl4=(TH1F *)f.Get("StripA Ch3 Arm1 Pl4");

TH1F *hStripBCh3Arm0Pl0=(TH1F *)f.Get("StripB Ch3 Arm0 Pl0");
TH1F *hStripBCh3Arm0Pl1=(TH1F *)f.Get("StripB Ch3 Arm0 Pl1");
TH1F *hStripBCh3Arm0Pl2=(TH1F *)f.Get("StripB Ch3 Arm0 Pl2");
TH1F *hStripBCh3Arm0Pl3=(TH1F *)f.Get("StripB Ch3 Arm0 Pl3");
TH1F *hStripBCh3Arm0Pl4=(TH1F *)f.Get("StripB Ch3 Arm0 Pl4");
TH1F *hStripBCh3Arm1Pl0=(TH1F *)f.Get("StripB Ch3 Arm1 Pl0");
TH1F *hStripBCh3Arm1Pl1=(TH1F *)f.Get("StripB Ch3 Arm1 Pl1");
TH1F *hStripBCh3Arm1Pl2=(TH1F *)f.Get("StripB Ch3 Arm1 Pl2");
TH1F *hStripBCh3Arm1Pl3=(TH1F *)f.Get("StripB Ch3 Arm1 Pl3");
TH1F *hStripBCh3Arm1Pl4=(TH1F *)f.Get("StripB Ch3 Arm1 Pl4");

TH1F *hWireCh4Arm0Pl0=(TH1F *)f.Get("Wire Ch4 Arm0 Pl0");
TH1F *hWireCh4Arm0Pl1=(TH1F *)f.Get("Wire Ch4 Arm0 Pl1");
TH1F *hWireCh4Arm0Pl2=(TH1F *)f.Get("Wire Ch4 Arm0 Pl2");
TH1F *hWireCh4Arm0Pl3=(TH1F *)f.Get("Wire Ch4 Arm0 Pl3");
TH1F *hWireCh4Arm0Pl4=(TH1F *)f.Get("Wire Ch4 Arm0 Pl4");
TH1F *hWireCh4Arm1Pl0=(TH1F *)f.Get("Wire Ch4 Arm1 Pl0");
TH1F *hWireCh4Arm1Pl1=(TH1F *)f.Get("Wire Ch4 Arm1 Pl1");
TH1F *hWireCh4Arm1Pl2=(TH1F *)f.Get("Wire Ch4 Arm1 Pl2");
TH1F *hWireCh4Arm1Pl3=(TH1F *)f.Get("Wire Ch4 Arm1 Pl3");
TH1F *hWireCh4Arm1Pl4=(TH1F *)f.Get("Wire Ch4 Arm1 Pl4");

TH1F *hStripACh4Arm0Pl0=(TH1F *)f.Get("StripA Ch4 Arm0 Pl0");
TH1F *hStripACh4Arm0Pl1=(TH1F *)f.Get("StripA Ch4 Arm0 Pl1");
TH1F *hStripACh4Arm0Pl2=(TH1F *)f.Get("StripA Ch4 Arm0 Pl2");
TH1F *hStripACh4Arm0Pl3=(TH1F *)f.Get("StripA Ch4 Arm0 Pl3");
TH1F *hStripACh4Arm0Pl4=(TH1F *)f.Get("StripA Ch4 Arm0 Pl4");
TH1F *hStripACh4Arm1Pl0=(TH1F *)f.Get("StripA Ch4 Arm1 Pl0");
TH1F *hStripACh4Arm1Pl1=(TH1F *)f.Get("StripA Ch4 Arm1 Pl1");
TH1F *hStripACh4Arm1Pl2=(TH1F *)f.Get("StripA Ch4 Arm1 Pl2");
TH1F *hStripACh4Arm1Pl3=(TH1F *)f.Get("StripA Ch4 Arm1 Pl3");
TH1F *hStripACh4Arm1Pl4=(TH1F *)f.Get("StripA Ch4 Arm1 Pl4");

TH1F *hStripBCh4Arm0Pl0=(TH1F *)f.Get("StripB Ch4 Arm0 Pl0");
TH1F *hStripBCh4Arm0Pl1=(TH1F *)f.Get("StripB Ch4 Arm0 Pl1");
TH1F *hStripBCh4Arm0Pl2=(TH1F *)f.Get("StripB Ch4 Arm0 Pl2");
TH1F *hStripBCh4Arm0Pl3=(TH1F *)f.Get("StripB Ch4 Arm0 Pl3");
TH1F *hStripBCh4Arm0Pl4=(TH1F *)f.Get("StripB Ch4 Arm0 Pl4");
TH1F *hStripBCh4Arm1Pl0=(TH1F *)f.Get("StripB Ch4 Arm1 Pl0");
TH1F *hStripBCh4Arm1Pl1=(TH1F *)f.Get("StripB Ch4 Arm1 Pl1");
TH1F *hStripBCh4Arm1Pl2=(TH1F *)f.Get("StripB Ch4 Arm1 Pl2");
TH1F *hStripBCh4Arm1Pl3=(TH1F *)f.Get("StripB Ch4 Arm1 Pl3");
TH1F *hStripBCh4Arm1Pl4=(TH1F *)f.Get("StripB Ch4 Arm1 Pl4");

TH1F *hWireCh5Arm0Pl0=(TH1F *)f.Get("Wire Ch5 Arm0 Pl0");
TH1F *hWireCh5Arm0Pl1=(TH1F *)f.Get("Wire Ch5 Arm0 Pl1");
TH1F *hWireCh5Arm0Pl2=(TH1F *)f.Get("Wire Ch5 Arm0 Pl2");
TH1F *hWireCh5Arm0Pl3=(TH1F *)f.Get("Wire Ch5 Arm0 Pl3");
TH1F *hWireCh5Arm0Pl4=(TH1F *)f.Get("Wire Ch5 Arm0 Pl4");
TH1F *hWireCh5Arm1Pl0=(TH1F *)f.Get("Wire Ch5 Arm1 Pl0");
TH1F *hWireCh5Arm1Pl1=(TH1F *)f.Get("Wire Ch5 Arm1 Pl1");
TH1F *hWireCh5Arm1Pl2=(TH1F *)f.Get("Wire Ch5 Arm1 Pl2");
TH1F *hWireCh5Arm1Pl3=(TH1F *)f.Get("Wire Ch5 Arm1 Pl3");
TH1F *hWireCh5Arm1Pl4=(TH1F *)f.Get("Wire Ch5 Arm1 Pl4");

TH1F *hStripACh5Arm0Pl0=(TH1F *)f.Get("StripA Ch5 Arm0 Pl0");
TH1F *hStripACh5Arm0Pl1=(TH1F *)f.Get("StripA Ch5 Arm0 Pl1");
TH1F *hStripACh5Arm0Pl2=(TH1F *)f.Get("StripA Ch5 Arm0 Pl2");
TH1F *hStripACh5Arm0Pl3=(TH1F *)f.Get("StripA Ch5 Arm0 Pl3");
TH1F *hStripACh5Arm0Pl4=(TH1F *)f.Get("StripA Ch5 Arm0 Pl4");
TH1F *hStripACh5Arm1Pl0=(TH1F *)f.Get("StripA Ch5 Arm1 Pl0");
TH1F *hStripACh5Arm1Pl1=(TH1F *)f.Get("StripA Ch5 Arm1 Pl1");
TH1F *hStripACh5Arm1Pl2=(TH1F *)f.Get("StripA Ch5 Arm1 Pl2");
TH1F *hStripACh5Arm1Pl3=(TH1F *)f.Get("StripA Ch5 Arm1 Pl3");
TH1F *hStripACh5Arm1Pl4=(TH1F *)f.Get("StripA Ch5 Arm1 Pl4");

TH1F *hStripBCh5Arm0Pl0=(TH1F *)f.Get("StripB Ch5 Arm0 Pl0");
TH1F *hStripBCh5Arm0Pl1=(TH1F *)f.Get("StripB Ch5 Arm0 Pl1");
TH1F *hStripBCh5Arm0Pl2=(TH1F *)f.Get("StripB Ch5 Arm0 Pl2");
TH1F *hStripBCh5Arm0Pl3=(TH1F *)f.Get("StripB Ch5 Arm0 Pl3");
TH1F *hStripBCh5Arm0Pl4=(TH1F *)f.Get("StripB Ch5 Arm0 Pl4");
TH1F *hStripBCh5Arm1Pl0=(TH1F *)f.Get("StripB Ch5 Arm1 Pl0");
TH1F *hStripBCh5Arm1Pl1=(TH1F *)f.Get("StripB Ch5 Arm1 Pl1");
TH1F *hStripBCh5Arm1Pl2=(TH1F *)f.Get("StripB Ch5 Arm1 Pl2");
TH1F *hStripBCh5Arm1Pl3=(TH1F *)f.Get("StripB Ch5 Arm1 Pl3");
TH1F *hStripBCh5Arm1Pl4=(TH1F *)f.Get("StripB Ch5 Arm1 Pl4");

TH1F *hRecoXArm1Pl4=(TH1F *)f.Get("Reco X Arm 1 Pl 4");
TH1F *hRecoXArm1Pl3=(TH1F *)f.Get("Reco X Arm 1 Pl 3");
TH1F *hRecoXArm1Pl2=(TH1F *)f.Get("Reco X Arm 1 Pl 2");
TH1F *hRecoXArm1Pl1=(TH1F *)f.Get("Reco X Arm 1 Pl 1");
TH1F *hRecoXArm1Pl0=(TH1F *)f.Get("Reco X Arm 1 Pl 0");
TH1F *hRecoXarm0Pl4=(TH1F *)f.Get("Reco X arm 0 Pl 4");
TH1F *hRecoXarm0Pl3=(TH1F *)f.Get("Reco X arm 0 Pl 3");
TH1F *hRecoXarm0Pl2=(TH1F *)f.Get("Reco X arm 0 Pl 2");
TH1F *hRecoXarm0Pl1=(TH1F *)f.Get("Reco X arm 0 Pl 1");
TH1F *hRecoXarm0Pl0=(TH1F *)f.Get("Reco X arm 0 Pl 0");
TH1F *hRecoYArm1Pl4=(TH1F *)f.Get("Reco Y Arm 1 Pl 4");
TH1F *hRecoYArm1Pl3=(TH1F *)f.Get("Reco Y Arm 1 Pl 3");
TH1F *hRecoYArm1Pl2=(TH1F *)f.Get("Reco Y Arm 1 Pl 2");
TH1F *hRecoYArm1Pl1=(TH1F *)f.Get("Reco Y Arm 1 Pl 1");
TH1F *hRecoYArm1Pl0=(TH1F *)f.Get("Reco Y Arm 1 Pl 0");
TH1F *hRecoYarm0Pl4=(TH1F *)f.Get("Reco Y arm 0 Pl 4");
TH1F *hRecoYarm0Pl3=(TH1F *)f.Get("Reco Y arm 0 Pl 3");
TH1F *hRecoYarm0Pl2=(TH1F *)f.Get("Reco Y arm 0 Pl 2");
TH1F *hRecoYarm0Pl1=(TH1F *)f.Get("Reco Y arm 0 Pl 1");
TH1F *hRecoYarm0Pl0=(TH1F *)f.Get("Reco Y arm 0 Pl 0");

TH1F *hTrackEtaRec=(TH1F *)f.Get("TrackEtaRec");
TH1F *hTrackPhiRec=(TH1F *)f.Get("TrackPhiRec");
TH1F *hChiSquaredOverN=(TH1F *)f.Get("ChiSquaredOverN");
TH1F *hChiSquaredProb=(TH1F *)f.Get("ChiSquaredProb");
TH1F *hChiSquaredXOverN=(TH1F *)f.Get("ChiSquaredXOverN");
TH1F *hChiSquaredXProb=(TH1F *)f.Get("ChiSquaredXProb");
TH1F *hChiSquaredYOverN=(TH1F *)f.Get("ChiSquaredYOverN");
TH1F *hChiSquaredYProb=(TH1F *)f.Get("ChiSquaredYProb");

TCanvas *cSimu = new TCanvas("Simulation","Simulation",1000,800);
cSimu->Divide(4,2);
Int_t i = 1;
cSimu->cd(i); 
i++;
hSimEta ->Draw(); // (TH1F*)f.Get("SimEta");
cSimu->cd(i);
i++;
hSimPT ->Draw(); // (TH1F)f.Get("SimPT");
cSimu->cd(i);
i++;
hSimMomentum->Draw(); //()f.Get("SimMomentum");
cSimu->cd(i);
i++;
hSimMomentumLog->Draw(); //()f.Get("SimMomentumLog");
cSimu->cd(i);
i++;
hSimEnergyloss->Draw(); //()f.Get("SimEnergyloss");
cSimu->cd(i);
i++;
hSimType->Draw(); //()f.Get("SimType");
cSimu->cd(i);
i++;
hSimUnitID->Draw(); //()f.Get("SimUnitID");
cSimu->cd();

i=1;
TCanvas *cClus = new TCanvas("Clusterization","Clusterization",1000,800);
cClus->cd(1);
hClusterWidth->Draw();
cClus->cd();

i=1;
TCanvas *cDigi00 = new TCanvas("Digitization Pl0 Arm0","Digitization Pl0 Arm0",1000,800);
cDigi00->Divide(3,6);

cDigi00->cd(i);
hWireCh0Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripACh0Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripBCh0Arm0Pl0->Draw();
i++;
cDigi00->cd(i);


cDigi00->cd(i);
hWireCh1Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripACh1Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripBCh1Arm0Pl0->Draw();
i++;
cDigi00->cd(i);



cDigi00->cd(i);
hWireCh2Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripACh2Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripBCh2Arm0Pl0->Draw();
i++;
cDigi00->cd(i);



cDigi00->cd(i);
hWireCh3Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripACh3Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripBCh3Arm0Pl0->Draw();
i++;
cDigi00->cd(i);



cDigi00->cd(i);
hWireCh4Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripACh4Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripBCh4Arm0Pl0->Draw();
i++;
cDigi00->cd(i);



cDigi00->cd(i);
hWireCh5Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripACh5Arm0Pl0->Draw();
i++;
cDigi00->cd(i);
hStripBCh5Arm0Pl0->Draw();
i++;
cDigi00->cd();
i=1;

TCanvas *cDigi01 = new TCanvas("Digitization Pl0 Arm1","Digitization Pl0 Arm1",1000,800);
cDigi01->Divide(3,6);

cDigi01->cd(i);
hWireCh0Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripACh0Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripBCh0Arm1Pl0->Draw();
i++;
cDigi01->cd(i);


cDigi01->cd(i);
hWireCh1Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripACh1Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripBCh1Arm1Pl0->Draw();
i++;
cDigi01->cd(i);



cDigi01->cd(i);
hWireCh2Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripACh2Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripBCh2Arm1Pl0->Draw();
i++;
cDigi01->cd(i);



cDigi01->cd(i);
hWireCh3Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripACh3Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripBCh3Arm1Pl0->Draw();
i++;
cDigi01->cd(i);



cDigi01->cd(i);
hWireCh4Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripACh4Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripBCh4Arm1Pl0->Draw();
i++;
cDigi01->cd(i);



cDigi01->cd(i);
hWireCh5Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripACh5Arm1Pl0->Draw();
i++;
cDigi01->cd(i);
hStripBCh5Arm1Pl0->Draw();
i++;
cDigi01->cd();
i=1;


















TCanvas *cDigi10 = new TCanvas("Digitization Pl1 Arm0","Digitization Pl1 Arm0",1000,800);
cDigi10->Divide(3,6);

cDigi10->cd(i);
hWireCh0Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripACh0Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripBCh0Arm0Pl1->Draw();
i++;
cDigi10->cd(i);


cDigi10->cd(i);
hWireCh1Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripACh1Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripBCh1Arm0Pl1->Draw();
i++;
cDigi10->cd(i);



cDigi10->cd(i);
hWireCh2Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripACh2Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripBCh2Arm0Pl1->Draw();
i++;
cDigi10->cd(i);



cDigi10->cd(i);
hWireCh3Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripACh3Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripBCh3Arm0Pl1->Draw();
i++;
cDigi10->cd(i);



cDigi10->cd(i);
hWireCh4Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripACh4Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripBCh4Arm0Pl1->Draw();
i++;
cDigi10->cd(i);



cDigi10->cd(i);
hWireCh5Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripACh5Arm0Pl1->Draw();
i++;
cDigi10->cd(i);
hStripBCh5Arm0Pl1->Draw();
i++;
cDigi10->cd();
i=1;

TCanvas *cDigi11 = new TCanvas("Digitization Pl1 Arm1","Digitization Pl1 Arm1",1000,800);
cDigi11->Divide(3,6);

cDigi11->cd(i);
hWireCh0Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripACh0Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripBCh0Arm1Pl1->Draw();
i++;
cDigi11->cd(i);


cDigi11->cd(i);
hWireCh1Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripACh1Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripBCh1Arm1Pl1->Draw();
i++;
cDigi11->cd(i);



cDigi11->cd(i);
hWireCh2Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripACh2Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripBCh2Arm1Pl1->Draw();
i++;
cDigi11->cd(i);



cDigi11->cd(i);
hWireCh3Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripACh3Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripBCh3Arm1Pl1->Draw();
i++;
cDigi11->cd(i);



cDigi11->cd(i);
hWireCh4Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripACh4Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripBCh4Arm1Pl1->Draw();
i++;
cDigi11->cd(i);



cDigi11->cd(i);
hWireCh5Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripACh5Arm1Pl1->Draw();
i++;
cDigi11->cd(i);
hStripBCh5Arm1Pl1->Draw();
i++;
cDigi11->cd();
i=1;















TCanvas *cDigi20 = new TCanvas("Digitization Pl2 Arm0","Digitization Pl2 Arm0",1000,800);
cDigi20->Divide(3,6);

cDigi20->cd(i);
hWireCh0Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripACh0Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripBCh0Arm0Pl2->Draw();
i++;
cDigi20->cd(i);


cDigi20->cd(i);
hWireCh1Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripACh1Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripBCh1Arm0Pl2->Draw();
i++;
cDigi20->cd(i);



cDigi20->cd(i);
hWireCh2Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripACh2Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripBCh2Arm0Pl2->Draw();
i++;
cDigi20->cd(i);



cDigi20->cd(i);
hWireCh3Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripACh3Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripBCh3Arm0Pl2->Draw();
i++;
cDigi20->cd(i);



cDigi20->cd(i);
hWireCh4Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripACh4Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripBCh4Arm0Pl2->Draw();
i++;
cDigi20->cd(i);



cDigi20->cd(i);
hWireCh5Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripACh5Arm0Pl2->Draw();
i++;
cDigi20->cd(i);
hStripBCh5Arm0Pl2->Draw();
i++;
cDigi20->cd();
i=1;

TCanvas *cDigi21 = new TCanvas("Digitization Pl2 Arm1","Digitization Pl2 Arm1",1000,800);
cDigi21->Divide(3,6);

cDigi21->cd(i);
hWireCh0Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripACh0Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripBCh0Arm1Pl2->Draw();
i++;
cDigi21->cd(i);


cDigi21->cd(i);
hWireCh1Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripACh1Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripBCh1Arm1Pl2->Draw();
i++;
cDigi21->cd(i);



cDigi21->cd(i);
hWireCh2Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripACh2Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripBCh2Arm1Pl2->Draw();
i++;
cDigi21->cd(i);



cDigi21->cd(i);
hWireCh3Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripACh3Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripBCh3Arm1Pl2->Draw();
i++;
cDigi21->cd(i);



cDigi21->cd(i);
hWireCh4Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripACh4Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripBCh4Arm1Pl2->Draw();
i++;
cDigi21->cd(i);



cDigi21->cd(i);
hWireCh5Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripACh5Arm1Pl2->Draw();
i++;
cDigi21->cd(i);
hStripBCh5Arm1Pl2->Draw();
i++;
cDigi21->cd();
i=1;











TCanvas *cDigi30 = new TCanvas("Digitization Pl3 Arm0","Digitization Pl3 Arm0",1000,800);
cDigi30->Divide(3,6);

cDigi30->cd(i);
hWireCh0Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripACh0Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripBCh0Arm0Pl3->Draw();
i++;
cDigi30->cd(i);


cDigi30->cd(i);
hWireCh1Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripACh1Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripBCh1Arm0Pl3->Draw();
i++;
cDigi30->cd(i);



cDigi30->cd(i);
hWireCh2Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripACh2Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripBCh2Arm0Pl3->Draw();
i++;
cDigi30->cd(i);



cDigi30->cd(i);
hWireCh3Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripACh3Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripBCh3Arm0Pl3->Draw();
i++;
cDigi30->cd(i);



cDigi30->cd(i);
hWireCh4Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripACh4Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripBCh4Arm0Pl3->Draw();
i++;
cDigi30->cd(i);



cDigi30->cd(i);
hWireCh5Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripACh5Arm0Pl3->Draw();
i++;
cDigi30->cd(i);
hStripBCh5Arm0Pl3->Draw();
i++;
cDigi30->cd();
i=1;

TCanvas *cDigi31 = new TCanvas("Digitization Pl3 Arm1","Digitization Pl3 Arm1",1000,800);
cDigi31->Divide(3,6);

cDigi31->cd(i);
hWireCh0Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripACh0Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripBCh0Arm1Pl3->Draw();
i++;
cDigi31->cd(i);


cDigi31->cd(i);
hWireCh1Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripACh1Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripBCh1Arm1Pl3->Draw();
i++;
cDigi31->cd(i);



cDigi31->cd(i);
hWireCh2Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripACh2Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripBCh2Arm1Pl3->Draw();
i++;
cDigi31->cd(i);



cDigi31->cd(i);
hWireCh3Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripACh3Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripBCh3Arm1Pl3->Draw();
i++;
cDigi31->cd(i);



cDigi31->cd(i);
hWireCh4Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripACh4Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripBCh4Arm1Pl3->Draw();
i++;
cDigi31->cd(i);



cDigi31->cd(i);
hWireCh5Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripACh5Arm1Pl3->Draw();
i++;
cDigi31->cd(i);
hStripBCh5Arm1Pl3->Draw();
i++;
cDigi31->cd();
i=1;










TCanvas *cDigi40 = new TCanvas("Digitization Pl4 Arm0","Digitization Pl4 Arm0",1000,800);
cDigi40->Divide(3,6);

cDigi40->cd(i);
hWireCh0Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripACh0Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripBCh0Arm0Pl4->Draw();
i++;
cDigi40->cd(i);


cDigi40->cd(i);
hWireCh1Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripACh1Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripBCh1Arm0Pl4->Draw();
i++;
cDigi40->cd(i);



cDigi40->cd(i);
hWireCh2Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripACh2Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripBCh2Arm0Pl4->Draw();
i++;
cDigi40->cd(i);



cDigi40->cd(i);
hWireCh3Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripACh3Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripBCh3Arm0Pl4->Draw();
i++;
cDigi40->cd(i);



cDigi40->cd(i);
hWireCh4Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripACh4Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripBCh4Arm0Pl4->Draw();
i++;
cDigi40->cd(i);



cDigi40->cd(i);
hWireCh5Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripACh5Arm0Pl4->Draw();
i++;
cDigi40->cd(i);
hStripBCh5Arm0Pl4->Draw();
i++;
cDigi40->cd();
i=1;

TCanvas *cDigi41 = new TCanvas("Digitization Pl4 Arm1","Digitization Pl4 Arm1",1000,800);
cDigi41->Divide(3,6);

cDigi41->cd(i);
hWireCh0Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripACh0Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripBCh0Arm1Pl4->Draw();
i++;
cDigi41->cd(i);


cDigi41->cd(i);
hWireCh1Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripACh1Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripBCh1Arm1Pl4->Draw();
i++;
cDigi41->cd(i);



cDigi41->cd(i);
hWireCh2Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripACh2Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripBCh2Arm1Pl4->Draw();
i++;
cDigi41->cd(i);



cDigi41->cd(i);
hWireCh3Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripACh3Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripBCh3Arm1Pl4->Draw();
i++;
cDigi41->cd(i);



cDigi41->cd(i);
hWireCh4Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripACh4Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripBCh4Arm1Pl4->Draw();
i++;
cDigi41->cd(i);



cDigi41->cd(i);
hWireCh5Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripACh5Arm1Pl4->Draw();
i++;
cDigi41->cd(i);
hStripBCh5Arm1Pl4->Draw();
i++;
cDigi41->cd();
i=1;




TCanvas *cReco = new TCanvas("Hit Reconstruction","Hit Reconstruction",1000,800);
cReco->Divide(5,4);

cReco->cd(i);

hRecoXArm1Pl4->Draw(); 
i++; cReco->cd(i); //f.Get("Reco X Arm1 Pl4");
hRecoXArm1Pl3->Draw(); i++; cReco->cd(i); //f.Get("Reco X Arm1 Pl3");
hRecoXArm1Pl2->Draw(); i++; cReco->cd(i); //f.Get("Reco X Arm1 Pl2");
hRecoXArm1Pl1->Draw(); i++; cReco->cd(i); //f.Get("Reco X Arm1 Pl1");
hRecoXArm1Pl0->Draw(); i++; cReco->cd(i); //f.Get("Reco X Arm1 Pl0");
hRecoXarm0Pl4->Draw(); i++; cReco->cd(i); //f.Get("Reco X arm0 Pl4");
hRecoXarm0Pl3->Draw(); i++; cReco->cd(i); //f.Get("Reco X arm0 Pl3");
hRecoXarm0Pl2->Draw(); i++; cReco->cd(i); //f.Get("Reco X arm0 Pl2");
hRecoXarm0Pl1->Draw(); i++; cReco->cd(i); //f.Get("Reco X arm0 Pl1");
hRecoXarm0Pl0->Draw(); i++; cReco->cd(i); //f.Get("Reco X arm0 Pl0");
hRecoYArm1Pl4->Draw(); i++; cReco->cd(i); //f.Get("Reco Y Arm1 Pl4");
hRecoYArm1Pl3->Draw(); i++; cReco->cd(i); //f.Get("Reco Y Arm1 Pl3");
hRecoYArm1Pl2->Draw(); i++; cReco->cd(i); //f.Get("Reco Y Arm1 Pl2");
hRecoYArm1Pl1->Draw(); i++; cReco->cd(i); //f.Get("Reco Y Arm1 Pl1");
hRecoYArm1Pl0->Draw(); i++; cReco->cd(i); //f.Get("Reco Y Arm1 Pl0");
hRecoYarm0Pl4->Draw(); i++; cReco->cd(i); //f.Get("Reco Y arm0 Pl4");
hRecoYarm0Pl3->Draw(); i++; cReco->cd(i); //f.Get("Reco Y arm0 Pl3");
hRecoYarm0Pl2->Draw(); i++; cReco->cd(i); //f.Get("Reco Y arm0 Pl2");
hRecoYarm0Pl1->Draw(); i++; cReco->cd(i); //f.Get("Reco Y arm0 Pl1");
hRecoYarm0Pl0->Draw(); i++; cReco->cd(i); //f.Get("Reco Y arm0 Pl0");

i=1;

cReco->cd();




TCanvas *cTrack = new TCanvas("Tracks","Tracks",1000,800);
cTrack->Divide(4,2);
cTrack->cd(i);


hTrackEtaRec->Draw(); i++; cTrack->cd(i); //f.Get("TrackEtaRec");
hTrackPhiRec->Draw(); i++; cTrack->cd(i); //f.Get("TrackPhiRec");
hChiSquaredOverN->Draw(); i++; cTrack->cd(i); //f.Get("ChiSquaredOverN");
hChiSquaredProb->Draw(); i++; cTrack->cd(i); //f.Get("ChiSquaredProb");
hChiSquaredXOverN->Draw(); i++; cTrack->cd(i); //f.Get("ChiSquaredXOverN");
hChiSquaredXProb->Draw(); i++; cTrack->cd(i); //f.Get("ChiSquaredXProb");
hChiSquaredYOverN->Draw(); i++; cTrack->cd(i); //f.Get("ChiSquaredYOverN");
hChiSquaredYProb->Draw(); i++; cTrack->cd(i); //f.Get("ChiSquaredYProb");

cTrack->cd();

}
