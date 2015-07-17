

// system include files
#include <memory>
//include <cmath>

#include "CLHEP/Vector/ThreeVector.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDigCluster.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
//include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include <DataFormats/T2DetId/interface/T2DetId.h>
#include <DataFormats/T1Road/interface/T1RecHitGlobal.h>
#include <DataFormats/T2Hit/interface/T2Hit.h>

//Data Formats
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
//include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"

#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
//include "TotemAnalysis/TotemNtuplizer/interface/TriggerDataFormats.h"


#include "RecoTotemRP/RPRecoDataFormats/interface/RPMulFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPMulTrackCandidateCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonCollection.h"

#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/TotemRPDetTopology/interface/RPTopology.h"

#include "DataFormats/TotemL1Trigger/interface/RPCCBits.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/RPRecoProtMADXVariables.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"


#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <map>

//define FOdebug 1

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
//include "TH1I.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TFile.h"
#include "TVector2.h"
#include "TDirectory.h"

using namespace std;

typedef std::vector<RPCCBits> CCVect;
typedef std::map<int, int> MapType;
typedef std::multimap<double, double> DMapType;
//typedef std::auto_ptr<TH1D> point1DH;


/// a copy of Totem::TriggerData from TotemRawData/DataFormats/interface/RawEvent.h
struct TriggerData
{
  unsigned char type;
  unsigned int event_num, bunch_num, src_id;
  unsigned int orbit_num;
  unsigned char revision_num;
  unsigned int run_num, trigger_num, inhibited_triggers_num, input_status_bits;
  
  TriggerData() {}
  virtual ~TriggerData() {}
};


//
// class decleration
//

class RPDataReduction : public edm::EDAnalyzer {
   public:
      explicit RPDataReduction(const edm::ParameterSet&);
      ~RPDataReduction();


   private:
      edm::InputTag rPReconstructedProtonCollectionLabel;
      edm::InputTag rawEventLabel;
      edm::InputTag detSetVectorRPDigClusterLabel;
      edm::InputTag rPMulFittedTrackCollectionLabel;
      edm::InputTag rPFittedTrackCollectionLabel;
      edm::InputTag rPRecognizedPatternsCollectionLabel;
      edm::InputTag rPTrackCandidateCollectionLabel;
      edm::InputTag rPMulTrackCandidateCollectionLabel;

      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void beginRun(edm::Run const&, edm::EventSetup const& es);
      virtual void endJob() ;

      // ----------member data ---------------------------
  int totalClu[2][4][8];
  int cluSize[2][4][8];
  int partialEvts[2][4][8];
  int perPlaneEvts[1280];
  float SOccup[2][4][8];
  int totalEvts;
  bool readMultiTrk;
  bool readT2;
  bool readT1;
  bool readLoNeg;
  bool readClusters;
  bool readRecoProt;
  unsigned int diagonal;
  unsigned int potsD2[9][4];


  MapType potAll4On,potAll4Off,pots3On1Off,pots3On1Half,potAll4Trk,potTrk3S1M,potTrk2S2M,pots3On1HalfAll;
  MapType pot3Plus1Medium,pot3Plus1Max,pot2Plus2Medium,pot2T2Empty;
  MapType TwoOnTwoOff,TwoOnTwoHalf,TwoOnTwoVMul,TwoOnTwoHalfAll;
  MapType TwoOn2ATwoOff,TwoOn2ATwoHalf,TwoOn2ATwoVMul;
  MapType TwoOn2ATwoHalfAll,pot2A2Plus2Medium,pot2ATrk2S2M;
  MapType pot3On1HHu,pot3On1HHv,OneDiag00SomepotOn;
  MapType potTrk3S1NoClu,potTrk2S2NoClu, pot2ATrk2S2NoClu;
  // not in use yet : pot2On2HHu,pot2On2HHv,
  MapType trig220cross,trig220h,trig220v,trigT2,trigT1,trigBx;
  MapType rpSTrT2NoTr,rpT2TrSame,rpT2TrOppos,rpT2TrBoth,badRpT2TrSame,badRpT2TrOppos,badRpT2TrBoth,rpBadTrT2NoTr;
  MapType AllPot00,AllPotUorV;
  MapType rpSTrT2NoTrFC,rpT2TrSameFC,rpT2TrOpposFC,rpT2TrBothFC; // Few clusters : demand <10 clusters/pot in all 
                                                                 //6 pots in other arm from the detected proton
  MapType rpTmcT2NoTr,rpTmcT2TrSame,rpTmcT2TrOppo,rpTmcT2TrBoth; // TMC on in two neighbour pots only, but no STrk 
  MapType rpTmcT2NoTrT,rpTmcT2TrSameT,rpTmcT2TrOppoT,rpTmcT2TrBothT; // TMC on in two neighbour pots only, but no STrk .AND. T2 triggered
  MapType rpTmcT2NoTrRT,rpTmcT2TrSameRT,rpTmcT2TrOppoRT,rpTmcT2TrBothRT; // TMC on in two neighbour pots only, but no STrk .AND. T2 triggered
  MapType rpSTrT2NoTrT2T,rpT2TrSameT2T,rpT2TrOpposT2T,rpT2TrBothT2T; // event was triggered by T2, according to LoNeg
  MapType rp2T2n[6],rp2T2s[6],rp2T2o[6],rp2T2b[6];

  float refEvtsOn,refEvtsOff,refEvts;
  float checkEvtsOn,checkEvtsOff,checkEvts;
  float noChkTrkOn,noChkTrkOff,refEvtsVale;
  float MTrCand01AnyPot,TrCand01AnyPot,MTrC01_diag31,TrC01_diag31,MTrC01_diag22,TrC01_diag22;
  float AnyPot01,MTr01AnyPot,STr01AnyPot,mptyEvts;
  //  float TwoDiag00,TwoDiagUorV; //rpSTrT2NoHit not yet used
  float potDiag3STrk,  potDiag2STrk, potDiag1STrk,potDiag0STrk;
  std::auto_ptr<TProfile> prfNClu;
  std::auto_ptr<TProfile> prfSClu;
  std::auto_ptr<TProfile> prfOccup;
  std::auto_ptr<TH1D> trgU;
  std::auto_ptr<TH1D> trgV;
  std::auto_ptr<TH1D> srcX;
  std::auto_ptr<TH1D> trkTx;
  std::auto_ptr<TH1D> lonegBunch[7];
  std::auto_ptr<TH1D> sdTrigs;
  std::auto_ptr<TH1D> spectroElastRefX;
  std::auto_ptr<TH1D> trkTy;
  std::auto_ptr<TH1D> TwoRPTrkSumCluRest[4];
  std::auto_ptr<TH1D> pullX;
  std::auto_ptr<TH2D> pullXpY;
  std::auto_ptr<TH2D> refTSvsExtrap;
  std::auto_ptr<TH2D> refTSvsExtrapU;
  std::auto_ptr<TH2D> refTSvsCheckU;
  std::auto_ptr<TH2D> refTSvsCheckV;
  std::auto_ptr<TH2D> pullXXc;
  std::auto_ptr<TH2D> recoProtXiVsN;
  std::auto_ptr<TH2D> recoProtXiVsErr;
  std::auto_ptr<TH1D> TwoRPTrkNumProt;
  std::auto_ptr<TH2D> pullYVsSpec;
  std::auto_ptr<TH2D> UvsTMCu;
  std::auto_ptr<TH1D> spectro45tp;
  std::auto_ptr<TH1D> spectro45bt;
  std::auto_ptr<TH1D> spectro56tp;
  std::auto_ptr<TH1D> spectro56bt;
  std::auto_ptr<TH2D> UvsTMCv;
  std::auto_ptr<TH2D> VvsTMCu;
  std::auto_ptr<TH2D> VvsTMCv;
  std::auto_ptr<TH2D> pullYYc;
  std::auto_ptr<TH2D> spectroElastUvsV;


  std::auto_ptr<TH2D> trkMisPointXY;
 
  std::auto_ptr<TH2D> refUVsCheckU;
  std::auto_ptr<TH2D> refVVsCheckV;
  std::auto_ptr<TH2D> refEvtXY;
  std::auto_ptr<TProfile> elasticRefVsCheckY;
  std::auto_ptr<TH2D> chkEvtXY;
  std::auto_ptr<TH1D> pullY;
  std::auto_ptr<TH1D> pullYElast;
  //  std::auto_ptr<TH3D> pullYvsXY;
  std::auto_ptr<TH1D> numHitsInCheckPointing;
  std::auto_ptr<TH1D> srcY;
  std::auto_ptr<TH1D> chkY;
  std::auto_ptr<TH1D> srcYHit;
  std::auto_ptr<TH1D> chkYHit;
  std::auto_ptr<TH1D> srcU;
  std::auto_ptr<TH1D> srcV;
  std::auto_ptr<TH1D> srcUpos;
  std::auto_ptr<TH1D> srcVpos;
  std::auto_ptr<TH1D> trgUexPos;
  std::auto_ptr<TH1D> trgVexPos;
  std::auto_ptr<TH1D> trksWOneTmcBitOn;
  std::auto_ptr<TH1D> NotrksWOneTmcBitOn;
  std::auto_ptr<TProfile> prfMultipl;
  std::auto_ptr<TProfile> prfMul[12];
  std::auto_ptr<TH1D> multiPot;
  std::auto_ptr<TH1D> mul01;
  std::auto_ptr<TH1D> sing01;
  std::auto_ptr<TH1D> num01;
  std::auto_ptr<TH1D> mul01cand;
  std::auto_ptr<TH1D> sing01cand;
  std::auto_ptr<TH1D> mul01diag3p1;
  std::auto_ptr<TH1D> mul01diag2p2;
  std::auto_ptr<TH1D> oneOriOnly3p1;
  std::auto_ptr<TH1D> oneOriOnly2p2;
  std::auto_ptr<TH1D> oneOriOnly3p011S;
  std::auto_ptr<TH1D> uv11Times3Plus001S;


  std::auto_ptr<TH1D> outputPercentagesD;
  std::auto_ptr<TH1D> outputFulfillEvtsD;
  std::auto_ptr<TH1D> outputPercentagesA;
  std::auto_ptr<TH1D> outputFulfillEvtsA;
  std::auto_ptr<TH1D> outputRefEvts;
  std::auto_ptr<TH1D> outputPerPot[12];
 

  std::auto_ptr<TH1D> outputPerBunchD[10000];
  std::auto_ptr<TH1D> percentagesPerBunchD[10000];
  std::auto_ptr<TH1D> outputPerBunchA[10000];
  std::auto_ptr<TH1D> percentagesPerBunchA[10000];


  std::auto_ptr<TH2D> T2FitEtaVs000Eta;
  std::auto_ptr<TH1D> T2Eta000Dispersion;
  std::auto_ptr<TH2D> T1FitEtaVs000Eta;
  std::auto_ptr<TH1D> T1Eta000Dispersion;

  MapType refEvtsBunches;
  MapType mp3Trk,mp2Trk,mp1Trk,mp0Trk;

  //  MapType mp411,mp400,mp311t100,mp311t101;
 
  std::auto_ptr<TH1D> SDTrkWhichPot[4];
  std::auto_ptr<TH1D> SDBadTrkWhichPot[4];
  std::auto_ptr<TH2D> SDXNearVsFar[4][6];
  std::auto_ptr<TH2D> SDYNearVsFar[4][6];
  std::auto_ptr<TH1D> SDXNearMinusFar[4][6];
  std::auto_ptr<TH1D> SDYNearMinusFar[4][6];

  std::auto_ptr<TH2D> otherSideCluVsMaxTMC[4];
  std::auto_ptr<TH2D> RPTrigOtherSideCluVsMaxTMC[4];
  std::auto_ptr<TH2D> otherSideCluVsCluPerPot[4];
  std::auto_ptr<TH2D> RPTrigOtherSideCluVsCluPerPot[4];
  std::auto_ptr<TH2D> BNon0OtherSideClu[4];
  std::auto_ptr<TH2D> BNon0RPTrigOtherSideClu[4];


  std::auto_ptr<TH1D> SDtrkX[4];
  std::auto_ptr<TH1D> SDtrkY[4];
  std::auto_ptr<TH1D> SDT2Multi[4];
  std::auto_ptr<TH2D> SDXvsY[4][4];
  std::auto_ptr<TH2D> SD2pAngleXvsX[4];
  std::auto_ptr<TH2D> SD2pAngleYvsY[4];

  std::auto_ptr<TH2D> okRgAngleXvsX;
  std::auto_ptr<TH2D> badRgAngleXvsX;
  std::auto_ptr<TH2D> FiverAngleXvsX;
  std::auto_ptr<TH2D> FiverThetaYMineVsMadx;
  std::auto_ptr<TH2D> FiverThetaXMineVsMadx;
  std::auto_ptr<TH2D> FiverXvsY;
  std::auto_ptr<TH1D> T1T2numTrk[4];
  std::auto_ptr<TH2D> rgCollimatorXY[2];

  std::auto_ptr<TH2D> T2TrigSD2pAngleXvsX[4];
  std::auto_ptr<TH2D> T2noRPTrigSD2pAngleXvsX[4];
  std::auto_ptr<TH2D> T2SpectroQuarter[4][6];
  std::auto_ptr<TH2D> T2SpectroQuarterBB[4][6];
  std::auto_ptr<TH2D> T2SpectroQuarterOB[4][6];

  std::auto_ptr<TH2D> T2TrigSD2pAngleXvsXBB[4];
  std::auto_ptr<TH2D> T2TrigSD2pAngleXvsXOB[4];


  std::auto_ptr<TH1D> SDCluNoTrk[4];
  std::auto_ptr<TH1D> rpTrgSDCluNoTrk[4];

  std::auto_ptr<TH1D> rpTee[4];
  std::auto_ptr<TH1D> t2rpTee[4];

  std::auto_ptr<TH1D> RgXiOkrpTee[2][2];
 
  std::auto_ptr<TH1D> t2rpTeePair[4][6];
  std::auto_ptr<TH1D> SDt2rpTeePair[4][6];

  std::auto_ptr<TH2D> xiRpVsXiRG[4];
  std::auto_ptr<TH2D> xiRGT2oppo[3];
  std::auto_ptr<TH2D> xiRGT2both[3];
  std::auto_ptr<TH1D> RGTrkVsRGxiT1BothT2Oppo;
  std::auto_ptr<TH2D> RGTrk2DRGxiT1BothT2Oppo;
  std::auto_ptr<TH2D> RGTrk2DRGxiT1T2Both;

  std::auto_ptr<TH2D> YcutXiRGT2oppo[3][6];
  std::auto_ptr<TH2D> noYcutXiRGT2oppo[3][6];

  std::auto_ptr<TH2D> xiBadOkRGT2both[2][6];

  std::auto_ptr<TH2D> t2rpTeeVsXi[4];
  std::auto_ptr<TH2D> ElastT2rpTeeVsXi[4];
  std::auto_ptr<TH1D> ElastT2rpTee[4];
  std::auto_ptr<TH1D> SdT2rpTee[4];

  std::auto_ptr<TH2D> allTrkCollimatorXY[4];

  std::auto_ptr<TH2D> funcTeeVsElastTee[4];
  //  std::auto_ptr<TH2D> funcTeeMinusETeeVsXi[4];
  std::auto_ptr<TH2D> funcTeeRelErrVsXi[4];
  //  std::auto_ptr<TH2D> ipAngleXMadxVsMine[4];
  //  std::auto_ptr<TH2D> ipAngleYMadxVsMine[4];
  std::auto_ptr<TH1D> ipAngleXRelErr[4];
  std::auto_ptr<TH1D> ipAngleYRelErr[4];


  std::auto_ptr<TH2D> ipAnglXREvsX[4];
  std::auto_ptr<TH2D> ipAnglYREvsY[4];

  std::auto_ptr<TH1D> ElaIpAngleXRelErr[4];
  std::auto_ptr<TH1D> ElaIpAngleYRelErr[4];
  std::auto_ptr<TH1D> T2ipAngleX[4];
  std::auto_ptr<TH1D> T2ipAngleY[4];
  std::auto_ptr<TH1D> ElaT2ipAngleX[4];
  std::auto_ptr<TH1D> ElaT2ipAngleY[4];
  std::auto_ptr<TH1D> SdT2ipAngleX[4];
  std::auto_ptr<TH1D> SdT2ipAngleY[4];
  std::auto_ptr<TH1D> numRecoPInSD;

  std::auto_ptr<TH1D> rpvBunchTrigs,rphBunchTrigs,t2BunchTrigs,bxBunchTrigs,sdEvBunches;
 
  int singleTOnly[128];
  int noTrk[128];
  int noTrkEff[128];
  int evtsEff[128];
  int noTrkSomeClus[128];
  //  int totalEvts;
  int multiTrk[128];
  int falseBunch,bigBunch;
  unsigned int nPl,nMaxPrPl,refPot,checkPot,whichVerb,oppositePot,modeF,diagPot1,diagPot2;
  double sigmas,elCut56,errX,errY,elCut45t,elCut45b,exxLow,exxHi,optLy,optDLx,optEb,cutX,cutAngleX,optVx,optDVx,optLx,sigmaXi,sigmaRg;
  std::string outFile,tmcModule,tmcProd,t2trModule,t2trProd,t1trModule,t1trProd,t2padModule,t2padProd,bunchText;

  MapType rpvTrig;
  MapType rphTrig;
  MapType bxTrig;
  MapType t2Trig;
  MapType sdBunchEv;

  TriggerData dataT;

  BeamOpticsParams optObj;

  //  static

};

//
// constants, enums and typedefs
//

  enum countT2{kT2none,kT2same,kT2oppo,kT2both,kT2NA};
  enum countT1{kT1SameT2Oppo,kT1T2Oppo,kNoT1};
  enum countT1b{kT1both,kT1either,kT1neither};
  const double nksi=4 ;


// static data member definitions
//
