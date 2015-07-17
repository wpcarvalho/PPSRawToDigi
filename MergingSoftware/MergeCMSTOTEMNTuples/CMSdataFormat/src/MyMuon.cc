#include "MyMuon.h"
#include <iostream>

using namespace std;

ClassImp(MyMuon)

MyMuon::MyMuon():MyPart(){
  this->Reset();
}

MyMuon::~MyMuon(){}


void MyMuon::Reset(){

  this->MyPart::Reset();
  globalTrack.Reset();
  innerTrack.Reset();
  outerTrack.Reset();

  isoR03sumPt    = 0 ;
  isoR03emEt     = 0 ;
  isoR03hadEt    = 0 ;
  isoR03hoEt     = 0 ;
  isoR03nTracks  = 0 ;
  isoR03nJets    = 0 ;
  
  isoR05sumPt    = 0 ;
  isoR05emEt     = 0 ;
  isoR05hadEt    = 0 ;
  isoR05hoEt     = 0 ;
  isoR05nTracks  = 0 ;
  isoR05nJets    = 0 ;
  
  calEnergyEm    = 0 ;
  calEnergyHad   = 0 ;
  calEnergyHo    = 0 ;
  calEnergyEmS9  = 0 ;
  calEnergyHadS9 = 0 ;
  calEnergyHoS9  = 0 ;
  
  IsGlobalMuon	   = 0 ;
  IsTrackerMuon    = 0 ;
  IsStandaloneMuon = 0 ;
  IsCaloMuon	   = 0 ;
  
  AllGlobalMuons 		         = 0 ;
  AllStandAloneMuons		         = 0 ;
  AllTrackerMuons		         = 0 ;
  TrackerMuonArbitrated  	         = 0 ;
  AllArbitrated  		         = 0 ;
  GlobalMuonPromptTight  	         = 0 ;
  TMLastStationLoose		         = 0 ;
  TMLastStationTight		         = 0 ;
  TM2DCompatibilityLoose 	         = 0 ;
  TM2DCompatibilityTight 	         = 0 ;
  TMOneStationLoose		         = 0 ;
  TMOneStationTight		         = 0 ;
  TMLastStationOptimizedLowPtLoose       = 0 ;
  TMLastStationOptimizedLowPtTight       = 0 ;
  GMTkChiCompatibility		         = 0 ;
  GMStaChiCompatibility  	         = 0 ;
  GMTkKinkTight  		         = 0 ;
  TMLastStationAngLoose  	         = 0 ;
  TMLastStationAngTight  	         = 0 ;
  TMOneStationAngLoose		         = 0 ;
  TMOneStationAngTight		         = 0 ;
  TMLastStationOptimizedBarrelLowPtLoose = 0 ;
  TMLastStationOptimizedBarrelLowPtTight = 0 ;
 
}

void MyMuon::Print(){

  this->MyPart::Print();
  globalTrack.Print();
  innerTrack.Print();
  outerTrack.Print();

  cout << "isoR03sumPt    : " << isoR03sumPt    << endl ;
  cout << "isoR03emEt	  : " << isoR03emEt     << endl ;
  cout << "isoR03hadEt    : " << isoR03hadEt    << endl ;
  cout << "isoR03hoEt	  : " << isoR03hoEt     << endl ;
  cout << "isoR03nTracks  : " << isoR03nTracks  << endl ;
  cout << "isoR03nJets    : " << isoR03nJets    << endl ;
  
  cout << "isoR05sumPt    : " << isoR05sumPt    << endl ;
  cout << "isoR05emEt	  : " << isoR05emEt     << endl ;
  cout << "isoR05hadEt    : " << isoR05hadEt    << endl ;
  cout << "isoR05hoEt	  : " << isoR05hoEt     << endl ;
  cout << "isoR05nTracks  : " << isoR05nTracks  << endl ;
  cout << "isoR05nJets    : " << isoR05nJets    << endl ;
  
  cout << "calEnergyEm    : " << calEnergyEm    << endl ;
  cout << "calEnergyHad   : " << calEnergyHad   << endl ;
  cout << "calEnergyHo    : " << calEnergyHo    << endl ;
  cout << "calEnergyEmS9  : " << calEnergyEmS9  << endl ;
  cout << "calEnergyHadS9 : " << calEnergyHadS9 << endl ;
  cout << "calEnergyHoS9  : " << calEnergyHoS9  << endl ;
  
  cout << "IsGlobalMuon     : " << IsGlobalMuon     << endl ;
  cout << "IsTrackerMuon    : " << IsTrackerMuon    << endl ;
  cout << "IsStandaloneMuon : " << IsStandaloneMuon << endl ;
  cout << "IsCaloMuon	    : " << IsCaloMuon	    << endl ;
  
  cout << "AllGlobalMuons			  : " << AllGlobalMuons 		        << endl ;
  cout << "AllStandAloneMuons			  : " << AllStandAloneMuons		        << endl ;
  cout << "AllTrackerMuons			  : " << AllTrackerMuons		        << endl ;
  cout << "TrackerMuonArbitrated		  : " << TrackerMuonArbitrated  	        << endl ;
  cout << "AllArbitrated			  : " << AllArbitrated  		        << endl ;
  cout << "GlobalMuonPromptTight		  : " << GlobalMuonPromptTight  	        << endl ;
  cout << "TMLastStationLoose			  : " << TMLastStationLoose		        << endl ;
  cout << "TMLastStationTight			  : " << TMLastStationTight		        << endl ;
  cout << "TM2DCompatibilityLoose		  : " << TM2DCompatibilityLoose 	        << endl ;
  cout << "TM2DCompatibilityTight		  : " << TM2DCompatibilityTight 	        << endl ;
  cout << "TMOneStationLoose			  : " << TMOneStationLoose		        << endl ;
  cout << "TMOneStationTight			  : " << TMOneStationTight		        << endl ;
  cout << "TMLastStationOptimizedLowPtLoose	  : " << TMLastStationOptimizedLowPtLoose       << endl ;
  cout << "TMLastStationOptimizedLowPtTight	  : " << TMLastStationOptimizedLowPtTight       << endl ;
  cout << "GMTkChiCompatibility 		  : " << GMTkChiCompatibility		        << endl ;
  cout << "GMStaChiCompatibility		  : " << GMStaChiCompatibility  	        << endl ;
  cout << "GMTkKinkTight			  : " << GMTkKinkTight  		        << endl ;
  cout << "TMLastStationAngLoose		  : " << TMLastStationAngLoose  	        << endl ;
  cout << "TMLastStationAngTight		  : " << TMLastStationAngTight  	        << endl ;
  cout << "TMOneStationAngLoose 		  : " << TMOneStationAngLoose		        << endl ;
  cout << "TMOneStationAngTight 		  : " << TMOneStationAngTight		        << endl ;
  cout << "TMLastStationOptimizedBarrelLowPtLoose : " << TMLastStationOptimizedBarrelLowPtLoose << endl ;
  cout << "TMLastStationOptimizedBarrelLowPtTight : " << TMLastStationOptimizedBarrelLowPtTight << endl ;

}
