// -*- C++ -*-
//
// Package:    RPDataReduction
// Class:      RPDataReduction
// 
/**\class RPDataReduction RPDataReduction.cc RecoTotemRP/RPDataReduction/src/RPDataReduction.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Fredrik Oljemark
//         Created:  Wed Nov 10 09:40:03 CET 2010
// $Id$
//
//


#include "RecoTotemRP/RPDataReduction/interface/RPDataReduction.h"

//ClassImp(TriggerData)

//
// constructors and destructor
//
RPDataReduction::RPDataReduction(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  rPReconstructedProtonCollectionLabel = iConfig.getParameter<edm::InputTag>("RPReconstructedProtonCollectionLabel");
  rawEventLabel = iConfig.getParameter<edm::InputTag>("RawEventLabel");  
  detSetVectorRPDigClusterLabel = iConfig.getParameter<edm::InputTag>("DetSetVectorRPDigClusterLabel");
  rPMulFittedTrackCollectionLabel = iConfig.getParameter<edm::InputTag>("RPMulFittedTrackCollectionLabel");
  rPFittedTrackCollectionLabel = iConfig.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");  
  rPRecognizedPatternsCollectionLabel = iConfig.getParameter<edm::InputTag>("RPRecognizedPatternsCollectionLabel");
  rPTrackCandidateCollectionLabel = iConfig.getParameter<edm::InputTag>("RPTrackCandidateCollectionLabel");
  rPMulTrackCandidateCollectionLabel = iConfig.getParameter<edm::InputTag>("RPMulTrackCandidateCollectionLabel");

  nPl=iConfig.getParameter<unsigned int> ("nPl");
  nMaxPrPl=iConfig.getParameter<unsigned int> ("nMaxPrPl");
  refPot=iConfig.getParameter<unsigned int> ("refPot");
  oppositePot=iConfig.getParameter<unsigned int> ("oppositePot");
  checkPot=iConfig.getParameter<unsigned int> ("checkPot");
  modeF=iConfig.getParameter<unsigned int> ("andNotOr");
  diagonal=iConfig.getParameter<unsigned int> ("whichDiag");
  readMultiTrk=iConfig.getParameter<bool> ("readMultiTrk");
  readT2=iConfig.getParameter<bool> ("readT2");
  readT1=iConfig.getParameter<bool> ("readT1");
  readLoNeg=iConfig.getParameter<bool> ("readLoNeg");
  readClusters=iConfig.getParameter<bool> ("readClusters");
  readRecoProt=iConfig.getParameter<bool> ("readRecoProt");
  falseBunch=iConfig.getParameter<int> ("whichBunchAll");
  bigBunch=iConfig.getParameter<int> ("bigBunch");
  unsigned int tempF2[9][4]={{20,24,121,125},{21,25,120,124},{20,24,120,124},{22,23,122,123},
			     {21,25,121,125},{20,24,122,123},{21,25,122,123},{22,23,120,124},{22,23,121,125}};
  for (int a=0;a<9;a++)
    {
      for (int b=0;b<4;b++)
	{
	  potsD2[a][b]=tempF2[a][b];
	}
    }

  sigmas=iConfig.getParameter<double> ("sigmaCut");
  elCut56=iConfig.getParameter<double> ("elCut56");
  elCut45t=iConfig.getParameter<double> ("elCut45t");
  elCut45b=iConfig.getParameter<double> ("elCut45b");
  exxLow=iConfig.getUntrackedParameter<double> ("xLow");
  exxHi=iConfig.getUntrackedParameter<double> ("xHigh");

  errX=iConfig.getParameter<double> ("errX");
  errY=iConfig.getParameter<double> ("errY");
  whichVerb=iConfig.getParameter<unsigned int> ("verb");
  outFile=iConfig.getParameter<std::string> ("fileName");
  bunchText=iConfig.getParameter<std::string> ("bunchText");
  tmcModule=iConfig.getParameter<std::string> ("tmcModule");
  tmcProd=iConfig.getParameter<std::string> ("tmcProd");
  t2trModule=iConfig.getParameter<std::string> ("t2trModule");
  t2trProd=iConfig.getParameter<std::string> ("t2trProd");
  t1trModule=iConfig.getParameter<std::string> ("t1trModule");
  t1trProd=iConfig.getParameter<std::string> ("t1trProd");
  t2padModule=iConfig.getParameter<std::string> ("t2padModule");
  t2padProd=iConfig.getParameter<std::string> ("t2padProd");

  optDLx=iConfig.getParameter<double> ("dLxDs");
  optLy=iConfig.getParameter<double> ("Ly");
  optLx=iConfig.getParameter<double> ("Lx");
  optEb=iConfig.getParameter<double> ("Eb");
  cutX=iConfig.getParameter<double> ("ElaSigmaX");
  cutAngleX=iConfig.getParameter<double> ("ElaSigmaAngleX");
  optDVx=iConfig.getParameter<double> ("dVxDs");
  optVx=iConfig.getParameter<double> ("Vx");


  sigmaXi=iConfig.getParameter<double> ("sigmaXi");
  sigmaRg=iConfig.getParameter<double> ("sigmaRg");
 

 //  optObj= std::auto_ptr<BeamOpticsParams> (new BeamOpticsParams(iConfig));


  //  const 
  //enum countT2 {kT2none,kT2same,kT2oppo,kT2both,kT2NA};

  
  int delta1=0;
  int delta2=0;

  delta1=  ((refPot>100) ? -100 : 100);
  delta2=  ((refPot%2) ? -1 : 1);
  diagPot1=refPot+delta1+delta2;
  diagPot2=checkPot+delta1+delta2;

  Double_t CBins[100];
  for (int u=0;u<10;u++)
    {
      CBins[u]=u-0.5;
      CBins[10+u]=10+2.5*u;
      CBins[20+u]=35+7*u;
      CBins[30+u]=105+18*u;
      CBins[40+u]=285+45*u;
      CBins[50+u]=735+94.5*u;
      CBins[60+u]=1680+252*u;
      CBins[70+u]=4200+750*u;
      CBins[80+u]=11700+2000*u;
      CBins[90+u]=31700+2000*u;
    }
  /*
  Double_t TBins[90]={-10000,-1000,-500,-200,-100};
  for (int tee=0;tee<5;tee++)
    {
      TBins[5+tee]=-80+10*tee;
      TBins[10+tee]=-35+5*tee;
      TBins[15+tee]=-12.5+2.5*tee;
      TBins[20+tee]=-2+0.25*tee;
      TBins[25+tee]=-0.9+0.1*tee;
      TBins[30+tee]=-0.45+0.05*tee;
      TBins[35+tee]=-0.2+0.02*tee;
      TBins[40+tee]=-0.11+0.01*tee;
      TBins[45+tee]=-0.06+0.01*tee;
      TBins[50+tee]=-0.01+0.01*tee;
      TBins[55+tee]=0.04+0.1*tee;
      TBins[60+tee]=0.5+0.5*tee;
      TBins[65+tee]=3+2.5*tee;
      TBins[70+tee]=15+15*tee;
      TBins[75+tee]=100+450*tee;
      TBins[80+tee]=2200+1000*tee;
      TBins[85+tee]=7000+3000*tee;

    }
  std::cout<<"Tbins : ";
  for (int tee=0;tee<90;tee++)
    std::cout<<TBins[tee]<<"  ";
  std::cout<<std::endl;

  */

  Double_t CluBins[15]={-0.5,0.5,5.5,10.5,20.5,30.5,50.5,75.5,125.5,250.5,500.5,1000.5,1500.5,3000.5,6000.5};
  Double_t CluBins2[18]={-0.5,0.5,2.5,4.5,6.5,8.5,10.5,15.5,20.5,25.5,30.5,35.5,45.5,65.5,90.5,180.5,720.5,6000.5};


  if (whichVerb>10)
    {
      for (int yy=0;yy<100;yy++)
	{
	  std::cout<<CBins[yy]<<" ";
	  if (CBins[yy]>CBins[yy-1])
	    std::cout<<"!!";
	  if (CBins[yy]==CBins[yy-1])
	    std::cout<<"equal!!";
	}
      std::cout<<std::endl;

    }
  trgU=std::auto_ptr<TH1D>(new TH1D("trgU","trigBit U in target pot",17,-0.5,16.5));
  trgV=std::auto_ptr<TH1D>(new TH1D("trgV","trigBit V in target pot",17,-0.5,16.5));
  srcU=std::auto_ptr<TH1D>(new TH1D("srcU","trigBit U in reference pot",17,-0.5,16.5));
  srcV=std::auto_ptr<TH1D>(new TH1D("srcV","trigBit V in reference pot",17,-0.5,16.5));



  mul01=std::auto_ptr<TH1D>(new TH1D("mul01","The pot with a multitrack and (U,V)=(0,1) or (1,0)",129,-0.5,128.5));
  mul01cand=std::auto_ptr<TH1D>(new TH1D("mul01cand","The pot with a multitrack candidate and (U,V)=(0,1) or (1,0)",129,-0.5,128.5));
  sing01=std::auto_ptr<TH1D>(new TH1D("sing01","The pot with a nonParallel track and (U,V)=(0,1) or (1,0)",129,-0.5,128.5));
  sing01cand=std::auto_ptr<TH1D>(new TH1D("sing01cand","The pot with a nonPara track candidate and (U,V)=(0,1) or (1,0)",129,-0.5,128.5));
  num01=std::auto_ptr<TH1D>(new TH1D("num01","The pot with TMC (U,V)=(0,1) or (1,0)",129,-0.5,128.5));

  std::string whatCombi[9]={ "(45tp*56bt)","(45bt*56tp)", "(45tp*56tp)","(45hr*56hr)", "(45bt*56bt)", "(45tp*56hr)","(45bt*56hr)", "(45hr*56tp)","(45hr*56bt)"};

  std::string hisTitleA;
  std::string hisTitle="The pot with a multitrack and (U,V)=(0,1) : 1*(1,0)+3*(1,1) ";
  hisTitle.append(whatCombi[diagonal]);
  mul01diag3p1=std::auto_ptr<TH1D>(new TH1D("mul01diag3p1",hisTitle.c_str(),129,-0.5,128.5));
  hisTitle=  "The pots with a multitrack and (U,V)=(0,1): 2*(1,0)+2*(1,1) ";
  hisTitle.append(whatCombi[diagonal]);
  mul01diag2p2=std::auto_ptr<TH1D>(new TH1D("mul01diag2p2",hisTitle.c_str(),129,-0.5,128.5));
  hisTitle="The pot with a multitrack but no singleTrk in events with 3S+1M "; 
  hisTitle.append(whatCombi[diagonal]);
  multiPot=std::auto_ptr<TH1D>(new TH1D("multiPot",hisTitle.c_str(),129,-0.5,128.5));
  hisTitle="The pot with (U,V)=(0,1) or vice versa in evts with 3(1,1)+1(0,1) ";
  hisTitle.append(whatCombi[diagonal]);
  oneOriOnly3p1=std::auto_ptr<TH1D>(new TH1D("oneOriOnly3p1",hisTitle.c_str(),129,-0.5,128.5));
  hisTitle="The pots with (U,V)=(0,1) or vice versa in evts with 2(1,1)+2(0,1) ";
  hisTitle.append(whatCombi[diagonal]);
  oneOriOnly2p2=std::auto_ptr<TH1D>(new TH1D("oneOriOnly2p2",hisTitle.c_str(),129,-0.5,128.5));
  hisTitle="The pot with (U,V)=(0,1) or vice versa .AND. a singleTrack(nonpar.) in evts with 3(1,1)+1(0,1) ";
  hisTitle.append(whatCombi[diagonal]);
  oneOriOnly3p011S=std::auto_ptr<TH1D>(new TH1D("oneOriOnly3p011S",hisTitle.c_str(),129,-0.5,128.5));
  hisTitle="The pot with (U,V)=(0,0) .AND. a singleTrack(nonpar.) in evts with 3(1,1)+1(0,0) ";
  hisTitle.append(whatCombi[diagonal]);
  uv11Times3Plus001S=std::auto_ptr<TH1D>(new TH1D("uv11Times3Plus001S",hisTitle.c_str(),129,-0.5,128.5));

  TwoRPTrkNumProt=std::auto_ptr<TH1D>(new TH1D("TwoRPTrkNumProt","Number of reconstructed protons in events with trk in 2/12 RP (nr+fr,no other cuts)",10,-1.5,8.5));

/*
  if (diagonal)
    {
      mul01diag3p1=std::auto_ptr<TH1D>(new TH1D("mul01diag3p1","The pot with a multitrack and (U,V)=(0,1) in diag 45bt*56tp 1*(1,0)+3*(1,1)",129,-0.5,128.5));
      mul01diag2p2=std::auto_ptr<TH1D>(new TH1D("mul01diag2p2","The pots with a multitrack and (U,V)=(0,1) in diag 45bt*56tp 2*(1,0)+2*(1,1)",129,-0.5,128.5));
      multiPot=std::auto_ptr<TH1D>(new TH1D("multiPot","The pot with a multitrack but no singleTrk in events with 3S+1M (in 45bt*56tp)",129,-0.5,128.5));
      oneOriOnly3p1=std::auto_ptr<TH1D>(new TH1D("oneOriOnly3p1","The pot with (U,V)=(0,1) or vice versa in evts with 3(1,1)+1(0,1) (in 45bt*56tp)",129,-0.5,128.5));
      oneOriOnly2p2=std::auto_ptr<TH1D>(new TH1D("oneOriOnly2p2","The pots with (U,V)=(0,1) or vice versa in evts with 2(1,1)+2(0,1) (in 45bt*56tp)",129,-0.5,128.5));
    }
  else
    {
      mul01diag3p1=std::auto_ptr<TH1D>(new TH1D("mul01diag3p1","The pot with a multitrack and (U,V)=(0,1) in diag 45tp*56bt 1*(1,0)+3*(1,1)",129,-0.5,128.5));
      mul01diag2p2=std::auto_ptr<TH1D>(new TH1D("mul01diag2p2","The pots with a multitrack and (U,V)=(0,1) in diag 45tp*56bt 2*(1,0)+2*(1,1)",129,-0.5,128.5));


      multiPot=std::auto_ptr<TH1D>(new TH1D("multiPot","The pot with a multitrack but no singleTrk in events with 3S+1M (in 45tp*56bt)",129,-0.5,128.5));
      oneOriOnly3p1=std::auto_ptr<TH1D>(new TH1D("oneOriOnly3p1","The pot with (U,V)=(0,1) or vice versa in evts with 3(1,1)+1(0,1) (in 45tp*56bt)",129,-0.5,128.5));
      oneOriOnly2p2=std::auto_ptr<TH1D>(new TH1D("oneOriOnly2p2","The pots with (U,V)=(0,1) or vice versa in evts with 2(1,1)+2(0,1) (in 45tp*56bt)",129,-0.5,128.5));
    }
  */
  outputPerPot[0]=std::auto_ptr<TH1D>(new TH1D("outputPerPot020","Track & trigger : 1-pot (020) number of events",33,-0.5,32.5));
  outputPerPot[1]=std::auto_ptr<TH1D>(new TH1D("outputPerPot021","Track & trigger : 1-pot (021) number of events",33,-0.5,32.5));
  outputPerPot[2]=std::auto_ptr<TH1D>(new TH1D("outputPerPot022","Track & trigger : 1-pot (022) number of events",33,-0.5,32.5));
  outputPerPot[3]=std::auto_ptr<TH1D>(new TH1D("outputPerPot023","Track & trigger : 1-pot (023) number of events",33,-0.5,32.5));
  outputPerPot[4]=std::auto_ptr<TH1D>(new TH1D("outputPerPot024","Track & trigger : 1-pot (024) number of events",33,-0.5,32.5));
  outputPerPot[5]=std::auto_ptr<TH1D>(new TH1D("outputPerPot025","Track & trigger : 1-pot (025) number of events",33,-0.5,32.5));
  outputPerPot[6]=std::auto_ptr<TH1D>(new TH1D("outputPerPot120","Track & trigger : 1-pot (120) number of events",33,-0.5,32.5));
  outputPerPot[7]=std::auto_ptr<TH1D>(new TH1D("outputPerPot121","Track & trigger : 1-pot (121) number of events",33,-0.5,32.5));
  outputPerPot[8]=std::auto_ptr<TH1D>(new TH1D("outputPerPot122","Track & trigger : 1-pot (122) number of events",33,-0.5,32.5));
  outputPerPot[9]=std::auto_ptr<TH1D>(new TH1D("outputPerPot123","Track & trigger : 1-pot (123) number of events",33,-0.5,32.5));
  outputPerPot[10]=std::auto_ptr<TH1D>(new TH1D("outputPerPot124","Track & trigger : 1-pot (124) number of events",33,-0.5,32.5));
  outputPerPot[11]=std::auto_ptr<TH1D>(new TH1D("outputPerPot125","Track & trigger : 1-pot (125) number of events",33,-0.5,32.5));

  const char *labls[]={"tm(11)","tm(00)","tm(01)","tm(2-3,2-3)","tm(>4,>4)","tm(01)1S","tm(01)1M","tm(01)m!s","tm(00)1S","","","","","","","","1S","0S","0S1M","0S,>1M","0S0M","line u,!v","line v,!u","SCand","MCand","1M",">1M","1S1M","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""};
  int foundBin=0;

  /*
  if (!readMultiTrk)
    {
      labls[20]="";
      labls[21]="";
    }
  */

  for (int tt=0;tt<12;tt++)
    {
      for (int r=0;r<32;r++)
	{
	  foundBin=outputPerPot[tt]->FindBin((float) r);
	  outputPerPot[tt]->GetXaxis()->SetBinLabel(foundBin,labls[r]);
	}
    }


  hisTitle="Track & trigger output percentages "; 
  hisTitle.append(whatCombi[diagonal]);
  hisTitleA="Track & trigger output percentages "; 
  outputPercentagesD=std::auto_ptr<TH1D>(new TH1D("outputPercentagesDiag",hisTitle.c_str(),65,-0.5,64.5));
  outputPercentagesA=std::auto_ptr<TH1D>(new TH1D("outputPercentagesAllEvts",hisTitleA.c_str(),65,-0.5,64.5));
  hisTitle="Track & trigger output (numerator: events fulfilling criteria) "; 
  hisTitle.append(whatCombi[diagonal]);
  hisTitleA="Track & trigger output (numerator: events fulfilling criteria) "; 
  outputFulfillEvtsD=std::auto_ptr<TH1D>(new TH1D("outputFulfillEvtsDiag",hisTitle.c_str(),65,-0.5,64.5));
  outputFulfillEvtsA=std::auto_ptr<TH1D>(new TH1D("outputFulfillEvtsAllEvts",hisTitle.c_str(),65,-0.5,64.5));
  hisTitle="Track & trigger output (denominator: reference events) "; 
  hisTitle.append(whatCombi[diagonal]);
  outputRefEvts=std::auto_ptr<TH1D>(new TH1D("outputRefEvts",hisTitle.c_str(),65,-0.5,64.5));
  /*
  if (diagonal)
    {
      outputPercentages=std::auto_ptr<TH1D>(new TH1D("outputPercentages","Track & trigger output percentages (45bt*56tp)",33,-0.5,32.5));
      outputFulfillEvts=std::auto_ptr<TH1D>(new TH1D("outputFulfillEvts","Track & trigger output (numerator: events fulfilling criteria - 45bt*56tp)",33,-0.5,32.5));
      outputRefEvts=std::auto_ptr<TH1D>(new TH1D("outputRefEvts","Track & trigger output (denominator: reference events - 45bt*56tp)",33,-0.5,32.5));
    }
  else
    {
      outputPercentages=std::auto_ptr<TH1D>(new TH1D("outputPercentages","Track & trigger output percentages (45tp*56bt)",33,-0.5,32.5));
      outputFulfillEvts=std::auto_ptr<TH1D>(new TH1D("outputFulfillEvts","Track & trigger output (numerator: events fulfilling criteria - 45tp*56bt)",33,-0.5,32.5));
      outputRefEvts=std::auto_ptr<TH1D>(new TH1D("outputRefEvts","Track & trigger output (denominator: reference events - 45tp*56bt)",33,-0.5,32.5));
      }
*/
  //  std::auto_ptr<TH1D> ;
  //  std::auto_ptr<TH1D> ;
  //  std::auto_ptr<TH1D> outputRefEvts;


  trksWOneTmcBitOn=std::auto_ptr<TH1D>(new TH1D("trksWOneTmcBitOn","track in reference pot with only TMC bit #x on",16,-0.5,15.5));
  NotrksWOneTmcBitOn=std::auto_ptr<TH1D>(new TH1D("NotrksWOneTmcBitOn","No track in reference pot with only TMC bit #x on",16,-0.5,15.5));

  srcX=std::auto_ptr<TH1D>(new TH1D("srcX","Elastic (TMC) source pot track x0-parameter",1000,-50.,50.));
  srcY=std::auto_ptr<TH1D>(new TH1D("srcY","source pot track y0-parameter",1000,-50.,50.));
  chkY=std::auto_ptr<TH1D>(new TH1D("chkY","checkPot track y0-parameter",1000,-50.,50.));

  std::string relPots[]={"_45_nr","_45_fr","_56_nr","_56_fr"};
  std::string theTitle,theName;
  std::string suffix[]={"T2none","T2same","T2oppo","T2both"};
  std::string suffixT1[]={"T2oppoT1both","T2T1oppo","T2oppoNoT1"};
  //      std::string pots[]={"_45_tp","_45_hr","_45_bt","_56_tp","_56_hr","_56_bt"};
  std::string pots[]={"_45_tp","_45_bt","_45_hr","_56_tp","_56_bt","_56_hr"};
  //  std::string suffixOk[]={"XiRgNONCompat","XiRgCompat"};
  std::string nameYc[]={"noYcut","Ycut"};
  std::string titleYc[]={" - no Y cut - "," - 30>abs(y)>9 mm - "};
   
  
  if (readLoNeg)
    {
      lonegBunch[0]=std::auto_ptr<TH1D>(new TH1D("lonegBunchAny","LoNeg bunch number (no trigger demanded)",4000,-100.5,3999.5));
      lonegBunch[1]=std::auto_ptr<TH1D>(new TH1D("lonegBunchRPV","LoNeg bunch number (rp220v trigger demanded)",4000,-100.5,3999.5));
      lonegBunch[2]=std::auto_ptr<TH1D>(new TH1D("lonegBunchRPH","LoNeg bunch number (rp220h trigger demanded)",4000,-100.5,3999.5));
      lonegBunch[3]=std::auto_ptr<TH1D>(new TH1D("lonegBunchT1","LoNeg bunch number (t1 trigger demanded)",4000,-100.5,3999.5));
      lonegBunch[4]=std::auto_ptr<TH1D>(new TH1D("lonegBunchT2","LoNeg bunch number (t2 trigger demanded)",4000,-100.5,3999.5));
      lonegBunch[5]=std::auto_ptr<TH1D>(new TH1D("lonegBunchBx","LoNeg bunch number (bx trigger demanded)",4000,-100.5,3999.5));
      lonegBunch[6]=std::auto_ptr<TH1D>(new TH1D("lonegBunchRealSD","LoNeg bunch number (SD=2/12(nr+fr) RP trk+T2 tracks on opposite side demanded)",4000,-100.5,3999.5));


      refEvtsBunches.clear();
      
      //      outputPerBunch.clear();
      //      mp4Trk.clear();
      
      mp3Trk.clear();
      
      mp2Trk.clear();
      
      mp1Trk.clear();
      
      mp0Trk.clear();

      //      mp411.clear();
      //      mp400.clear();
      //      mp311t100.clear();
      //      mp311t101.clear();
      
 
      if (readT2)
	{
	  sdTrigs=std::auto_ptr<TH1D>(new TH1D("sdTrigs","LoNeg trigger bits for SD events (2/12 RP Trk (nr+fr)+Opposite side T2 Trk)",10,-0.5,9.5));
	  const char *lablsSD[]={"220v","220h","220x","T2","T1","BX","","","","","","","","","","","","","","","","","","","","","","","","","",""};
	  int foundBinSd=0;
	  for (int r=0;r<10;r++)
	    {
	      foundBinSd=sdTrigs->FindBin((float) r);
	      sdTrigs->GetXaxis()->SetBinLabel(foundBinSd,lablsSD[r]);
	    }


	  //std::map<int,int> 
	  rpvTrig.clear();
	  //std::map<int,int> 
	  rphTrig.clear();
	  //std::map<int,int>
	  bxTrig.clear();
	  //	  std::map<int,int> 
	  t2Trig.clear();
	  //	  std::map<int,int> 
	  sdBunchEv.clear();


	}


    }
				      
  if (readRecoProt)
    {



      Double_t xiBins[105]={-0.04,-0.03,-0.02,-0.01,0};
      double logFactor=1.148153621; //10^(0.06);
      double temp5=0.000001;
      for (int u=5;u<105;u++)
	{
	  xiBins[u]=temp5;
	  temp5=temp5*logFactor;
	}
  
      recoProtXiVsN =std::auto_ptr<TH2D>(new TH2D("recoProtXiVsN","Events with reconstructed protons, numProtons vs Xi (fractional momentum loss);num;Xi",20,-0.5,19.5,1002,-1.005,0.005));
      recoProtXiVsErr =std::auto_ptr<TH2D>(new TH2D("recoProtXiVsErr","Events with reconstructed protons, Xi (fractional momentum loss) vs sqrt(autocovariance(xi,xi)), or -1% if covar<0;Xi;Xi error",425,-0.1,0.07,60,-0.02,0.1));
      
      numRecoPInSD=std::auto_ptr<TH1D>(new TH1D("numRecoPInSD","Number of RecoProtons - Tvals.size (2/12 RP tracks)",11,-0.5,10.5));


	   if (readT2)
	     {
	     
	       //xiRpVsXiRG
	       if (readT1)
		 {
		   xiRpVsXiRG[kT2none]=std::auto_ptr<TH2D>(new TH2D("xiRpVsXiRGT2none","Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) (trk only in T1);MADX;rap.gap",104,xiBins,104,xiBins));
		   xiRpVsXiRG[kT2same]=std::auto_ptr<TH2D>(new TH2D("xiRpVsXiRGT2same","Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) (T2 same-side trks);MADX;rap.gap",104,xiBins,104,xiBins));
		   xiRpVsXiRG[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("xiRpVsXiRGT2oppo","Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) (T2 opposite-side trks);MADX;rap.gap",104,xiBins,104,xiBins));
		   xiRpVsXiRG[kT2both]=std::auto_ptr<TH2D>(new TH2D("xiRpVsXiRGT2both","Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) (T2 both sides trk);MADX;rap.gap",104,xiBins,104,xiBins));
		   
		   
		   RGTrkVsRGxiT1BothT2Oppo=std::auto_ptr<TH1D>(new TH1D("RGTrkVsRGxiT1BothT2Oppo","Rapidity gap between T1 track and proton, minus implied RG from Energy loss xi=-exp(-RG) of each RecoProton (MadX) (T2 opposite-side trks only, but T1 both arms);RG(tr)-RG(xi)",100,-20,20));
		   RGTrk2DRGxiT1BothT2Oppo=std::auto_ptr<TH2D>(new TH2D("RGTrk2DRGxiT1BothT2Oppo","Rapidity gap between T1 track and proton, versus implied RG from Energy loss xi=-exp(-RG) of each RecoProton (MadX) (T2 opposite-side trks only, T1 both arms, 9<abs(y)<30);RG(trT1T2);RG(xiRP)",500,0.,10.,250,0.,20.));
		   
		   RGTrk2DRGxiT1T2Both=std::auto_ptr<TH2D>(new TH2D("RGTrk2DRGxiT1T2Both","Rapidity gap between T1 track and proton, versus implied RG from Energy loss xi=-exp(-RG) of each RecoProton (MadX) (T2 & T1 both arms, 9<abs(y)<30);RG(trT1T2);RG(xiRP)",500,0.,10.,250,0.,20.));
		   
		   for (int u=0;u<3;u++)
		     {
		       for (int v=0;v<6;v++)
			 {
			   theName="YcutXiVsRgXi";
			   theName.append(suffixT1[u]);
			   theName.append(pots[v]);
			   theTitle="Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) - RP trk 30>abs(y)>9mm- T2&T1 condition: ";
			   theTitle.append(suffixT1[u]);
			   theTitle.append(pots[v]);
			   theTitle.append(";MADX;rap.gap");
			   //			   if (!u)
			   //			     xiBadRGT2oppo[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),104,xiBins,104,xiBins));
			   //			   else
			   YcutXiRGT2oppo[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),350,-0.1,0.25,200,-0.001,0.1));

			   theName="noYcutXiVsRgXi";
			   theName.append(suffixT1[u]);
			   theName.append(pots[v]);
			   theTitle="Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) - no Y cut - T2&T1 condition: ";
			   theTitle.append(suffixT1[u]);
			   theTitle.append(pots[v]);
			   theTitle.append(";MADX;rap.gap");
			   //			   if (!u)
			   //			     xiOkRGT2oppo[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),104,xiBins,104,xiBins));
			   //			   else
			   noYcutXiRGT2oppo[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),350,-0.1,0.25,200,-0.001,0.1));


			 }
		     }
		 
		   for (int u=0;u<2;u++)
		     {
		       for (int v=0;v<6;v++)
			 {
			   theName="xiVsRgXiT2T1both";
			   theName.append(nameYc[u]);
			   theName.append(pots[v]);
			   theTitle="Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) - T2&T1 both arms";
			   theTitle.append(titleYc[u]);
			   theTitle.append(pots[v]);
			   theTitle.append(";MADX;rap.gap");
			   //			   if (!u)
			   xiBadOkRGT2both[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),350,-0.1,0.25,200,-0.001,0.25));
			 }
		     }
		   
		   xiRGT2oppo[kT1SameT2Oppo]=std::auto_ptr<TH2D>(new TH2D("xiRGT2oppoButT1Both","Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) (T2 opposite-side trks only, but T1 both arms);MADX;rap.gap",104,xiBins,104,xiBins));
		   xiRGT2oppo[kT1T2Oppo]=std::auto_ptr<TH2D>(new TH2D("xiRGT1T2oppo","Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) (T2 & T1 opposite-side trks only);MADX;rap.gap",200,-0.05,0.05,200,-0.05,0.05));
		   xiRGT2oppo[kNoT1]=std::auto_ptr<TH2D>(new TH2D("xiRGT2oppoNoT1","Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) (T2 opposite-side trks only, no T1 trk);MADX;rap.gap",200,-0.05,0.05,200,-0.05,0.05));
		   
		   xiRGT2both[kT1both]=std::auto_ptr<TH2D>(new TH2D("xiRGT1T2both","Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) (T2 opposite-side trks only, but T1 also/only on p-side);MADX;rap.gap",400,-0.25,0.25,200,-0.25,0.25));
		   xiRGT2both[kT1either]=std::auto_ptr<TH2D>(new TH2D("xiRGT2bothButT1one","Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) (T2 & T1 opposite-side trks only);MADX;rap.gap",400,-0.25,0.25,200,-0.25,0.25));
		   xiRGT2both[kT1neither]=std::auto_ptr<TH2D>(new TH2D("xiRGT2bothNoT1","Energy loss -xi of each RecoProton (MadX) vs -xi(Rapidity gap in T1&T2) (T2 opposite-side trks only, no T1 trk);MADX;rap.gap",400,-0.25,0.25,200,-0.25,0.25));
		   
		 }

	       rpTee[kT2none]=std::auto_ptr<TH1D>(new TH1D("rpTeeT2none","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 no tracks)",1000,0.,10.));
	       rpTee[kT2same]=std::auto_ptr<TH1D>(new TH1D("rpTeeT2same","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trk same side)",1000,0.,10.));
	       rpTee[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("rpTeeT2oppo","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trk opposite side)",1000,0.,10.));
	       rpTee[kT2both]=std::auto_ptr<TH1D>(new TH1D("rpTeeT2both","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trks both sides)",1000,0.,10.));

	       funcTeeVsElastTee[kT2none]=std::auto_ptr<TH2D>(new TH2D("rpTeeVsCalcTeeT2none","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 no tracks) : MadX-function vs 0-xi nominal 90m optics approx;MADX;approx",1000,-0.05,9.95,1000,-0.05,9.95));
	       funcTeeVsElastTee[kT2same]=std::auto_ptr<TH2D>(new TH2D("rpTeeVsCalcTeeT2same","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trk same side) : MadX-function vs 0-xi nominal 90m optics approx;MADX;approx",1000,-0.05,9.95,1000,-0.05,9.95));
	       funcTeeVsElastTee[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("rpTeeVsCalcTeeT2oppo","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trk opposite side) : MadX-function vs 0-xi nominal 90m optics approx;MADX;approx",1000,-0.05,9.95,1000,-0.05,9.95));
	       funcTeeVsElastTee[kT2both]=std::auto_ptr<TH2D>(new TH2D("rpTeeVsCalcTeeT2both","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trks both sides) : MadX-function vs 0-xi nominal 90m optics approx;MADX;approx",1000,-0.05,9.95,1000,-0.05,9.95));

	       funcTeeRelErrVsXi[kT2none]=std::auto_ptr<TH2D>(new TH2D("rpTeeRelErrVsXiT2none","Momentum transfer -t relative error of each RecoProton (2/12 RP tracks+T2 no tracks) : xi vs (MadX-t minus nominal-90m-optics-0xi-approx-t/MadX-t;xi;MADX-approx (%)",520,-1,0.04,100,-100.,100.));
	       funcTeeRelErrVsXi[kT2same]=std::auto_ptr<TH2D>(new TH2D("rpTeeRelErrVsXiT2same","Momentum transfer -t relative error of each RecoProton (2/12 RP tracks+T2 same side tracks) : xi vs (MadX-t minus nominal-90m-optics-0xi-approx-t/MadX-t;xi;MADX-approx (%)",520,-1,0.04,100,-100.,100.));
	       funcTeeRelErrVsXi[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("rpTeeRelErrVsXiT2oppo","Momentum transfer -t relative error of each RecoProton (2/12 RP tracks+T2 opposite side tracks) : xi vs (MadX-t minus nominal-90m-optics-0xi-approx-t/MadX-t;xi;MADX-approx (%)",520,-1,0.04,100,-100.,100.));
	       funcTeeRelErrVsXi[kT2both]=std::auto_ptr<TH2D>(new TH2D("rpTeeRelErrVsXiT2both","Momentum transfer -t relative error of each RecoProton (2/12 RP tracks+T2 both sides tracks) : xi vs (MadX-t minus nominal-90m-optics-0xi-approx-t/MadX-t;xi;MADX-approx (%)",520,-1,0.04,100,-100.,100.));

	       //	       numRecoPInSD[kT2none]=std::auto_ptr<TH1D>(new TH1D("numRecoPInSDT2none","Number of RecoProtons (2/12 RP tracks+T2 no tracks)",11,-0.5,10.5));
	       //	       numRecoPInSD[kT2same]=std::auto_ptr<TH1D>(new TH1D("numRecoPInSDT2same","Number of RecoProtons (2/12 RP tracks+T2 trk same side)",11,-0.5,10.5));
	       //	       numRecoPInSD[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("numRecoPInSDT2oppo","Number of RecoProtons (2/12 RP tracks+T2 trk opposite side)",11,-0.5,10.5));
	       //	       numRecoPInSD[kT2both]=std::auto_ptr<TH1D>(new TH1D("numRecoPInSDT2both","Number of RecoProtons (2/12 RP tracks+T2 trks both sides)",11,-0.5,10.5));

	     

	       ipAngleXRelErr[kT2none]=std::auto_ptr<TH1D>(new TH1D("ipAngleXRecProtRelErrT2none","Theta_X_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 no tracks) : (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-300.,300.));
	       ipAngleXRelErr[kT2same]=std::auto_ptr<TH1D>(new TH1D("ipAngleXRecProtRelErrT2same","Theta_X_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 same side tracks): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-300.,300.));
	       ipAngleXRelErr[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("ipAngleXRecProtRelErrT2oppo","Theta_X_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 opposite side tracks): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-300.,300.));
	       ipAngleXRelErr[kT2both]=std::auto_ptr<TH1D>(new TH1D("ipAngleXRecProtRelErrT2both","Theta_X_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 both sides tracks): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-300.,300.));

	       ipAngleYRelErr[kT2none]=std::auto_ptr<TH1D>(new TH1D("ipAngleYRecProtRelErrT2none","Theta_Y_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 no tracks) : (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-30.,30.));
	       ipAngleYRelErr[kT2same]=std::auto_ptr<TH1D>(new TH1D("ipAngleYRecProtRelErrT2same","Theta_Y_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 same side tracks): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-30.,30.));
	       ipAngleYRelErr[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("ipAngleYRecProtRelErrT2oppo","Theta_Y_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 opposite side tracks): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-30.,30.));
	       ipAngleYRelErr[kT2both]=std::auto_ptr<TH1D>(new TH1D("ipAngleYRecProtRelErrT2both","Theta_Y_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 both sides tracks): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-30.,30.));

	       ElaIpAngleXRelErr[kT2none]=std::auto_ptr<TH1D>(new TH1D("ElaIpAngleXRecProtRelErrT2none","Theta_X_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 no tracks+elastic-box spectroCut) : (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-300.,300.));
	       ElaIpAngleXRelErr[kT2same]=std::auto_ptr<TH1D>(new TH1D("ElaIpAngleXRecProtRelErrT2same","Theta_X_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 same side tracks+elastic-box spectroCut): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-300.,300.));
	       ElaIpAngleXRelErr[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("ElaIpAngleXRecProtRelErrT2oppo","Theta_X_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 opposite side tracks+elastic-box spectroCut): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-300.,300.));
	       ElaIpAngleXRelErr[kT2both]=std::auto_ptr<TH1D>(new TH1D("ElaIpAngleXRecProtRelErrT2both","Theta_X_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 both sides tracks+elastic-box spectroCut): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-300.,300.));

	       ElaIpAngleYRelErr[kT2none]=std::auto_ptr<TH1D>(new TH1D("ElaIpAngleYRecProtRelErrT2none","Theta_Y_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 no tracks+elastic-box spectroCut) : (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-30.,30.));
	       ElaIpAngleYRelErr[kT2same]=std::auto_ptr<TH1D>(new TH1D("ElaIpAngleYRecProtRelErrT2same","Theta_Y_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 same side tracks+elastic-box spectroCut): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-30.,30.));
	       ElaIpAngleYRelErr[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("ElaIpAngleYRecProtRelErrT2oppo","Theta_Y_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 opposite side tracks+elastic-box spectroCut): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-30.,30.));
	       ElaIpAngleYRelErr[kT2both]=std::auto_ptr<TH1D>(new TH1D("ElaIpAngleYRecProtRelErrT2both","Theta_Y_IP for RecoProton : abs(Xi)<0.5% (2/12 RP tracks+T2 both sides tracks+elastic-box spectroCut): (MadX - nominal optics 0-xi-approximation)/MadX;Relative Error (%)",3000,-30.,30.));

	    
	     

	       if (readLoNeg)
		 {

		   if (readT1)
		     {
		     
		       RgXiOkrpTee[0][0]=std::auto_ptr<TH1D>(new TH1D("okRgT2oppoTeeSpectr","-t-value of the proton (RapGap==Xi+-2sigma (consistent), 2/12 RP tracks+T2 opposite side trks & T1 both arms trks, T2 Trig)",1000,0.,10.));
		       RgXiOkrpTee[1][0]=std::auto_ptr<TH1D>(new TH1D("badRgT2oppoTeeSpectr","-t-value of the proton (RapGap=!=Xi+-2sigma (NOT consistent), 2/12 RP tracks+T2 opposite side trks & T1 both arms trks, T2 Trig)",1000,0.,10.));
		       RgXiOkrpTee[0][1]=std::auto_ptr<TH1D>(new TH1D("okRgT2bothTeeSpectr","-t-value of the proton (RapGap==Xi+-2sigma (consistent), 2/12 RP tracks+T2 both arms trks & T1 both arms trks, T2 Trig)",1000,0.,10.));
		       RgXiOkrpTee[1][1]=std::auto_ptr<TH1D>(new TH1D("badRgT2bothTeeSpectr","-t-value of the proton (RapGap=!=Xi+-2sigma (NOT consistent), 2/12 RP tracks+T2 both arms trks & T1 both arms trks, T2 Trig)",1000,0.,10.));
		       okRgAngleXvsX=std::auto_ptr<TH2D>(new TH2D("okRgAngleXvsX","Spectrometer plot (RapGap==Xi+-2sigma): Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near  for RecoProton (2/12 RP tracks+T2 oppo+trig+T1 both arms tracks;Angle;x0",500,-0.001,0.001,800,-40.,40.));

		       badRgAngleXvsX=std::auto_ptr<TH2D>(new TH2D("badRgAngleXvsX","Spectrometer plot (RapGap=!=Xi+-2sigma): Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near  for RecoProton (2/12 RP tracks+T2 oppo+trig+T1 both arms tracks;Angle;x0",500,-0.001,0.001,800,-40.,40.));
		       FiverAngleXvsX=std::auto_ptr<TH2D>(new TH2D("FiverAngleXvsX","Spectrometer plot (xi_MADX > +4.8%): Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near  for RecoProton (2/12 RP tracks+T2 oppo+trig+T1 both arms tracks;Angle;x0",500,-0.001,0.001,800,-40.,40.));

		       FiverXvsY=std::auto_ptr<TH2D>(new TH2D("FiverXvsY","(x,y) for misreconstructed protons (xi_MADX > +4.8%): Events with two tracks, X0_nearvs Y0_near  for RecoProton (2/12 RP tracks+T2 oppo+trig+T1 both arms tracks;x0(mm);y0(mm)",800,-40.,40.,800,-40.,40.));

		       T1T2numTrk[0]=std::auto_ptr<TH1D>(new TH1D("okRgT1numTrk","Number of T1 Tracks(RapGap==Xi+-2sigma) (2/12 RP tracks+T2 opposite side trks & T1 both arms trks, T2 Trig)",80,-0.5,79.5));
		       T1T2numTrk[1]=std::auto_ptr<TH1D>(new TH1D("badRgT1numTrk","Number of T1 Tracks(RapGap=!=Xi+-2sigma) (2/12 RP tracks+T2 opposite side trks & T1 both arms trks, T2 Trig)",80,-0.5,79.5));
		       T1T2numTrk[2]=std::auto_ptr<TH1D>(new TH1D("okRgT2numTrk","Number of T2 Tracks(RapGap==Xi+-2sigma) (2/12 RP tracks+T2 opposite side trks & T1 both arms trks, T2 Trig)",80,-0.5,79.5));
		       T1T2numTrk[3]=std::auto_ptr<TH1D>(new TH1D("badRgT2numTrk","Number of T2 Tracks(RapGap=!=Xi+-2sigma) (2/12 RP tracks+T2 opposite side trks & T1 both arms trks, T2 Trig)",80,-0.5,79.5));
		       rgCollimatorXY[0]=std::auto_ptr<TH2D>(new TH2D("okRgCollimatorXY","Extrapolated the RP 1-pot track backward 20m from near pot towards IP, to presumptive circular apperture (RapGap==Xi+-2sigma) (last collimator, T2 oppos side trks,T1 both);x(mm);y(mm)",300,-60.,60.,300,-60.,60.));
		       rgCollimatorXY[1]=std::auto_ptr<TH2D>(new TH2D("badRgCollimatorXY","Extrapolated the RP 1-pot track backward 20m from near pot towards IP, to presumptive circular apperture (RapGap=!=Xi+-2sigma) (last collimator, T2 oppos side trks,T1 both);x(mm);y(mm)",300,-60.,60.,300,-60.,60.));
		  
		     }

		   FiverThetaYMineVsMadx=std::auto_ptr<TH2D>(new TH2D("FiverThetaYMineVsMadx","Theta_Y_IP for RecoProton;mine;MadX",800,-0.0012,0.0012,800,-0.0004,0.0004));
		   FiverThetaXMineVsMadx=std::auto_ptr<TH2D>(new TH2D("FiverThetaXMineVsMadx","Theta_X_IP for RecoProton;mine;MadX",800,-0.001,0.002,800,-0.001,0.002));
		 
		   T2ipAngleX[kT2none]=std::auto_ptr<TH1D>(new TH1D("t2ipAngleXT2none","Theta_X_IP for RecoProton (2/12 RP tracks+T2 no tracks, but T2 Trig)",8000,-0.004,0.004));
		   T2ipAngleX[kT2same]=std::auto_ptr<TH1D>(new TH1D("t2ipAngleXT2same","Theta_X_IP for RecoProton (2/12 RP tracks+T2 same side tracks and T2 Trig)",8000,-0.004,0.004));
		   T2ipAngleX[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("t2ipAngleXT2oppo","Theta_X_IP for RecoProton (2/12 RP tracks+T2 opposite side tracks and T2 Trig)",8000,-0.004,0.004));
		   T2ipAngleX[kT2both]=std::auto_ptr<TH1D>(new TH1D("t2ipAngleXT2both","Theta_X_IP for RecoProton (2/12 RP tracks+T2 both sides tracks and T2 Trig)",8000,-0.004,0.004));
		   T2ipAngleY[kT2none]=std::auto_ptr<TH1D>(new TH1D("t2ipAngleYT2none","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 no tracks, but T2 Trig)",8000,-0.004,0.004));
		   T2ipAngleY[kT2same]=std::auto_ptr<TH1D>(new TH1D("t2ipAngleYT2same","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 same side tracks and T2 Trig)",8000,-0.004,0.004));
		   T2ipAngleY[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("t2ipAngleYT2oppo","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 opposite side tracks and T2 Trig)",8000,-0.004,0.004));
		   T2ipAngleY[kT2both]=std::auto_ptr<TH1D>(new TH1D("t2ipAngleYT2both","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 both sides tracks and T2 Trig)",8000,-0.004,0.004));
		   
		   ElaT2ipAngleX[kT2none]=std::auto_ptr<TH1D>(new TH1D("t2ElaIpAngleXT2none","Theta_X_IP for RecoProton (2/12 RP tracks+T2 no tracks, but T2 Trig&Elastic spectroBox)",8000,-0.004,0.004));
		   ElaT2ipAngleX[kT2same]=std::auto_ptr<TH1D>(new TH1D("t2ElaIpAngleXT2same","Theta_X_IP for RecoProton (2/12 RP tracks+T2 same side tracks and T2 Trig&Elastic spectroBox)",8000,-0.004,0.004));
		   ElaT2ipAngleX[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("t2ElaIpAngleXT2oppo","Theta_X_IP for RecoProton (2/12 RP tracks+T2 opposite side tracks and T2 Trig&Elastic spectroBox)",8000,-0.004,0.004));
		   ElaT2ipAngleX[kT2both]=std::auto_ptr<TH1D>(new TH1D("t2ElaIpAngleXT2both","Theta_X_IP for RecoProton (2/12 RP tracks+T2 both sides tracks and T2 Trig&Elastic spectroBox)",8000,-0.004,0.004));
		   ElaT2ipAngleY[kT2none]=std::auto_ptr<TH1D>(new TH1D("t2ElaIpAngleYT2none","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 no tracks, but T2 Trig&Elastic spectroBox)",8000,-0.004,0.004));
		   ElaT2ipAngleY[kT2same]=std::auto_ptr<TH1D>(new TH1D("t2ElaIpAngleYT2same","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 same side tracks and T2 Trig&Elastic spectroBox)",8000,-0.004,0.004));
		   ElaT2ipAngleY[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("t2ElaIpAngleYT2oppo","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 opposite side tracks and T2 Trig&Elastic spectroBox)",8000,-0.004,0.004));
		   ElaT2ipAngleY[kT2both]=std::auto_ptr<TH1D>(new TH1D("t2ElaIpAngleYT2both","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 both sides tracks and T2 Trig&Elastic spectroBox)",8000,-0.004,0.004));
		   
		   SdT2ipAngleX[kT2none]=std::auto_ptr<TH1D>(new TH1D("t2SdIpAngleXT2none","Theta_X_IP for RecoProton (2/12 RP tracks+T2 no tracks, but T2 Trig&SD spectroBox)",8000,-0.004,0.004));
		   SdT2ipAngleX[kT2same]=std::auto_ptr<TH1D>(new TH1D("t2SdIpAngleXT2same","Theta_X_IP for RecoProton (2/12 RP tracks+T2 same side tracks and T2 Trig&SD spectroBox)",8000,-0.004,0.004));
		   SdT2ipAngleX[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("t2SdIpAngleXT2oppo","Theta_X_IP for RecoProton (2/12 RP tracks+T2 opposite side tracks and T2 Trig&SD spectroBox)",8000,-0.004,0.004));
		   SdT2ipAngleX[kT2both]=std::auto_ptr<TH1D>(new TH1D("t2SdIpAngleXT2both","Theta_X_IP for RecoProton (2/12 RP tracks+T2 both sides tracks and T2 Trig&SD spectroBox)",8000,-0.004,0.004));
		   SdT2ipAngleY[kT2none]=std::auto_ptr<TH1D>(new TH1D("t2SdIpAngleYT2none","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 no tracks, but T2 Trig&SD spectroBox)",8000,-0.004,0.004));
		   SdT2ipAngleY[kT2same]=std::auto_ptr<TH1D>(new TH1D("t2SdIpAngleYT2same","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 same side tracks and T2 Trig&SD spectroBox)",8000,-0.004,0.004));
		   SdT2ipAngleY[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("t2SdIpAngleYT2oppo","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 opposite side tracks and T2 Trig&SD spectroBox)",8000,-0.004,0.004));
		   SdT2ipAngleY[kT2both]=std::auto_ptr<TH1D>(new TH1D("t2SdIpAngleYT2both","Theta_Y_IP for RecoProton (2/12 RP tracks+T2 both sides tracks and T2 Trig&SD spectroBox)",8000,-0.004,0.004));
		   
		   
		   t2rpTee[kT2none]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2none","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 no tracks but trig)",5000,0.,10.));
		   t2rpTee[kT2same]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2same","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trk same side and trig)",5000,0.,10.));
		   t2rpTee[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2oppo","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trk opposite side and trig)",5000,0.,10.));
		   t2rpTee[kT2both]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2both","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trks both sides and trig)",5000,0.,10.));
		   ElastT2rpTee[kT2none]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2noneElBox","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 no tracks but trig AND (x,angle) within elastic/low-xi box)",5000,0.,10.));
		   ElastT2rpTee[kT2same]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2sameElBox","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trk same side and trig AND (x,angle) within elastic/low-xi box)",5000,0.,10.));
		   ElastT2rpTee[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2oppoElBox","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trk opposite side and trig AND (x,angle) within elastic/low-xi box)",5000,0.,10.));
		   ElastT2rpTee[kT2both]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2bothElBox","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trks both sides and trig AND (x,angle) within elastic/low-xi box)",5000,0.,10.));
		   
		   SdT2rpTee[kT2none]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2noneSdBox","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 no tracks but trig AND x GT -2sigma=low-xi SD limit)",5000,0.,10.));
		   SdT2rpTee[kT2same]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2sameSdBox","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trk same side and trig AND x GT -2sigma=low-xi SD limit)",5000,0.,10.));
		   SdT2rpTee[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2oppoSdBox","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trk opposite side and trig AND x GT -2sigma=low-xi SD limit)",5000,0.,10.));
		   SdT2rpTee[kT2both]=std::auto_ptr<TH1D>(new TH1D("T2TrgRpTeeT2bothSdBox","Momentum transfer -t of each RecoProton (2/12 RP tracks+T2 trks both sides and trig AND x GT -2sigma=low-xi SD limit)",5000,0.,10.));


		   for (int u=0;u<4;u++)
		     {
		       for (int v=0;v<6;v++)
			 {
			   theName="t2Tee";
			   theName.append(suffix[u]);
			   theName.append(pots[v]);
			   theTitle="Momentum transfer -t (track in one pair near+far only AND T2 trig) - T2&RP condition: ";
			   theTitle.append(suffix[u]);
			   theTitle.append(pots[v]);
			   theTitle.append(";-t(GeV**2)");
			   t2rpTeePair[u][v]=std::auto_ptr<TH1D>(new TH1D(theName.c_str(),theTitle.c_str(),5000,0.,10.));

			   theName="SDt2Tee";
			   theName.append(suffix[u]);
			   theName.append(pots[v]);
			   theTitle="Momentum transfer -t (track in one pair near+far only, AND T2 Trig AND SD-x-cut applied) - T2&RP condition: ";
			   theTitle.append(suffix[u]);
			   theTitle.append(pots[v]);
			   theTitle.append(";-t(GeV**2)");
			   SDt2rpTeePair[u][v]=std::auto_ptr<TH1D>(new TH1D(theName.c_str(),theTitle.c_str(),5000,0.,10.));


			   //" (T2 empty, T2 triggered);Angle;x0",5000,-0.01,0.01,1000,-50.,50.));
			   theName="T2SDSpectroQuart";
			   theName.append(suffix[u]);
			   theName.append(pots[v]);
			   theTitle="RP Delta X0 (far - near) / abs delta Z0 vs X0_near (track in one pair near+far only, AND T2 Trig applied) - T2&RP condition: ";
			   theTitle.append(suffix[u]);
			   theTitle.append(pots[v]);
			   theTitle.append(";xAngle(rad);x0(mm)");
			   T2SpectroQuarter[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),500,-0.001,0.001,1000,-50.,50));

			   theName="bbT2SDSpectroQuart";
			   theName.append(suffix[u]);
			   theName.append(pots[v]);
			   theTitle="RP Delta X0 (far - near) / abs delta Z0 vs X0_near (track in one pair near+far only, AND T2 Trig applied + big bunch) - T2&RP condition: ";
			   theTitle.append(suffix[u]);
			   theTitle.append(pots[v]);
			   theTitle.append(";xAngle(rad);x0(mm)");
			   T2SpectroQuarterBB[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),500,-0.001,0.001,1000,-50.,50));


			   theName="obT2SDSpectroQuart";
			   theName.append(suffix[u]);
			   theName.append(pots[v]);
			   theTitle="RP Delta X0 (far - near) / abs delta Z0 vs X0_near (track in one pair near+far only, AND T2 Trig applied + not big bunch) - T2&RP condition: ";
			   theTitle.append(suffix[u]);
			   theTitle.append(pots[v]);
			   theTitle.append(";xAngle(rad);x0(mm)");
			   T2SpectroQuarterOB[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),500,-0.001,0.001,1000,-50.,50));
			 }
		     }



		   ElastT2rpTeeVsXi[kT2none]=std::auto_ptr<TH2D>(new TH2D("T2TrgRpTeeVsXiT2noneElBox","Momentum transfer -t of each RecoProton vs its xi (2/12 RP tracks+T2 no tracks but trig AND (x,angle) within elastic/low-xi box);t;xi",1000,0.,100.,520,-1.,0.04));
		   ElastT2rpTeeVsXi[kT2same]=std::auto_ptr<TH2D>(new TH2D("T2TrgRpTeeVsXiT2sameElBox","Momentum transfer -t of each RecoProton vs its xi (2/12 RP tracks+T2 trk same side and trig AND (x,angle) within elastic/low-xi box);t;xi",1000,0.,100.,520,-1.,0.04));
		   ElastT2rpTeeVsXi[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("T2TrgRpTeeVsXiT2oppoElBox","Momentum transfer -t of each RecoProton vs its xi (2/12 RP tracks+T2 trk opposite side and trig AND (x,angle) within elastic/low-xi box);t;xi",1000,0.,100.,520,-1.,0.04));
		   ElastT2rpTeeVsXi[kT2both]=std::auto_ptr<TH2D>(new TH2D("T2TrgRpTeeVsXiT2bothElBox","Momentum transfer -t of each RecoProton vs its xi (2/12 RP tracks+T2 trks both sides and trig AND (x,angle) within elastic/low-xi box);t;xi",1000,0.,100.,520,-1.,0.04));

		   t2rpTeeVsXi[kT2none]=std::auto_ptr<TH2D>(new TH2D("T2TrgRpTeeVsXiT2none","Momentum transfer -t of each RecoProton vs its xi (2/12 RP tracks+T2 no tracks but trig;t;xi",1000,0.,100.,520,-1.,0.04));
		   t2rpTeeVsXi[kT2same]=std::auto_ptr<TH2D>(new TH2D("T2TrgRpTeeVsXiT2same","Momentum transfer -t of each RecoProton vs its xi (2/12 RP tracks+T2 trk same side and trig;t;xi",1000,0.,100.,520,-1.,0.04));
		   t2rpTeeVsXi[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("T2TrgRpTeeVsXiT2oppo","Momentum transfer -t of each RecoProton vs its xi (2/12 RP tracks+T2 trk opposite side and trig;t;xi",1000,0.,100.,520,-1.,0.04));
		   t2rpTeeVsXi[kT2both]=std::auto_ptr<TH2D>(new TH2D("T2TrgRpTeeVsXiT2both","Momentum transfer -t of each RecoProton vs its xi (2/12 RP tracks+T2 trks both sides and trig;t;xi",1000,0.,100.,520,-1.,0.04));

		 }
	     }

    }


  if (readT2)
    {
      SDtrkX[kT2none]=std::auto_ptr<TH1D>(new TH1D("sdXT2none","2-pot track x0-parameter (T2 empty)",800,-40.,40.));
      SDtrkX[kT2same]=std::auto_ptr<TH1D>(new TH1D("sdXT2same","2-pot track x0-parameter (T2 trk same side)",800,-40.,40.));
      SDtrkX[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("sdXT2oppo","2-pot track x0-parameter (T2 trk opposite side)",800,-40.,40.));
      SDtrkX[kT2both]=std::auto_ptr<TH1D>(new TH1D("sdXT2both","2-pot track x0-parameter (T2 trks both sides)",800,-40.,40.));

      SDtrkY[kT2none]=std::auto_ptr<TH1D>(new TH1D("sdYT2none","2-pot track y0-parameter (T2 empty)",800,-40.,40.));
      SDtrkY[kT2same]=std::auto_ptr<TH1D>(new TH1D("sdYT2same","2-pot track y0-parameter (T2 trk same side)",800,-40.,40.));
      SDtrkY[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("sdYT2oppo","2-pot track y0-parameter (T2 trk opposite side)",800,-40.,40.));
      SDtrkY[kT2both]=std::auto_ptr<TH1D>(new TH1D("sdYT2both","2-pot track y0-parameter (T2 trks both sides)",800,-40.,40.));

      SDT2Multi[0]=std::auto_ptr<TH1D>(new TH1D("sdMulT2same","T2 multiplicity, same side",500,0.,500.));
      SDT2Multi[1]=std::auto_ptr<TH1D>(new TH1D("sdMulT2oppo","T2 multiplicity, opposite side",500,0.,500.));
      SDT2Multi[2]=std::auto_ptr<TH1D>(new TH1D("sdMulT2bothO","T2 multiplicity, both sides (opposite side)",500,0.,500.));
      SDT2Multi[3]=std::auto_ptr<TH1D>(new TH1D("sdMulT2bothS","T2 multiplicity, both sides (same side)",500,0.,500.));


      allTrkCollimatorXY[kT2none]=std::auto_ptr<TH2D>(new TH2D("allTrkCollimatorXYNoT2","Extrapolated the RP 1-pot track backward 20m from near pot towards IP, to presumptive circular apperture (last collimator, T2 no trk);x(mm);y(mm)",300,-60.,60.,300,-60.,60.));
      allTrkCollimatorXY[kT2same]=std::auto_ptr<TH2D>(new TH2D("allTrkCollimatorXYSameT2","Extrapolated the RP 1-pot track backward 20m from near pot towards IP, to presumptive circular apperture (last collimator, T2 same side trks);x(mm);y(mm)",300,-60.,60.,300,-60.,60.));
      allTrkCollimatorXY[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("allTrkCollimatorXYOppoT2","Extrapolated the RP 1-pot track backward 20m from near pot towards IP, to presumptive circular apperture (last collimator, T2 oppos side trks);x(mm);y(mm)",300,-60.,60.,300,-60.,60.));
      allTrkCollimatorXY[kT2both]=std::auto_ptr<TH2D>(new TH2D("allTrkCollimatorXYBothT2","Extrapolated the RP 1-pot track backward 20m from near pot towards IP, to presumptive circular apperture (last collimator, T2 both sides trks);x(mm);y(mm)",300,-60.,60.,300,-60.,60.));

      T2FitEtaVs000Eta=std::auto_ptr<TH2D>(new TH2D("T2FitEtaVs000Eta","T1T2Track: average eta calculated from T2Hits vs Eta() from track fit;T2Hit eta;Method Eta()",100,-10.,10.,100,-10.,10.));
      T2Eta000Dispersion=std::auto_ptr<TH1D>(new TH1D("T2Eta000Dispersion","T1T2Track: eta calculated from each T2Hit minus Eta() from track fit;eta dispersion",200,-10.,10.));

      T1FitEtaVs000Eta=std::auto_ptr<TH2D>(new TH2D("T1FitEtaVs000Eta","T1T2Track: average eta calculated from T1RecHitGlobals vs Eta() from track fit;T1RecHitGlobal eta;Method Eta()",100,-10.,10.,100,-10.,10.));
      T1Eta000Dispersion=std::auto_ptr<TH1D>(new TH1D("T1Eta000Dispersion","T1T2Track: eta calculated from each T1RecHitGlobal minus Eta() from track fit;eta dispersion",200,-10.,10.));



      if (readClusters)
	{
	  TwoRPTrkSumCluRest[kT2none]=std::auto_ptr<TH1D>(new TH1D("TwoRPTrkSumCluRestT2None","Sum of clusters in other 10 pots, when exactly 2 pots have tracks, near+far (T2 no tracks)",5000,-0.5,4999.5));
	  TwoRPTrkSumCluRest[kT2same]=std::auto_ptr<TH1D>(new TH1D("TwoRPTrkSumCluRestT2Same","Sum of clusters in other 10 pots, when exactly 2 pots have tracks, near+far (T2 tracks only on same side)",5000,-0.5,4999.5));
	  TwoRPTrkSumCluRest[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("TwoRPTrkSumCluRestT2Opp","Sum of clusters in other 10 pots, when exactly 2 pots have tracks, near+far (T2 tracks only on opposite side)",5000,-0.5,4999.5));
	  TwoRPTrkSumCluRest[kT2both]=std::auto_ptr<TH1D>(new TH1D("TwoRPTrkSumCluRestT2Both","Sum of clusters in other 10 pots, when exactly 2 pots have tracks, near+far (T2 tracks both sides)",5000,-0.5,4999.5));


	  if (readLoNeg)
	    {
	      SDCluNoTrk[kT2none]=std::auto_ptr<TH1D>(new TH1D("sd10potCluT2none","10 pot without a track : sum of clusters (+=same side as proton, -=other side,+/-0=+/-0.1) (T2 empty)",50000,-5000.,5000.));
	      SDCluNoTrk[kT2same]=std::auto_ptr<TH1D>(new TH1D("sd10potCluT2same","10 pot without a track : sum of clusters (+=same side as proton, -=other side,+/-0=+/-0.1) (T2 trk only same side)",50000,-5000.,5000.));
	      SDCluNoTrk[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("sd10potCluT2oppo","10 pot without a track : sum of clusters (+=same side as proton, -=other side,+/-0=+/-0.1) (T2 trk only opposite side)",50000,-5000.,5000.));
	      SDCluNoTrk[kT2both]=std::auto_ptr<TH1D>(new TH1D("sd10potCluT2both","10 pot without a track : sum of clusters (+=same side as proton, -=other side,+/-0=+/-0.1) (T2 trk both sides)",50000,-5000.,5000.));

	      rpTrgSDCluNoTrk[kT2none]=std::auto_ptr<TH1D>(new TH1D("sd10potCluRPTrigT2none","RP-V triggered events : 10 pot without a track : sum of clusters (+=same side as proton, -=other side,+/-0=+/-0.1) (T2 empty)",50000,-5000.,5000.));
	      rpTrgSDCluNoTrk[kT2same]=std::auto_ptr<TH1D>(new TH1D("sd10potCluRPTrigT2same","RP-V triggered events : 10 pot without a track : sum of clusters (+=same side as proton, -=other side,+/-0=+/-0.1) (T2 trk only same side)",50000,-5000.,5000.));
	      rpTrgSDCluNoTrk[kT2oppo]=std::auto_ptr<TH1D>(new TH1D("sd10potCluRPTrigT2oppo","RP-V triggered events : 10 pot without a track : sum of clusters (+=same side as proton, -=other side,+/-0=+/-0.1) (T2 trk only opposite side)",50000,-5000.,5000.));
	      rpTrgSDCluNoTrk[kT2both]=std::auto_ptr<TH1D>(new TH1D("sd10potCluRPTrigT2both","RP-V triggered events : 10 pot without a track : sum of clusters (+=same side as proton, -=other side,+/-0=+/-0.1) (T2 trk both sides)",50000,-5000.,5000.));


	      otherSideCluVsMaxTMC[kT2none] =std::auto_ptr<TH2D>(new TH2D("otherSideCluVsMaxTMCT2none","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs U,V-max out of 6 pots in other arm (T2 no trk);U/V;clusters",17,-0.5,16.5,14,CluBins));
	      otherSideCluVsMaxTMC[kT2same] =std::auto_ptr<TH2D>(new TH2D("otherSideCluVsMaxTMCT2same","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs U,V-max out of 6 pots in other arm (T2 trk same side);U/V;clusters",17,-0.5,16.5,14,CluBins));
	      otherSideCluVsMaxTMC[kT2oppo] =std::auto_ptr<TH2D>(new TH2D("otherSideCluVsMaxTMCT2oppo","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs U,V-max out of 6 pots in other arm (T2 trk opposite side);U/V;clusters",17,-0.5,16.5,14,CluBins));
	      otherSideCluVsMaxTMC[kT2both] =std::auto_ptr<TH2D>(new TH2D("otherSideCluVsMaxTMCT2both","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs U,V-max out of 6 pots in other arm (T2 trk both sides);U/V;clusters",17,-0.5,16.5,14,CluBins));
	      otherSideCluVsCluPerPot[kT2none] =std::auto_ptr<TH2D>(new TH2D("otherSideCluVsCluPerPotT2none","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 no trk);clu/pot;clusters",17,CluBins2,14,CluBins));
	      otherSideCluVsCluPerPot[kT2same] =std::auto_ptr<TH2D>(new TH2D("otherSideCluVsCluPerPotT2same","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk same side);clu/pot;clusters",17,CluBins2,14,CluBins));
	      otherSideCluVsCluPerPot[kT2oppo] =std::auto_ptr<TH2D>(new TH2D("otherSideCluVsCluPerPotT2oppo","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk opposite side);clu/pot;clusters",17,CluBins2,14,CluBins));
	      otherSideCluVsCluPerPot[kT2both] =std::auto_ptr<TH2D>(new TH2D("otherSideCluVsCluPerPotT2both","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk both sides);clu/pot;clusters",17,CluBins2,14,CluBins));
	      BNon0OtherSideClu[kT2none] =std::auto_ptr<TH2D>(new TH2D("BNon0OtherSideCluVsCluPerPotT2none","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 no trk, bunch!=0);clu/pot;clusters",17,CluBins2,14,CluBins));
	      BNon0OtherSideClu[kT2same] =std::auto_ptr<TH2D>(new TH2D("BNon0OtherSideCluVsCluPerPotT2same","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk same side, bunch!=0);clu/pot;clusters",17,CluBins2,14,CluBins));
	      BNon0OtherSideClu[kT2oppo] =std::auto_ptr<TH2D>(new TH2D("BNon0OtherSideCluVsCluPerPotT2oppo","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk opposite side, bunch!=0);clu/pot;clusters",17,CluBins2,14,CluBins));
	      BNon0OtherSideClu[kT2both] =std::auto_ptr<TH2D>(new TH2D("BNon0OtherSideCluVsCluPerPotT2both","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk both sides, bunch!=0);clu/pot;clusters",17,CluBins2,14,CluBins));

	      RPTrigOtherSideCluVsCluPerPot[kT2none] =std::auto_ptr<TH2D>(new TH2D("RPTrigOtherSideCluVsCluPerPotT2none","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 no trk);clu/pot;clusters",17,CluBins2,14,CluBins));
	      RPTrigOtherSideCluVsCluPerPot[kT2same] =std::auto_ptr<TH2D>(new TH2D("RPTrigOtherSideCluVsCluPerPotT2same","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk same side);clu/pot;clusters",17,CluBins2,14,CluBins));
	      RPTrigOtherSideCluVsCluPerPot[kT2oppo] =std::auto_ptr<TH2D>(new TH2D("RPTrigOtherSideCluVsCluPerPotT2oppo","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk opposite side);clu/pot;clusters",17,CluBins2,14,CluBins));
	      RPTrigOtherSideCluVsCluPerPot[kT2both] =std::auto_ptr<TH2D>(new TH2D("RPTrigOtherSideCluVsCluPerPotT2both","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk both sides);clu/pot;clusters",17,CluBins2,14,CluBins));


	      BNon0RPTrigOtherSideClu[kT2none] =std::auto_ptr<TH2D>(new TH2D("BNon0RPTrigOtherSideCluVsCluPerPotT2none","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 no trk, bunch!=0);clu/pot;clusters",17,CluBins2,14,CluBins));
	      BNon0RPTrigOtherSideClu[kT2same] =std::auto_ptr<TH2D>(new TH2D("BNon0RPTrigOtherSideCluVsCluPerPotT2same","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk same side, bunch!=0);clu/pot;clusters",17,CluBins2,14,CluBins));
	      BNon0RPTrigOtherSideClu[kT2oppo] =std::auto_ptr<TH2D>(new TH2D("BNon0RPTrigOtherSideCluVsCluPerPotT2oppo","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk opposite side, bunch!=0);clu/pot;clusters",17,CluBins2,14,CluBins));
	      BNon0RPTrigOtherSideClu[kT2both] =std::auto_ptr<TH2D>(new TH2D("BNon0RPTrigOtherSideCluVsCluPerPotT2both","Demanded 2/12 pots with tracks, plotted total clusters in the other, empty, arm vs clusters per pot in other arm (T2 trk both sides, bunch!=0);clu/pot;clusters",17,CluBins2,14,CluBins));


	      RPTrigOtherSideCluVsMaxTMC[kT2none] =std::auto_ptr<TH2D>(new TH2D("RPTrigOtherSideCluVsMaxTMCT2none","Demanded 2/12 pots with tracks AND RP-V trigbit, plotted total clusters in the other, empty, arm vs U,V-max out of 6 pots in other arm (T2 no trk);U/V;clusters",17,-0.5,16.5,14,CluBins));
	      RPTrigOtherSideCluVsMaxTMC[kT2same] =std::auto_ptr<TH2D>(new TH2D("RPTrigOtherSideCluVsMaxTMCT2same","Demanded 2/12 pots with tracksAND RP-V trigbit, plotted total clusters in the other, empty, arm vs U,V-max out of 6 pots in other arm (T2 trk same side);U/V;clusters",17,-0.5,16.5,14,CluBins));
	      RPTrigOtherSideCluVsMaxTMC[kT2oppo] =std::auto_ptr<TH2D>(new TH2D("RPTrigOtherSideCluVsMaxTMCT2oppo","Demanded 2/12 pots with tracksAND RP-V trigbit, plotted total clusters in the other, empty, arm vs U,V-max out of 6 pots in other arm (T2 trk opposite side);U/V;clusters",17,-0.5,16.5,14,CluBins));
	      RPTrigOtherSideCluVsMaxTMC[kT2both] =std::auto_ptr<TH2D>(new TH2D("RPTrigOtherSideCluVsMaxTMCT2both","Demanded 2/12 pots with tracksAND RP-V trigbit, plotted total clusters in the other, empty, arm vs U,V-max out of 6 pots in other arm (T2 trk both sides);U/V;clusters",17,-0.5,16.5,14,CluBins));

	      
	      //	      const char *lablsUV[]={"00","01","11","02","12","22","03","13","23","33","","","","","","","1S","0S","0S1M","0S,>1M","0S0M","line u,!v","line v,!u","SCand","MCand","1M",">1M","1S1M","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""};
	      //	      for (int u=0;u<4;u++)


	    }
	}

 

      SD2pAngleXvsX[kT2none]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2none","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 empty);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      SD2pAngleXvsX[kT2same]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2same","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk same side);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      SD2pAngleXvsX[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2oppo","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk opposite side);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      SD2pAngleXvsX[kT2both]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2both","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk both sides);Angle;x0",500,-0.001,0.001,800,-40.,40.));
 

      T2TrigSD2pAngleXvsX[kT2none]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2noneT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 empty, T2 triggered);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2TrigSD2pAngleXvsX[kT2same]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2sameT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk same side, T2 triggered);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2TrigSD2pAngleXvsX[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2oppoT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk opposite side, T2 triggered);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2TrigSD2pAngleXvsX[kT2both]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2bothT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk both sides, T2 triggered);Angle;x0",500,-0.001,0.001,800,-40.,40.));

      T2TrigSD2pAngleXvsXBB[kT2none]=std::auto_ptr<TH2D>(new TH2D("BBsdAngleXXT2noneT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 empty, T2 triggered,big bunch);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2TrigSD2pAngleXvsXBB[kT2same]=std::auto_ptr<TH2D>(new TH2D("BBsdAngleXXT2sameT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk same side, T2 triggered,big bunch);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2TrigSD2pAngleXvsXBB[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("BBsdAngleXXT2oppoT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk opposite side, T2 triggered, big bunch);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2TrigSD2pAngleXvsXBB[kT2both]=std::auto_ptr<TH2D>(new TH2D("BBsdAngleXXT2bothT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk both sides, T2 triggered, big bunch);Angle;x0",500,-0.001,0.001,800,-40.,40.));

      T2TrigSD2pAngleXvsXOB[kT2none]=std::auto_ptr<TH2D>(new TH2D("OBsdAngleXXT2noneT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 empty, T2 triggered, not big bunch);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2TrigSD2pAngleXvsXOB[kT2same]=std::auto_ptr<TH2D>(new TH2D("OBsdAngleXXT2sameT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk same side, T2 triggered, not big bunch);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2TrigSD2pAngleXvsXOB[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("OBsdAngleXXT2oppoT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk opposite side, T2 triggered, not big bunch);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2TrigSD2pAngleXvsXOB[kT2both]=std::auto_ptr<TH2D>(new TH2D("OBsdAngleXXT2bothT2Trig","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk both sides, T2 triggered, not big bunch);Angle;x0",500,-0.001,0.001,800,-40.,40.));
 
    
      T2noRPTrigSD2pAngleXvsX[kT2none]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2noneT2TrigRPnoTrg","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 empty, T2 & !RP triggered);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2noRPTrigSD2pAngleXvsX[kT2same]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2sameT2TrigRPnoTrg","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk same side, T2 & !RP triggered);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2noRPTrigSD2pAngleXvsX[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2oppoT2TrigRPnoTrg","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk opposite side, T2 & !RP triggered);Angle;x0",500,-0.001,0.001,800,-40.,40.));
      T2noRPTrigSD2pAngleXvsX[kT2both]=std::auto_ptr<TH2D>(new TH2D("sdAngleXXT2bothT2TrigRPnoTrg","Events with two tracks, delta X0 (far - near) / abs delta Z0 vs X0_near (T2 trk both sides, T2 & !RP triggered);Angle;x0",500,-0.001,0.001,800,-40.,40.));
 
     
      SD2pAngleYvsY[kT2none]=std::auto_ptr<TH2D>(new TH2D("sdAngleYYT2none","Events with two tracks, delta Y0 (far - near) / abs delta Z0 vs Y0_near (T2 empty);Angle;y0",500,-0.001,0.001,800,-40.,40.));
      SD2pAngleYvsY[kT2same]=std::auto_ptr<TH2D>(new TH2D("sdAngleYYT2same","Events with two tracks, delta Y0 (far - near) / abs delta Z0 vs Y0_near (T2 trk same side);Angle;y0",500,-0.001,0.001,800,-40.,40.));
      SD2pAngleYvsY[kT2oppo]=std::auto_ptr<TH2D>(new TH2D("sdAngleYYT2oppo","Events with two tracks, delta Y0 (far - near) / abs delta Z0 vs Y0_near (T2 trk opposite side);Angle;y0",500,-0.001,0.001,800,-40.,40.));
      SD2pAngleYvsY[kT2both]=std::auto_ptr<TH2D>(new TH2D("sdAngleYYT2both","Events with two tracks, delta Y0 (far - near) / abs delta Z0 vs Y0_near (T2 trk both sides);Angle;y0",500,-0.001,0.001,800,-40.,40.));


      for (int u=0;u<4;u++)
	{
	  theName="sd2TrkWhichPot";
	  theName.append(suffix[u]);
	  theTitle="Which pots have a track (track in one pair near+far only) - T2 condition: ";
	  theTitle.append(suffix[u]);
	  SDTrkWhichPot[u]=std::auto_ptr<TH1D>(new TH1D(theName.c_str(),theTitle.c_str(),128,-0.5,127.5));
	  theName="sd1TrkWhichPot";
	  theName.append(suffix[u]);
	  theTitle="Which pots have a track (track in one pot only) - T2 condition: ";
	  theTitle.append(suffix[u]);
	  SDBadTrkWhichPot[u]=std::auto_ptr<TH1D>(new TH1D(theName.c_str(),theTitle.c_str(),128,-0.5,127.5));
	  for (int v=0;v<6;v++)
	    {
	      theName="sd2TrkXNearVsFar";
	      theName.append(suffix[u]);
	      theName.append(pots[v]);
	      theTitle="X_near vs X_far (track in one pair near+far only) - T2+RP condition: ";
	      theTitle.append(suffix[u]);
	      theTitle.append(pots[v]);
	      theTitle.append(";far;near");
	      SDXNearVsFar[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),800,-40.,40.,800,-40.,40.));

	      theName="sd2TrkYNearVsFar";
	      theName.append(suffix[u]);
	      theName.append(pots[v]);
	      theTitle="Y_near vs Y_far (track in one pair near+far only) - T2+RP condition: ";
	      theTitle.append(suffix[u]);
	      theTitle.append(pots[v]);
	      theTitle.append(";far;near");
	      SDYNearVsFar[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),800,-40.,40.,800,-40.,40.));

	      theName="sd2TrkXNearMinusFar";
	      theName.append(suffix[u]);
	      theName.append(pots[v]);
	      theTitle="X_near minus X_far (track in one pair near+far only) - T2+RP condition: ";
	      theTitle.append(suffix[u]);
	      theTitle.append(pots[v]);
	      SDXNearMinusFar[u][v]=std::auto_ptr<TH1D>(new TH1D(theName.c_str(),theTitle.c_str(),1000,-100.,100.));

	      theName="sd2TrkYNearMinusFar";
	      theName.append(suffix[u]);
	      theName.append(pots[v]);
	      theTitle="Y_near minus Y_far (track in one pair near+far only) - T2+RP condition: ";
	      theTitle.append(suffix[u]);
	      theTitle.append(pots[v]);
	      SDYNearMinusFar[u][v]=std::auto_ptr<TH1D>(new TH1D(theName.c_str(),theTitle.c_str(),1000,-100.,100.));
	    }
	}

 

      for (int u=0;u<4;u++)
	{
	  for (int v=0;v<4;v++)
	    {
	      theName="sdXY";
	      theName.append(suffix[u]);
	      theName.append(relPots[v]);
	      theTitle="X0 vs Y0 (track in one pair near+far only) - T2&RP condition: ";
	      theTitle.append(suffix[u]);
	      theTitle.append(relPots[v]);
	      theTitle.append(";x;y");
	      SDXvsY[u][v]=std::auto_ptr<TH2D>(new TH2D(theName.c_str(),theTitle.c_str(),800,-40.,40.,800,-40.,40.));
	    }
	}
      /*
      SDXvsY[0][0]=std::auto_ptr<TH2D>(new TH2D("sdXYT2none45","Events with two tracks in arm 45, X0 vs Y0 (T2 empty);x;y",1000,-50.,50.,1000,-50.,50.));
      SDXvsY[1][0]=std::auto_ptr<TH2D>(new TH2D("sdXYT2same45","Events with two tracks in arm 45, X0 vs Y0 (T2 trk same side);x;y",1000,-50.,50.,1000,-50.,50.));
      SDXvsY[2][0]=std::auto_ptr<TH2D>(new TH2D("sdXYT2oppo45","Events with two tracks in arm 45, X0 vs Y0 (T2 trk opposite side);x;y",1000,-50.,50.,1000,-50.,50.));
      SDXvsY[3][0]=std::auto_ptr<TH2D>(new TH2D("sdXYT2both45","Events with two tracks in arm 45, X0 vs Y0 (T2 trk both sides);x;y",1000,-50.,50.,1000,-50.,50.));


      SDXvsY[0][1]=std::auto_ptr<TH2D>(new TH2D("sdXYT2none56","Events with two tracks in arm 56, X0 vs Y0 (T2 empty);x;y",1000,-50.,50.,1000,-50.,50.));
      SDXvsY[1][1]=std::auto_ptr<TH2D>(new TH2D("sdXYT2same56","Events with two tracks in arm 56, X0 vs Y0 (T2 trk same side);x;y",1000,-50.,50.,1000,-50.,50.));
      SDXvsY[2][1]=std::auto_ptr<TH2D>(new TH2D("sdXYT2oppo56","Events with two tracks in arm 56, X0 vs Y0 (T2 trk opposite side);x;y",1000,-50.,50.,1000,-50.,50.));
      SDXvsY[3][1]=std::auto_ptr<TH2D>(new TH2D("sdXYT2both56","Events with two tracks in arm 56, X0 vs Y0 (T2 trk both sides);x;y",1000,-50.,50.,1000,-50.,50.));


      */
    }
  //
  spectro45tp=std::auto_ptr<TH1D>(new TH1D("spectro45tp","spectrometer cut arm 45_tp",12000,-30.,30.));
  spectro45bt=std::auto_ptr<TH1D>(new TH1D("spectro45bt","spectrometer cut arm 45_bt",12000,-30.,30.));
  spectro56tp=std::auto_ptr<TH1D>(new TH1D("spectro56tp","spectrometer cut arm 56_tp",12000,-30.,30.));
  spectro56bt=std::auto_ptr<TH1D>(new TH1D("spectro56bt","spectrometer cut arm 56_bt",12000,-30.,30.));

  srcYHit=std::auto_ptr<TH1D>(new TH1D("srcYHit","source pot track y0-parameter (IsHit)",1000,-50.,50.));
  chkYHit=std::auto_ptr<TH1D>(new TH1D("chkYHit","checkPot track y0-parameter (IsHit)",1000,-50.,50.));
  
 
  srcUpos=std::auto_ptr<TH1D>(new TH1D("srcUpos","Local coordinate U of track fixed point in reference pot",201,-100.5,100.5));
  srcVpos=std::auto_ptr<TH1D>(new TH1D("srcVpos","Local coordinate V of track fixed point in reference pot",201,-100.5,100.5));
  trgUexPos=std::auto_ptr<TH1D>(new TH1D("trgUexPos","Local coordinate U of track extrapolated point in checkPot",2010,-100.5,100.5));
  trgVexPos=std::auto_ptr<TH1D>(new TH1D("trgVexPos","Local coordinate V of track extrapolated point in checkPot",2010,-100.5,100.5));

  trkTx=std::auto_ptr<TH1D>(new TH1D("trkTx","Fitted track angle in x",2000,-0.05,0.05));
  trkTy=std::auto_ptr<TH1D>(new TH1D("trkTy","Fitted track angle in y",2000,-0.05,0.05));


  elasticRefVsCheckY=std::auto_ptr<TProfile>(new TProfile("elasticRefVsCheckY","For events with TMC U=V=12 AND spectrometer cuts=true;refTrk Y0 (mm);checkPot Trk Y0 (mm)",150,-15,15));


  trkMisPointXY=std::auto_ptr<TH2D>(new TH2D("trkMisPointXY","For elastic events with checkTrk&refTrk, extrapol minus check x & y;x;y",450,-45,45,1800,-45.,45.));

  refEvtXY=std::auto_ptr<TH2D>(new TH2D("refEvtXY","For events with deltaSectors.GT.0, refTrk x0 vs y0",900,-45,45,900,-45.,45.));
  chkEvtXY=std::auto_ptr<TH2D>(new TH2D("chkEvtXY","For events with deltaSectors.GT.0 .AND. checkPot triggered, refTrk x0 vs y0",900,-45,45,900,-45.,45.));
  UvsTMCu=std::auto_ptr<TH2D>(new TH2D("UvsTMCu","First local coord in refPot plane 5 vs TMC bit on in U",200,-20.,20.,17,-0.5,16.5));
  UvsTMCv=std::auto_ptr<TH2D>(new TH2D("UvsTMCv","First local coord in refPot plane 5 vs TMC bit on in V",200,-20.,20.,17,-0.5,16.5));
  VvsTMCu=std::auto_ptr<TH2D>(new TH2D("VvsTMCu","Second local coord in refPot plane 5 vs TMC bit on in U",200,-20.,20.,17,-0.5,16.5));
  VvsTMCv=std::auto_ptr<TH2D>(new TH2D("VvsTMCv","Second local coord in refPot plane 5 vs TMC bit on in V",200,-20.,20.,17,-0.5,16.5));
  refTSvsExtrap=std::auto_ptr<TH2D>(new TH2D("refTSvsExtrap","refTrack local coord in plane 5: TrigSector vs same extrapol to checkPot pl5",17,-0.5,16.5,17,-0.5,16.5));
  refTSvsExtrapU=std::auto_ptr<TH2D>(new TH2D("refTSvsExtrapU","refTrack local coord in plane 4: TrigSector vs same extrapol to checkPot pl4",17,-0.5,16.5,17,-0.5,16.5));
  refTSvsCheckU=std::auto_ptr<TH2D>(new TH2D("refTSvsCheckU","refTrack local coord in plane 4: TrigSector vs same for checkPot track plane 4",17,-0.5,16.5,17,-0.5,16.5));
  refTSvsCheckV=std::auto_ptr<TH2D>(new TH2D("refTSvsCheckV","refTrack local coord in plane 5: TrigSector vs same for checkPot track plane 5",17,-0.5,16.5,17,-0.5,16.5));


 spectroElastUvsV=std::auto_ptr<TH2D>(new TH2D("spectroElastUvsV","RefPot TMC singleU+singleV: TrigSector distribution for spectrometer cutted all-events",22,-0.5,21.5,22,-0.5,21.5));

  spectroElastRefX=std::auto_ptr<TH1D>(new TH1D("spectroElastRefX","RefTrack x0() distribution for spectrometer-elastic evts",800,-40,40));
  pullX=std::auto_ptr<TH1D>(new TH1D("pullX","Pull between near and far trk x",800,-100,100));
  pullY=std::auto_ptr<TH1D>(new TH1D("pullY","Pull between near and far trk y",800,-100,100));
  pullYElast=std::auto_ptr<TH1D>(new TH1D("pullYElast","Pull between near and far trk y (for elastic evts)",800,-100,100));
  pullXpY=std::auto_ptr<TH2D>(new TH2D("pullXpY","Pull between near and far trk x vs pull in y",800,-100,100,2000,-100.,100.));
  pullXXc=std::auto_ptr<TH2D>(new TH2D("pullXXc","Pull between near and far trk x vs x_check",800,-100,100,600,-100.,100.));
  pullYYc=std::auto_ptr<TH2D>(new TH2D("pullYYc","Pull between near and far trk y vs y_check (cut in x)",800,-100,100,900,-45.,45.));
  pullYVsSpec=std::auto_ptr<TH2D>(new TH2D("pullYVsSpec","Pull between near and far trk y vs spectrometer cut value (elastic events)",800,-50,50,500,-5.,5.));
  //  pullYvsXY=std::auto_ptr<TH3D>(new TH3D("pullYvsXY","Pull between near and far trk y vs x_ch AND y_check",800,-100,100,600,-45.,45.,200,-45.,45.));
  refUVsCheckU=std::auto_ptr<TH2D>(new TH2D("refUVsCheckU","U-coord in refPot plane 5 vs refTrk check-U in same plane",800,-100,100,800,-100.,100.));
  refVVsCheckV=std::auto_ptr<TH2D>(new TH2D("refVVsCheckV","V-coord in refPot plane 5 vs refTrk check-V in same plane",800,-100,100,800,-100.,100.));

  if (readClusters)
    {
      numHitsInCheckPointing=std::auto_ptr<TH1D>(new TH1D("numHitsInCheckPointing","Num hits in checkpot for evts w pointing refTrks and not chkTMC U.OR.V on",600,0,600));

      prfNClu=std::auto_ptr<TProfile>(new TProfile("prfNClu","Num of clusters",128,-0.5,127.5));
      prfSClu=std::auto_ptr<TProfile>(new TProfile("prfSClu","Cluster size",128,-0.5,127.5));
      prfOccup=std::auto_ptr<TProfile>(new TProfile("prfOccup","Strip occupancy (plane not empty)",1280,-0.5,1279.5));
      prfMultipl=std::auto_ptr<TProfile>(new TProfile("prfMultipl","Number of tracks versus number of clusters in all pots",99,CBins));

      prfMul[0]=std::auto_ptr<TProfile>(new TProfile("prfMul020","Number of tracks versus number of clusters in pot 020",99,CBins));
      prfMul[1]=std::auto_ptr<TProfile>(new TProfile("prfMul021","Number of tracks versus number of clusters in pot 021",99,CBins));
      prfMul[2]=std::auto_ptr<TProfile>(new TProfile("prfMul022","Number of tracks versus number of clusters in pot 022",99,CBins));
      prfMul[3]=std::auto_ptr<TProfile>(new TProfile("prfMul023","Number of tracks versus number of clusters in pot 023",99,CBins));
      prfMul[4]=std::auto_ptr<TProfile>(new TProfile("prfMul024","Number of tracks versus number of clusters in pot 024",99,CBins));
      prfMul[5]=std::auto_ptr<TProfile>(new TProfile("prfMul025","Number of tracks versus number of clusters in pot 025",99,CBins));
      prfMul[6]=std::auto_ptr<TProfile>(new TProfile("prfMul120","Number of tracks versus number of clusters in pot 120",99,CBins));
      prfMul[7]=std::auto_ptr<TProfile>(new TProfile("prfMul121","Number of tracks versus number of clusters in pot 121",99,CBins));
      prfMul[8]=std::auto_ptr<TProfile>(new TProfile("prfMul122","Number of tracks versus number of clusters in pot 122",99,CBins));
      prfMul[9]=std::auto_ptr<TProfile>(new TProfile("prfMul123","Number of tracks versus number of clusters in pot 123",99,CBins));
      prfMul[10]=std::auto_ptr<TProfile>(new TProfile("prfMul124","Number of tracks versus number of clusters in pot 124",99,CBins));
      prfMul[11]=std::auto_ptr<TProfile>(new TProfile("prfMul125","Number of tracks versus number of clusters in pot 125",99,CBins));
      
      for (int t=0;t<12;t++)
	{
	  prfMul[t]->SetDirectory(0); 
	  prfMul[t]->SetXTitle("Number of clusters");
	  prfMul[t]->SetYTitle("number of tracks");
      
	}
    

      prfNClu->SetDirectory(0); 

      prfSClu->SetDirectory(0); 

      prfSClu->SetXTitle("RP id");
      prfNClu->SetXTitle("RP id");

      prfOccup->SetDirectory(0); 
      
      prfOccup->SetXTitle("RP id x 10 plus plane nr");
      
      prfMultipl->SetDirectory(0); 
      
      prfMultipl->SetXTitle("Number of clusters");
      prfMultipl->SetYTitle("number of tracks");
    }
}


RPDataReduction::~RPDataReduction()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
RPDataReduction::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  const unsigned int potsAll[8]={20,21,24,25,120,121,124,125};
  const unsigned int potsAllHz[12]={20,21,22,23,24,25,120,121,122,123,124,125};
  //  const unsigned int potsD2[4]={21,25,120,124};

#ifdef FOdebug    
  std::cout<<"DEBUG: Begin analyze event"<<std::endl;
#endif

  TotemRPGeometry Totem_RP_geometry_;
 

  edm::ESHandle<TotemRPGeometry> Totem_RP_geometry;
  iSetup.get<RealGeometryRecord>().get(Totem_RP_geometry);
  if(!Totem_RP_geometry.isValid()) 
    {
      throw cms::Exception("RPDataReduction") << "edm::ESHandle<TotemRPGeometry> is invalid, RP finding impossible!";
    }
 
  if (refPot==checkPot)
    throw cms::Exception("RPDataReduction") << "Trying to extrapolate track from pot "<<refPot<<" to itself!";

  if ((refPot!=(checkPot+4))&&(refPot!=(checkPot-4)))
    throw cms::Exception("RPDataReduction") << "Extrapolating track from refPot "<<refPot<<" to a non-collinear pot! (checkPot:"<<checkPot<<")";

 
  for(unsigned int arm = 0 ; arm < 2 ; ++arm)                     // Loop over all of the arms.
    {
      std::set<unsigned int> detects = Totem_RP_geometry->DetsInRP(arm ? refPot : checkPot);
      //      for (std::set<unsigned int>::const_iterator st = detects.begin(); st != detects.end(); ++st)
      //	std::cout<<(arm ? " refPot: " : " checkPot: ")<<" detector "<<(*st);
      //      std::cout<<std::endl;
	// Loop over all the stations of the actual arm.
    }

  // Step B: Get Inputs
  edm::Handle< edm::DetSetVector<RPDigCluster> >  input;
  edm::Handle<RPMulFittedTrackCollection>  inputTrk;
  edm::Handle<RPFittedTrackCollection>  input1Trk;
  edm::Handle<RPRecognizedPatternsCollection>  inputLines;
  edm::Handle<RPTrackCandidateCollection>  inputPossTrk;
  edm::Handle<RPMulTrackCandidateCollection>  inputPossMTrk;
  edm::Handle<RPReconstructedProtonCollection> inputProtons;
  edm::Handle<T1T2TrackCollection>  inputT2Trk;
  edm::Handle<T1T2TrackCollection>  inputT1Trk;
  //  edm::Handle<T2PadClusterCollection> inputT2Pad;
  edm::Handle<CCVect>  inputTMC;
 
  edm::Handle< Totem::RawEvent> inputRaw;

  unsigned int trigStatus=0;

#ifdef FOdebug    
  std::cout<<"DEBUG: Declare input classes"<<std::endl;
#endif

 
  if (readRecoProt)
    iEvent.getByLabel(rPReconstructedProtonCollectionLabel, inputProtons);
    
   
  int bunchNumSD=-1;

  int T2TrigLN=-1;

  if (readLoNeg)
    {
      iEvent.getByLabel(rawEventLabel, inputRaw);


      dataT.type = inputRaw->triggerData.type;
      dataT.event_num = inputRaw->triggerData.event_num;
      dataT.bunch_num = inputRaw->triggerData.bunch_num;
      bunchNumSD=dataT.bunch_num;

      if (falseBunch==bunchNumSD)
	throw cms::Exception("RPDataReduction") << "This event has a real bunch number "<<bunchNumSD<<" equal to the chosen value for the parameter whichBunchAll, CHANGE THE LATTER! (false bunch number used for bunch-nondiscriminating sums & counts.)";



      lonegBunch[0]->Fill(dataT.bunch_num );


      if (!refEvtsBunches.count(bunchNumSD))
	{
	  refEvtsBunches[bunchNumSD]=1;
	}
      else
	{
	  refEvtsBunches[bunchNumSD]++;
	}


      dataT.src_id = inputRaw->triggerData.src_id;
      dataT.orbit_num = inputRaw->triggerData.orbit_num;
      dataT.revision_num = inputRaw->triggerData.revision_num;
      dataT.run_num = inputRaw->triggerData.run_num;
      dataT.trigger_num = inputRaw->triggerData.trigger_num;
      dataT.inhibited_triggers_num = inputRaw->triggerData.inhibited_triggers_num;
      dataT.input_status_bits = inputRaw->triggerData.input_status_bits;
      trigStatus=dataT.input_status_bits;
      if (whichVerb>15)
	std::cout<<"dataT.type:"<<(unsigned int) dataT.type<<"  dataT.event_num:"<<dataT.event_num<<" dataT.input_status_bits:"<<dataT.input_status_bits<<std::endl;
      if (trigStatus%2)
	{
	  trig220v[falseBunch]++;

	  if (!trig220v.count(bunchNumSD))
	    {
	      trig220v[bunchNumSD]=1;
	    }
	  else
	    {
	      trig220v[bunchNumSD]++;
	    }


	  lonegBunch[1]->Fill(dataT.bunch_num );
	}
      if ((trigStatus/2)%2)
	{
	  trig220h[falseBunch]++;


	  if (!trig220h.count(bunchNumSD))
	    {
	      trig220h[bunchNumSD]=1;
	    }
	  else
	    {
	      trig220h[bunchNumSD]++;
	    }

	  lonegBunch[2]->Fill(dataT.bunch_num );
	}
      if ((trigStatus/4)%2)
	{
	  trig220cross[falseBunch]++;

	  if (!trig220cross.count(bunchNumSD))
	    {
	      trig220cross[bunchNumSD]=1;
	    }
	  else
	    {
	      trig220cross[bunchNumSD]++;
	    }

	}


      T2TrigLN=(trigStatus/64)%2;


      if ((trigStatus/64)%2)
	{
	  trigT2[falseBunch]++;

	  if (!trigT2.count(bunchNumSD))
	    {
	      trigT2[bunchNumSD]=1;
	    }
	  else
	    {
	      trigT2[bunchNumSD]++;
	    }


	  lonegBunch[4]->Fill(dataT.bunch_num );
      
	}
      if ((trigStatus/256)%2)
	{
	  trigT1[falseBunch]++;

	  if (!trigT1.count(bunchNumSD))
	    {
	      trigT1[bunchNumSD]=1;
	    }
	  else
	    {
	      trigT1[bunchNumSD]++;
	    }


	  lonegBunch[3]->Fill(dataT.bunch_num );

	}
      if ((trigStatus/512)%2)
	{
	  trigBx[falseBunch]++;


	  if (!trigBx.count(bunchNumSD))
	    {
	      trigBx[bunchNumSD]=1;
	    }
	  else
	    {
	      trigBx[bunchNumSD]++;
	    }


	  lonegBunch[5]->Fill(dataT.bunch_num );

	}
    }

#ifdef FOdebug    
  std::cout<<"DEBUG: Have read LoNeg class"<<std::endl;
#endif


  //e.getByLabel(digiProducer_,digiLabel_,input);  //FIXME: fix this label
  if (readClusters)
    iEvent.getByLabel(detSetVectorRPDigClusterLabel, input);  //FIXME: fix this label
  if (readMultiTrk)
    iEvent.getByLabel(rPMulFittedTrackCollectionLabel, inputTrk);  //FIXME: fix this label
  iEvent.getByLabel(rPFittedTrackCollectionLabel, input1Trk);  //FIXME: fix this label


 
  //  iEvent.getByType(inputTMC);  //FIXME: fix this label
  //  iEvent.getByLabel("RPDataCCProducer","RPCCRawBits",inputTMC);  //FIXME: fix this label
  iEvent.getByLabel(tmcModule.c_str(),tmcProd.c_str(),inputTMC);  //FIXME: fix this label

  if (readT2)
    {
      iEvent.getByLabel(t2trModule.c_str(),t2trProd.c_str(),inputT2Trk);  //FIXME: fix this label
      //      iEvent.getByLabel(t2padModule.c_str(),t2padProd.c_str(),inputT2Pad);  //FIXME: fix this label
    }



  if (readT1)
    iEvent.getByLabel(t1trModule.c_str(),t1trProd.c_str(),inputT1Trk);  //FIXME: fix this label
    


  iEvent.getByLabel(rPTrackCandidateCollectionLabel, inputPossTrk);  //FIXME: fix this label
  if (readMultiTrk)
    iEvent.getByLabel(rPMulTrackCandidateCollectionLabel, inputPossMTrk);  //FIXME: fix this label
  iEvent.getByLabel(rPRecognizedPatternsCollectionLabel, inputLines);  //FIXME: fix this label
 
  totalEvts++;


  //loop on all detset inside the input collection
  edm::DetSetVector<RPDigCluster>::const_iterator DSViter;
  if (readClusters)
    DSViter=(*input).begin();
  //  int aa=0;


  int OneEvtTrkSum=0;
  int OneEvtCluSum=0;

  int OnePotTrkSum[128]={0};
  int OnePotCluSum[128]={0};

  std::map<RPId, std::vector<RPFittedTrack> >::const_iterator TrkIter;
  std::map<RPId, std::vector<RPTrackCandidate> >::const_iterator PossMTrkIter;

  std::vector<RPReconstructedProton> ::const_iterator protonIter;

  if (readRecoProt)
    protonIter=(*inputProtons).begin();
    
  //  std::map<T2DetId, std::vector<T2Cluster> >::const_iterator t2PadIter;
  std::vector<T1T2Track>::const_iterator t2TrkIter; 
  std::vector<T1T2Track>::const_iterator t1TrkIter; 
  if (readT2)
    {
      t2TrkIter=(*inputT2Trk).begin();
      //     t2PadIter=(*inputT2Pad).begin();
    }

  if (readT1)
    t1TrkIter=(*inputT1Trk).begin();

  if (readMultiTrk)
    TrkIter=(*inputTrk).begin();
  std::map<RPId, RPFittedTrack>::const_iterator OneTrkIter=(*input1Trk).begin();
  std::map<RPId, RPTrackCandidate>::const_iterator PossTrkIter=(*inputPossTrk).begin();
  if (readMultiTrk)
    PossMTrkIter=(*inputPossMTrk).begin();


  std::map<RPId, RPRecognizedPatterns>::const_iterator LineIter=(*inputLines).begin();
  
  std::vector<RPCCBits>::const_iterator TMCIter=(*inputTMC).begin();
 
  bool multiT[128];
  for (int y=0;y<128;y++)
    multiT[y]=false;
  bool noTr[128];
  for (int y=0;y<128;y++)
    noTr[y]=false;
  bool siTr[128];
  for (int y=0;y<128;y++)
    siTr[y]=false;

  bool potHasMTrk1[4]={false};

  bool potHasMTrkCand[8]={false};
  bool potHasSTrkCand[8]={false};
  bool AllPotMTrk[8]={false};
  bool AllPotSTrk[8]={false};

  int HzAllPotClu[12]={0};
  int HzAllPotTmcU[12]={0};
  int HzAllPotTmcV[12]={0};


  bool HzAllPot1MTrk[12]={false};
  bool HzAllPotMoreMTrk[12]={false};
  bool HzAllPotSTrk[12]={false};
  bool HzLineInUOnly[12]={false};
  bool HzLineInVOnly[12]={false};

  bool lineInUOnly[4]={false};
  bool lineInVOnly[4]={false};

  bool T2minusTrk=false;
  bool T2plusTrk=false;
  //  bool T2minusPad=false;
  //  bool T2plusPad=false;
  int T2PlusMulti=0;
  int T2MinusMulti=0;

  bool T1minusTrk=false;
  bool T1plusTrk=false;
  int T1PlusMulti=0;
  int T1MinusMulti=0;

  double T2maxEta=5.3;
  double T2minEta=-5.3;
  double T2maxTrk=-20;
  double T2minTrk=20;
  double T1T2maxTrk=-20;
  double T1T2minTrk=20;

  double xiErr,xiErr2;

  DMapType Tvals;
  std::vector<double> ThetaXes,ThetaYs;
  Tvals.clear();
  ThetaXes.clear();
  ThetaYs.clear();

  int numProt=-1;


#ifdef FOdebug    
  std::cout<<"DEBUG: Next read proton class"<<std::endl;
#endif


  if (readRecoProt)
    {
      RPRecoProtMADXVariables findT;
      double tee,xiProt,thetaX,thetaY;
      for (;protonIter!=inputProtons->end();protonIter++)
	{
	  numProt=inputProtons->size();
	  xiProt=protonIter->Ksi();
	  xiErr2=protonIter->CovarianceMartixElement(nksi,nksi);
	  if (xiErr2>0)
	    xiErr=sqrt(xiErr2);
	  else
	    {
	      if (whichVerb>10)
		std::cout<<"xi auto-covariance =<0 : "<<xiErr2<<" , set equal to 1%"<<std::endl;
	      xiErr=-0.01;
	    }
	  thetaX=protonIter->Theta_x();
	  thetaY=protonIter->Theta_y();
	  ThetaXes.push_back(thetaX);
	  ThetaYs.push_back(thetaY);

	  recoProtXiVsN->Fill(numProt,xiProt);

	  recoProtXiVsErr->Fill(xiProt,xiErr);

	  findT=protonIter->GetMADXVariables();

	  tee=optObj.MADXCanonicalVariablesTot(findT);
	  Tvals.insert(pair<double,double>(tee,xiProt));

	  if (whichVerb>10)
	    std::cout<<"Proton t-value is: "<<tee<<" (GeV^2)"<<std::endl;

	}
    }

#ifdef FOdebug    
  std::cout<<"DEBUG: Have read proton class, next T1 tracks"<<std::endl;
#endif



  if (readT1&&readT2)
    {
      T1RecHitGlobal tempHit;
      float hitEta,avgHitEta;
      for (;t1TrkIter!=(*inputT1Trk).end();t1TrkIter++)
	{
	  avgHitEta=0;
	  float eta=(*t1TrkIter).Eta();
	  int hits=t1TrkIter->GetHitEntries();
	  for (int u=0;u<hits;u++)
	    {
	      tempHit=t1TrkIter->GetHitT1(u);
	      hitEta=tempHit.eta();
	      avgHitEta+=hitEta/hits;
	      T1Eta000Dispersion->Fill(hitEta-eta);
	    }

	  T1FitEtaVs000Eta->Fill(avgHitEta,eta);

	  if (avgHitEta>T1T2maxTrk)
	    T1T2maxTrk=avgHitEta;
	  if (avgHitEta<T1T2minTrk)
	    T1T2minTrk=avgHitEta;


	  if (avgHitEta>0)
	    {
	      T1plusTrk=true;
	      T1PlusMulti++;
	    }
	  if (avgHitEta<0)
	    {
	      T1minusTrk=true;
	      T1MinusMulti++;
	    }

	}
    }

#ifdef FOdebug    
  std::cout<<"DEBUG: Have read T1 class, next T2 tracks"<<std::endl;
#endif

  if (readT2)
    {
      T2Hit tempHit;
      float hitEta,r,z,c,avgHitEta;
      for (;t2TrkIter!=(*inputT2Trk).end();t2TrkIter++)
	{
	  avgHitEta=0;
	  float eta=(*t2TrkIter).Eta();
	  int hits=t2TrkIter->GetHitEntries();
	  for (int u=0;u<hits;u++)
	    {
	      tempHit=t2TrkIter->GetHitT2(u);
	      r=tempHit.GetHitR();
	      z=tempHit.GetHitZ();

	      //theta
	      if(z>0) c = atan(r/z);
	      else if(z<0) c = atan(r/z)+3.14159;
	      else {c = 3.14159;}
	      //pseudorapidity
	      hitEta = -log(tan(c/2.));

	      avgHitEta+=hitEta/hits;
	      T2Eta000Dispersion->Fill(hitEta-eta);
	    }
	  T2FitEtaVs000Eta->Fill(avgHitEta,eta);

	  if (avgHitEta>T2maxEta)
	    T2maxEta=avgHitEta;
	  if (avgHitEta<T2minEta)
	    T2minEta=avgHitEta;

	  if (avgHitEta>T2maxTrk)
	    T2maxTrk=avgHitEta;
	  if (avgHitEta<T2minTrk)
	    T2minTrk=avgHitEta;

	  if (avgHitEta>T1T2maxTrk)
	    T1T2maxTrk=avgHitEta;
	  if (avgHitEta<T1T2minTrk)
	    T1T2minTrk=avgHitEta;

	  if (avgHitEta>0)
	    {
	      T2plusTrk=true;
	      T2PlusMulti++;
	    }
	  if (avgHitEta<0)
	    {
	      T2minusTrk=true;
	      T2MinusMulti++;
	    }
	}
      //      for (;t2PadIter!=(*inputT2Pad).end();t2PadIter++)
      //	{
      //	}
    }

  for (; LineIter!=(*inputLines).end();LineIter++)
   {
     for (int u=0;u<4;u++)
       {
	 bool lineU=false;
	 bool lineV=false;
	 if (((*LineIter).first==potsD2[diagonal][u])&&((*LineIter).second.uLines.size()))
	   lineU=true;
	 if (((*LineIter).first==potsD2[diagonal][u])&&((*LineIter).second.vLines.size()))
	   lineV=true;
	 if (lineU!=lineV)
	   {
	     (lineU ?lineInUOnly[u] : lineInVOnly[u])=true;
	   }
	 
       }


     for (int u=0;u<12;u++)
       {
	 bool lineU=false;
	 bool lineV=false;
	 if (((*LineIter).first==potsAllHz[u])&&((*LineIter).second.uLines.size()))
	   lineU=true;
	 if (((*LineIter).first==potsAllHz[u])&&((*LineIter).second.vLines.size()))
	   lineV=true;
	 if (lineU!=lineV)
	   {
	     float onePot=(lineU ? 21. : 22.);
	     (lineU ? HzLineInUOnly[u] : HzLineInVOnly[u]) = true;
	     outputPerPot[u]->Fill(onePot);
	   }
	 
       }



   }

 

  for (; PossTrkIter!=(*inputPossTrk).end();PossTrkIter++)
    {

//     for (int u=0;u<4;u++)
//	{
//	  if (((*PossTrkIter).first==potsD2[diagonal][u])&&((*PossTrkIter).second.Fittable()))
//	    {
	      //	      RPFittedTrack testU=((*TrkIter).second).at(0);
	      //	      if (testU.IsValid())
	      //		potHasMTrk1[u]=true;
//	    }
//	}
      for (int u=0;u<8;u++)
	{
	  if (((*PossTrkIter).first==potsAll[u])&&((*PossTrkIter).second.Fittable()))
	    {
	      potHasSTrkCand[u]=true;
	      //	      RPFittedTrack testU=((*TrkIter).second).at(0);
	      //	      if (testU.IsValid())
	      //		potHasMTrk1[u]=true;
	    }
	}


      for (int u=0;u<12;u++)
	{
	  if (((*PossTrkIter).first==potsAllHz[u])&&((*PossTrkIter).second.Fittable()))
	    {
	     float onePot=23;
	     outputPerPot[u]->Fill(onePot);
	    }
	}





    }

#ifdef FOdebug    
  std::cout<<"DEBUG: Have read trackLines & possible track candidates class, multitracks next"<<std::endl;
#endif




  bool HzPotM0S[12]={false};

  if (readMultiTrk)
    {
      for (; PossMTrkIter!=(*inputPossMTrk).end();PossMTrkIter++)
	{
	  
	  for (int u=0;u<8;u++)
	    {
	      if (((*PossMTrkIter).first==potsAll[u])&&((*PossMTrkIter).second.size()))
		{
		  potHasMTrkCand[u]=true;
		}
	    }
	  
	  for (int u=0;u<12;u++)
	    {
	      if (((*PossMTrkIter).first==potsAllHz[u])&&((*PossMTrkIter).second.size()))
		{
		  float onePot=24;
		  outputPerPot[u]->Fill(onePot);
		}
	    }
	}
    

      for (; TrkIter!=(*inputTrk).end();TrkIter++)
	{
	  
	  for (int u=0;u<4;u++)
	    {
	      if (((*TrkIter).first==potsD2[diagonal][u])&&((*TrkIter).second.size()))
		{
		  RPFittedTrack testU=((*TrkIter).second).at(0);
		  if (testU.IsValid())
		    potHasMTrk1[u]=true;
		}
	    }
	  
	  for (int u=0;u<8;u++)
	    {
	      if (((*TrkIter).first==potsAll[u])&&((*TrkIter).second.size()))
		AllPotMTrk[u]=true;
	    }
	  
	  
	  bool potHasSTrk=false;
	  for (OneTrkIter=(*input1Trk).begin(); OneTrkIter!=(*input1Trk).end();OneTrkIter++)
	    {
	      if ((*OneTrkIter).first==(*TrkIter).first)
		{
		  potHasSTrk=true;
		}
	    }


	  for (int u=0;u<12;u++)
	    {
	      if ((*TrkIter).first==potsAllHz[u])
		{
		  int numM=(*TrkIter).second.size();
		  if (numM&&(!potHasSTrk))
		    {
		      HzPotM0S[u]=true;
		    }
		  if (numM==1)
		    {
		      HzAllPot1MTrk[u]=true;
		      outputPerPot[u]->Fill(25.);
		      outputPerPot[u]->Fill(potHasSTrk ? 27. : 18.);
		    }
		  if (numM>1)
		    {
		      HzAllPotMoreMTrk[u]=true;
		      outputPerPot[u]->Fill(26.);
		      if (!potHasSTrk)
			outputPerPot[u]->Fill(19.);
		    }
		}
	    }


	  if ((*TrkIter).second.empty())
	    {
	      noTr[(*TrkIter).first]=true;
	      std::cout<<"pot entry with no tracks!!!"<<std::endl;
	    }
	  //      if (!(totalEvts%300))
	  if (whichVerb>10)
	    {
	      if (potHasSTrk) 
		{
		  std::cout<<"Mul-trk & single trk both found in pot: "<<(*TrkIter).first<<", mulReco found "<<TrkIter->second.size()<<" tracks";
		  std::cout<<std::endl;
		}
	      else
		{
		  std::cout<<"Mul-trk but not single trk found in pot: "<<(*TrkIter).first<<", mulReco found "<<TrkIter->second.size()<<" tracks";
		  std::cout<<std::endl;
		}
	    }
	  
	  if ((*TrkIter).second.size()>1)
	    multiT[(*TrkIter).first]=true;
	  if ((*TrkIter).second.size()==1)
	    siTr[(*TrkIter).first]=true;
	  OneEvtTrkSum+=(*TrkIter).second.size();
	  OnePotTrkSum[(*TrkIter).first]+=(*TrkIter).second.size();
	}
      for (int u=0;u<128;u++)
	{
	  if (multiT[u])
	    multiTrk[u]++;
	  if (siTr[u])
	    singleTOnly[u]++;
	  if ((!multiT[u])&&(!siTr[u])&& ((u%100)<26)&&((u%100)>19))
	    {
	      noTrk[u]++;
	      noTr[u]=true;
	    }
	}
      
    }

  bool rp45Trk=false;
  bool rp56Trk=false;
  bool rp45Tmc=false;
  bool rp56Tmc=false;

  

  bool trackPoints2Check=false;
  bool checkTrack=false;
  bool refTrack=false;
  double xRefEx=0;
  double yRefEx=0;
  double xRefMid=0;
  //  double checkY0=0;
  //  double refY0=0;
  double yRefMid=0;
  double xCheck=0;
  double yCheck=0;
  double sigmaX=0;
  double sigmaY=0;
  //  bool trackWithin25mm=false;
  double cutoff=0;
  RPTopology foTopo;
  double fracEvt=0;
  double uRef=0;
  double vRef=0;
  double uCheck=0;
  double uReallyRef=0;
  double vReallyRef=0;
  double uRefU=0;
  double uCheckU=0;
  double uReallyRefU=0;

  double x0rp[12]={0};
  double y0rp[12]={0};
  double z0rp[12]={0};

  //diag 1
  const double cn024=-0.701656557291165384;
  const double cn020=0.712515316053142556;
  const double delta45tp=-0.0518657;
  const double lim45tp=3*0.0138554;
  const double cn125=-0.699685961249467803;
  const double cn121=0.714450527069865227;
  const double delta56bt=-0.000307259;
  const double lim56bt=3*0.0128446;

  //diag2
  const double cn124=-0.698033146310773733;
  const double cn120=0.716065448580981090;
  const double delta56tp=0.0161861;
  const double lim56tp=3*0.0142454;
  const double cn025=-0.702238331590702591;
  const double cn021=0.711941939798960255;
  const double delta45bt=-0.0563749;
  const double lim45bt=3*0.0132129;

  //  double sigmas=1;

#ifdef FOdebug    
  std::cout<<"DEBUG: Have read Multitrack class, singleTrack next"<<std::endl;
#endif


  bool multiSTrk[12]={false};
  bool potHasSTrk1[4]={false};

  for (OneTrkIter=(*input1Trk).begin(); OneTrkIter!=(*input1Trk).end();OneTrkIter++)
    {


     for (int u=0;u<8;u++)
	{
	  if ((*OneTrkIter).first==potsAll[u]) //&&((*OneTrkIter).second.size()))
	    AllPotSTrk[u]=true;
	}


      //      fracEvt=0;
      if ((*OneTrkIter).first==oppositePot)
        {
	  RPFittedTrack testO=(*OneTrkIter).second;
//	  if (testO.IsValid())
        }

    for (int u=0;u<12;u++)
	{
	  if ((*OneTrkIter).first==potsAllHz[u]) //&&((*OneTrkIter).second.size()))
	    {
	      RPFittedTrack testO=(*OneTrkIter).second;
	      if (testO.IsValid())
		{
		  if (multiSTrk[u])
		    throw cms::Exception("RPDataReduction") << "More than one singleTrk in pot "<<potsAllHz[u];

		  HzAllPotSTrk[u]=true;
		  x0rp[u]=testO.X0();
		  y0rp[u]=testO.Y0();
		  z0rp[u]=testO.Z0();
		  //		  if (potsAllHz[u]<100)
		  //		    {
		  //		      rp45Trk=true;
		  //		    }
		  //		  else
		  //		    {
		  //		      rp56Trk=true;
		  //		    }
		  outputPerPot[u]->Fill(16.);
		  multiSTrk[u]=true;
		}
	    }
	}

      for (int u=0;u<4;u++)
	{
	  if ((*OneTrkIter).first==potsD2[diagonal][u])
	    {
	      RPFittedTrack testU=(*OneTrkIter).second;
	      if (testU.IsValid())
		potHasSTrk1[u]=true;
	    }
	}

      if ((*OneTrkIter).first==checkPot)
	{
	  RPFittedTrack testR=(*OneTrkIter).second;
	  if (testR.IsValid())
	    {
	      checkTrack=true;
	      //	      TotRPDetId check = TotRPDetId(checkPot/100, (checkPot%100)/10, checkPot%10, 0);
	      //	      double detPos=Totem_RP_geometry->GetDetTranslation(check).z();
	      //	      TVector2 checkTrk=testR.GetTrackPoint(detPos);
	      //	      xCheck=checkTrk.X();
	      //	      yCheck=checkTrk.Y();
	      xCheck=testR.X0();
	      yCheck=testR.Y0();
	      chkY->Fill(testR.Y0());

	      TotRPDetId check2 = TotRPDetId(checkPot/100, (checkPot%100)/10, checkPot%10, 5);
	      TotRPDetId check2U = TotRPDetId(checkPot/100, (checkPot%100)/10, checkPot%10, 4);
	      CLHEP::Hep3Vector fixLocalPos3,fixPos3,fixLocalPos3U;
	      fixPos3.setX(testR.X0());
	      fixPos3.setY(testR.Y0());
	      fixPos3.setZ(testR.Z0());
	      fixLocalPos3=Totem_RP_geometry->GlobalToLocal(check2,fixPos3);
	      fixLocalPos3U=Totem_RP_geometry->GlobalToLocal(check2U,fixPos3);
	      uCheck=fixLocalPos3.x();
	      fixLocalPos3.y();
	      uCheckU=fixLocalPos3U.x();
	      fixLocalPos3U.y();

	      if (foTopo.IsHit(fixLocalPos3.x(),fixLocalPos3.y()))
		chkYHit->Fill(testR.Y0());

	    }
	}


      if ((*OneTrkIter).first==refPot)
	{
	  
	  if (whichVerb>5)
	    {
	      std::cout<<"Single-trk found in reference pot: "<<(*OneTrkIter).first;//<<", mulReco found "<<OneTrkIter->second.size()<<" tracks";
	      std::cout<<std::endl;
	    }
	  RPFittedTrack test=(*OneTrkIter).second;
	  if (test.IsValid())
	    {
	      refTrack=true;
	      srcY->Fill(test.Y0());
	      double detPos,errPX,errPY;
	      CLHEP::Hep3Vector extrapPos,localPos,fixPos,fixLocalPos,fixLocalPos2,fixLocalPosU4;
	      //	      std::vector<CLHEP::Hep3Vector> localPos;
	      //	      localPos.clear();

	      fracEvt=0;
	      //if (test.GetTx()>(1.00*32*0.066/5000.))
		//	continue; //protons are parallel to the beam in one dimension only? (x angle up to 300 microrad, y up to 15 microrad -> both less or equal than one trigger sector)
	      //	      if (test.GetTy()>(1.00*32*0.066/5000.))
	      //	continue;
	      if (whichVerb>10)
		std::cout<<"Parallel track found in reference pot: "<<std::endl;
	      TotRPDetId refIDu = TotRPDetId(refPot/100, (refPot%100)/10, refPot%10, 0);
	      TotRPDetId refIDv = TotRPDetId(refPot/100, (refPot%100)/10, refPot%10, 5);
	      TotRPDetId refIDu4 = TotRPDetId(refPot/100, (refPot%100)/10, refPot%10, 4);

	      fixPos.setX(test.X0());
	      fixPos.setY(test.Y0());
	      fixPos.setZ(test.Z0());
	      fixLocalPos=Totem_RP_geometry->GlobalToLocal(refIDu,fixPos);
	      fixLocalPos2=Totem_RP_geometry->GlobalToLocal(refIDv,fixPos);
	      fixLocalPosU4=Totem_RP_geometry->GlobalToLocal(refIDu4,fixPos);
	      uReallyRef=fixLocalPos2.x();
	      vReallyRef=fixLocalPos2.y();
	      uReallyRefU=fixLocalPosU4.x();
	      fixLocalPosU4.y();


	      if (foTopo.IsHit(fixLocalPos2.x(),fixLocalPos2.y()))
		srcYHit->Fill(test.Y0());

	      //if (whichVerb>10)
		
		  //		  std::cout<<"x0error: "<<test.X0Sigma()<<" y0error: "<<test.Y0Sigma()<<std::endl;
		  //		  std::cout<<"v_lower edge (0): "<<foTopo.GetHitPositionInReadoutDirection(0.)<<"v_upper edge (511)"<<foTopo.GetHitPositionInReadoutDirection(511.)<<std::endl;
		

	      std::set<unsigned int> detects = Totem_RP_geometry->DetsInRP(checkPot);
	      for (std::set<unsigned int>::const_iterator st = detects.begin(); st != detects.end(); ++st)
		{
		  TotRPDetId near = TotRPDetId(checkPot/100, (checkPot%100)/10, checkPot%10, ((*st)%10));
		  detPos=Totem_RP_geometry->GetDetTranslation(near).z();
		  test.Z0();
		  //		  std::cout<<" tot rp id is="<<near.DetectorDecId();
		  test.GetTxSigma();
		  test.GetTySigma();
		  errPX=fabs(errX);
		  errPY=fabs(errY);
		  //		  TVector2 extrapol=test.GetTrackPoint(detPos);

		  double xSlope,ySlope;
		  if (refPot<checkPot)
		    {
		      xSlope=0.87; 
		      //		      ySlope=1.09; //general
		      ySlope=1.02; //elastic
		      if (refPot==21)
			ySlope=1.01; //elastic
			
		    }
		  else
		    {
		      xSlope=1./0.87;
		      //ySlope=1./1.09; //general
		      ySlope=1./1.02; //elastic
		      if (refPot==25)
			ySlope=1./1.01; //elastic
		    }
		  TVector2 extrapol(xSlope*test.X0(),ySlope*test.Y0());//=test.GetTrackPoint(detPos);
		  //this is correct when going from nearPot=refPot to farPot=checkPot


		  extrapPos.setZ(detPos);
		 
		  if (((*st)%10)==5)
		    {
		      xRefMid=test.X0();
		      yRefMid=test.Y0();
		      xRefEx=extrapol.X();
		      yRefEx=extrapol.Y();
		      test.GetTxSigma();
		      test.GetTySigma();
		      sigmaX=fabs(errX/sigmas);
		      sigmaY=fabs(errY/sigmas);
		    }
		  

		  //no sigma checks
		  //add global to local transformation

		  int fracHit=0;
		  bool hitty=false;
		  double offset56=0;
		  double offset45tp=0;
		  double offset45bt=0;
		  if (refPot==120||refPot==121)
		    offset56=fabs(elCut56);
		  //		  if (refPot==20)
		  // 15.6. remember to change this back soonest
		//		  if (refPot==20||refPot==24)
		//  26.10 changed back
		  if (refPot==20) 
		    offset45tp=fabs(elCut45t);
		  if (refPot==21)
		    offset45bt=fabs(elCut45b);

		  // 5-6 tp & bottom near start at abs(Y0)~7.35 ... 7.40, far pots start at ~8.4mm +- 50mu
		  // 4-5 near tp starts at 7.82 +- 50mu, far tp at 7.92mm +- 80mu
		  // 4-5 near bt -7.78 +- 50mu, far at -7.99mm +- 50mu
		  // (last digit is in the above is a guess based on the bin contents for the last & next-to-last bin - bin size 100mu)
		  //		  if (refPot==121)
		  //		    offset56=-elCut56;

		  if (sigmas>0.01)
		    {
		      for (int u=0;u<9;u++)
			{
			  extrapPos.setX(extrapol.X()+errPX*(u/4.-1.));
			  for (int v=0;v<9;v++)
			    {
			      double sumOff=offset56+offset45bt+offset45tp;
			      //			      if (fabs(test.Y0())>offset56)
			      if (fabs(test.Y0())>sumOff)
				{
				  //			      double offset56=((refPot/100) ? 
				  extrapPos.setY(extrapol.Y()+errPY*(v/4.-1.));
				  localPos=Totem_RP_geometry->GlobalToLocal(near,extrapPos);
				  hitty=foTopo.IsHit(localPos.x(),localPos.y());
				  if ((v==4)&&(u==4)&&(((*st)%10)==5))
				    {
				      uRef=localPos.x(); 
				      vRef=localPos.y();
				    }

				  if ((v==4)&&(u==4)&&(((*st)%10)==4))
				    {
				      uRefU=localPos.x(); 
				      localPos.y();
				    }


				}
			      else
				{
				  hitty=false;
				}
			      
			      if (hitty) 
				{
				  //			      trackPoints2Check=true;
				  fracHit++;
				}
			    }
			}
		      cutoff=(float) fracHit/81.;
		      //		      if (cutoff>0.999)
			
			
		    }
		  else
		    {
		      if (whichVerb>5)
			std::cout<<"sigma: "<<sigmas<<" disregarded as too small"<<std::endl;
		      extrapPos.setX(extrapol.X());
		      if (fabs(test.Y0())>offset56)
			{
			  extrapPos.setY(extrapol.Y());
			  localPos=Totem_RP_geometry->GlobalToLocal(near,extrapPos);
			  hitty=foTopo.IsHit(localPos.x(),localPos.y());
			}
		      else
			{
			  hitty=false;
			}
		      //double a=fabs(localPos.x());
		      //double b=fabs(localPos.y());
		      //cutoff= ((a<25)&&(b<25)) ? 1 : 0;
			
		      if (((*st)%10)==5)
			{
			  uRef=localPos.x(); 
			  vRef=localPos.y();
			}
		      if (((*st)%10)==4)
			{
			  uRefU=localPos.x(); 
			  localPos.y();
			}

		      cutoff= (hitty ? 1. : 0.);
		      //		      if (hitty)
			
			

		    }
		
		  if (whichVerb>10)
		    std::cout<<" in det: "<<(*st)<<" u:"<<localPos.x()<<"+- "<<errPX<<" v:"<<localPos.y()<<" +-"<<errPY<<" tx:"<<test.GetTx()<<" ty:"<<test.GetTy();

		 
		  fracEvt+=cutoff/10.;
		  if (whichVerb>10)
		    std::cout<<"track from refPot hits checkpot detector "<<(*st)<<" :"<<((float) 100.*fracHit/81.)<< "% of the time"<<std::endl;
		  if (whichVerb>10)
		    std::cout<<std::endl;
		}

	      if (fracEvt>0.9999)
		{
		  srcUpos->Fill(fixLocalPos.x());
		  srcVpos->Fill(fixLocalPos.y());
		  srcUpos->Fill(fixLocalPos2.y());
		  srcVpos->Fill(fixLocalPos2.x());

		  trkTx->Fill(test.GetTx());
		  trkTy->Fill(test.GetTy());

		  trackPoints2Check=true;
		  if (fracEvt>1.001)
		    std::cout<<"unphysical fracEvt above 100% !!!"<<std::endl;
		}
	      else
		{
		  trackPoints2Check=false;
		}
	      //	      if (trackPoints2Check)
	      //	      	refEvts++;
	      //	      //	      refEvts+=fracEvt;



	    }
	}
    }
  for (int u=0;u<12;u++)
    {
      if (!HzAllPotSTrk[u])
	{
	  outputPerPot[u]->Fill(17.);
	  if ((!HzAllPotMoreMTrk[u])&&(!HzAllPot1MTrk[u]))
	    outputPerPot[u]->Fill(20.);
	    
	}
      
    }
 
  /*
  
  //  if (trackPoints2Check&&checkTrack)
  if (refTrack&&checkTrack&&trackPoints2Check)
    {
      if (whichVerb>7)
	{
	  std::cout<<"PULL VALUES : sigma x: "<<sigmaX<<" , sigma y: "<<sigmaY<<" , xR:"<<xRefEx<<" , xC:"<<xCheck;
	  std::cout<<" , yR:"<<yRefEx<<" , yC:"<<yCheck<<std::endl;
	}
      //hardcoded number of sigmas in X & Y, will soon add param(?)
      double xSig=(xRefEx-xCheck)*2./(fabs(errX)+0.0000000001);
      double ySig=(yRefEx-yCheck)*3./(fabs(errY)+0.0000000001);
      if (errX)
	pullX->Fill(xSig);
	
      //	pullX->Fill((xRefEx-xCheck)/sigmaX);
      if (errY)
	pullY->Fill(ySig);
      if (errY&&errX&&(fabs(xCheck)<0.5))
	{
	  
	  pullXpY->Fill(xSig,ySig);
	  pullXXc->Fill(xSig,xCheck);

	  if (fabs(xCheck)<0.5)
	    pullYYc->Fill(ySig,yCheck);
	  pullYvsXY->Fill(ySig,xCheck,yCheck);
	}
      //	pullY->Fill((yRefEx-yCheck)/sigmaY);

      //      refEvts++;
      //      if ((fabs(xSig)>2)||(fabs(ySig)>3))
      //      if (((fabs(xSig)>2)&&(fabs(xSig)<3)&&(fabs(ySig)<3)) ||  ((fabs(ySig)>3) && (fabs(ySig)<5)&&(fabs(xSig)<2))) //||((fabs(ySig)>3)
      //	checkEvts++;
    }

  */

#ifdef FOdebug    
  std::cout<<"DEBUG: Have read single-track class, next form various boolean combinations"<<std::endl;
#endif




  bool all4=true;
  int potS=0;
  int potM=0;
  int potNoS=0;
  unsigned int vilkenMultiPot=1212;
  unsigned int vilkenNoSPot[4]={1212};
  for (int u=0;u<4;u++)
    {
      all4=all4&&potHasSTrk1[u];
      if (potHasSTrk1[u])
	potS++;
      if (!potHasSTrk1[u])
	{
	  vilkenNoSPot[potNoS]=potsD2[diagonal][u];
	  potNoS++;
	}
      if ((!potHasSTrk1[u])&&potHasMTrk1[u]&&readMultiTrk)
	{
	  potM++;
	  vilkenMultiPot=potsD2[diagonal][u];
	}

    }
  //  if (all4)
  if (potS==4)
    {
      potAll4Trk[falseBunch]++;
      if (!potAll4Trk.count(bunchNumSD))
	{
	  potAll4Trk[bunchNumSD]=1;
	}
      else
	{
	  potAll4Trk[bunchNumSD]++;
	}

    }

  if ((potS==3)&&(potM==1)&&readMultiTrk)
    {
      potTrk3S1M[falseBunch]++;
 
      if (!potTrk3S1M.count(bunchNumSD))
	{
	  potTrk3S1M[bunchNumSD]=1;
	}
      else
	{
	  potTrk3S1M[bunchNumSD]++;
	}

     multiPot->Fill(vilkenMultiPot);
    }

  if ((potS==2)&&(potM==2)&&readMultiTrk)
    {
      potTrk2S2M[falseBunch]++;


      if (!potTrk2S2M.count(bunchNumSD))
	{
	  potTrk2S2M[bunchNumSD]=1;
	}
      else
	{
	  potTrk2S2M[bunchNumSD]++;
	}



      int armA=vilkenNoSPot[0]/100;
      int armB=vilkenNoSPot[1]/100;
      if (armA+armB==1)
	{
	  pot2ATrk2S2M[falseBunch]++;
	  if (!pot2ATrk2S2M.count(bunchNumSD))
	    {
	      pot2ATrk2S2M[bunchNumSD]=1;
	    }
	  else
	    {
	      pot2ATrk2S2M[bunchNumSD]++;
	    }
	  

	}
    }

  bool xCutPass=false;
  if (refTrack&&trackPoints2Check)
    {
      if ((xRefMid>exxLow)&&(xRefMid<exxHi))
	xCutPass=true;
    }

  bool willCrossOneTS=false;

  if (refTrack&&trackPoints2Check)
    {
      refUVsCheckU->Fill(uReallyRef,uRef);
      refVVsCheckV->Fill(vReallyRef,vRef);

      const float ishitBorder=2*16.9; //small overestimate
      const float tmcSectWidth=ishitBorder/16;
      int extrapSector=(-uRef+ishitBorder/2.)/tmcSectWidth;
      int extrapSectorU=(-uRefU+ishitBorder/2.)/tmcSectWidth;
      int refSectorV=(-uReallyRef+ishitBorder/2.)/tmcSectWidth;
      int refSectorU=(-uReallyRefU+ishitBorder/2.)/tmcSectWidth;
      if (checkTrack)
	{
	  int chkSectorV=(-uCheck+ishitBorder/2.)/tmcSectWidth;
	  int chkSectorU=(-uCheckU+ishitBorder/2.)/tmcSectWidth;
	  refTSvsCheckU->Fill(refSectorU,chkSectorU);
	  refTSvsCheckV->Fill(refSectorV,chkSectorV);
	}
      if (extrapSector!=refSectorV)
	{
	  refTSvsExtrap->Fill(refSectorV,extrapSector);
	  willCrossOneTS=true;
	}
 
      if (extrapSectorU!=refSectorU)
	{
	  refTSvsExtrapU->Fill(refSectorU,extrapSectorU);
	  willCrossOneTS=true;
	}
     
      //plane 5 is a V-plane, so the first local coordinate is I think V, so ask for delta"V"=delta"1st"=deltaU=uRef-uCheck
      //abs larger than 32 strips of 66mu each -> 
      //      if (fabs(uRef-uReallyRef)>(32*0.066))
      //	willCrossOneTS=true;

    }

  if (trackPoints2Check&&willCrossOneTS)
    {
      //      refEvts++;
      refEvtXY->Fill(xRefMid,yRefMid);
    }

  bool UTrig=false;
  bool VTrig=false;
  int SrcSectorU=20;
  int SrcSectorV=20;
  int nonPointTSU=20;
  int nonPointTSV=20;
  bool trigUSrc=false;
  bool trigVSrc=false;
  bool trigUtrg=false;
  bool trigVtrg=false;
  bool elastU=false;
  bool elastV=false;
  bool refTrkEvt=false;
  int UPotOn[4]={0};
  int VPotOn[4]={0};
  int UPotOff[4]={0};
  int VPotOff[4]={0};




  int allPotU[8]={44};
  int allPotV[8]={44};
  //MTrCand01AnyPot

  int HzAllPotU[12]={44};
  int HzAllPotV[12]={44};

  //bool checkOn=false;
  std::bitset<16> elBits,trgElBits;
  elBits.reset();
  elBits.set(12);
  if ((refPot==120)||(refPot==124))
    elBits.set(11);
  if (refPot==124)
    elBits.set(13);
  trgElBits.reset();
  for (int i=1;i<15;i++)
    {
      bool left1=elBits.test(i-1);
      bool middle1=elBits.test(i);
      bool right1=elBits.test(i+1);
      if (left1||middle1||right1)
	trgElBits.set(i);
    }
  //  std::cout<<"trigger target matching bits: "<<trgElBits<<std::endl;

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, loop twice over TMC bits"<<std::endl;
#endif




  for (int y=0;y<2;y++)
    {

      if (y&&trackPoints2Check&&elastU&&elastV&&xCutPass)
	{
	  if ((refPot!=124)&&(refPot!=120))
	    {
	      refTrkEvt=true;
	      //	      refEvts++;
	      if (whichVerb>5)
		std::cout<<"elastic reference event in refPot TMC: U="<<SrcSectorU<<" and V="<<SrcSectorV<<std::endl;
	    }
	  else
	    {
	      if (SrcSectorU+1==SrcSectorV)
		{
		  refTrkEvt=true;

		  //		  refEvts++;
		  if (whichVerb>5)
		    std::cout<<"elastic reference event in refPot TMC (U+1==V): U="<<SrcSectorU<<" and V="<<SrcSectorV<<std::endl;

		}
	    }
	}
    
      for ( TMCIter=(*inputTMC).begin(); TMCIter!=(*inputTMC).end();TMCIter++)
	{
	  RPCCId cc((*TMCIter).getId());
	  std::bitset<16> tmc=(*TMCIter).getBS();
	  unsigned int bitsetPot=100*cc.Arm()+10*cc.Station()+cc.RomanPot();
	  //	  if  ((!y)&&(tmc.count())&&(refPot/100==cc.Arm())&& (((refPot%100)/10)==cc.Station())&&(refPot%10==cc.RomanPot()))
	  tmc.count();
	  
	  if (!y)
	    {
	      for (int t=0;t<8;t++)
		{
		  if ((bitsetPot==potsAll[t])&&(!tmc.count()))
		    {
		      ( cc.IsStripsCoordinateUDirection() ? allPotU[t]=0 : allPotV[t]=0);
		    }
		  if ((bitsetPot==potsAll[t])&&(tmc.count()==1))
		    {
		      ( cc.IsStripsCoordinateUDirection() ? allPotU[t]=1 : allPotV[t]=1);
		    }
		}

	      for (int uu=0;uu<12;uu++)
		{
		  if (bitsetPot==potsAllHz[uu])
		    {
		      ( cc.IsStripsCoordinateUDirection() ? HzAllPotTmcU[uu] : HzAllPotTmcV[uu])=tmc.count();
		      if (!tmc.count())
			( cc.IsStripsCoordinateUDirection() ? HzAllPotU[uu]=0 : HzAllPotV[uu]=0);
		      if (tmc.count()==1)
			( cc.IsStripsCoordinateUDirection() ? HzAllPotU[uu]=1 : HzAllPotV[uu]=1);
		      if (tmc.count()>1)
			( cc.IsStripsCoordinateUDirection() ? HzAllPotU[uu]=2 : HzAllPotV[uu]=2);
		      if (tmc.count()>4)
			( cc.IsStripsCoordinateUDirection() ? HzAllPotU[uu]=4 : HzAllPotV[uu]=4);
		    }
		}


	    }

	  if (y&&(tmc.count()==1))//((tmc&trgElBits).count()))
	    {
	      for (int ii=0;ii<4;ii++)
		{
		  if (potsD2[diagonal][ii]==bitsetPot)
		    {
		      ( cc.IsStripsCoordinateUDirection() ? UPotOn[ii]=1 : VPotOn[ii]=1);
		    }
		}
	    }

	  if (y&&(!tmc.count()))//((tmc&trgElBits).count()))
	    {
	      for (int ii=0;ii<4;ii++)
		{
		  if (potsD2[diagonal][ii]==bitsetPot)
		    {
		      ( cc.IsStripsCoordinateUDirection() ? UPotOff[ii]=tmc.count() : VPotOff[ii]=tmc.count());
		    }
		}
	    }

	  if (y&&(tmc.count()>1))//((tmc&trgElBits).count()))
	    {
	      for (int ii=0;ii<4;ii++)
		{
		  if (potsD2[diagonal][ii]==bitsetPot)
		    {
		      int delta;
		      ((tmc.count()>4) ? delta=2 : delta=0);
		      ( cc.IsStripsCoordinateUDirection() ? UPotOn[ii]=2+delta : VPotOn[ii]=2+delta);
		    }
		}
	    }


	  //	  if  ((!y)&&trackPoints2Check&&(tmc.count()==1)&&(refPot/100==cc.Arm())&& (((refPot%100)/10)==cc.Station())&&(refPot%10==cc.RomanPot()))

	  if  ((!y)&&(tmc.count()==1)&&(refPot==bitsetPot)&&(cc.IsStripsCoordinateUDirection()))
	    {
	      for (int u=0;u<16;u++)
		{
		  if (tmc.test(u))
		    {
		      UvsTMCu->Fill(uReallyRef,u);
		      VvsTMCu->Fill(vReallyRef,u);
		      (refTrack ? trksWOneTmcBitOn : NotrksWOneTmcBitOn)->Fill(u);
		      //fill histo only for V-planes
		    }
		}
	    }

	  if  ((!y)&&(tmc.count()==1)&&(refPot==bitsetPot)&&(cc.IsStripsCoordinateVDirection()))
	    {
	      for (int u=0;u<16;u++)
		{
		  if (tmc.test(u))
		    {
		      UvsTMCv->Fill(uReallyRef,u);
		      VvsTMCv->Fill(vReallyRef,u);
		      //		      (refTrack ? trksWOneTmcBitOn : NotrksWOneTmcBitOn)->Fill(u);
		      //fill histo only for V-planes
		    }
		}
	    }




	  if  (y&&(tmc.count()==1)&&(refPot==bitsetPot)&&(tmc&elBits).count())
	    {
	      for (int u=0;u<16;u++)
		{
		  if (tmc.test(u))
		    ( cc.IsStripsCoordinateUDirection() ? nonPointTSU : nonPointTSV) = u;
		}  
		    
	    }


	  if  ((!y)&&trackPoints2Check&&(tmc.count()==1)&&(refPot==bitsetPot))
	    {
	      (cc.IsStripsCoordinateUDirection() ? trigUSrc : trigVSrc) = true;
	      if (whichVerb>5)
		std::cout<<"ref track trigger sector on: ";

	      if (whichVerb>5)
		std::cout<<"ref A:"<<cc.Arm()<<" S:"<<cc.Station()<<" RP:"<<cc.RomanPot()<<" dir:"<<cc.IsStripsCoordinateUDirection()<<std::endl;

	      if ((tmc&elBits).count())
		{
		  if ((refPot!=120)&&(refPot!=124))
		    {
		      ( cc.IsStripsCoordinateUDirection() ? elastU : elastV) = true; 
		    }
		  else
		    {
		      bool temptemp=tmc.test( cc.IsStripsCoordinateUDirection() ? 12 : 13);
		      ( cc.IsStripsCoordinateUDirection() ? elastU : elastV) = (tmc.test( cc.IsStripsCoordinateUDirection() ? 11 : 12)||(refPot==124 ? temptemp : false)); 
		      //		      if (refPot==124)
		      //			( cc.IsStripsCoordinateUDirection() ? elastU : elastV) = tmc.test( cc.IsStripsCoordinateUDirection() ? 12 : 13); 

		 
		    }
		}

	      for (int u=0;u<16;u++)
		{
		  if (tmc.test(u))
		    {
		      if (whichVerb>5) 
			std::cout<<" "<<u;
		      ( cc.IsStripsCoordinateUDirection() ? SrcSectorU : SrcSectorV) = u;
		      ( cc.IsStripsCoordinateUDirection() ? srcU->Fill(u) : srcV->Fill(u));
		      
		    }
		  
		}
	      if (whichVerb>5) 
		std::cout<<std::endl;
	    }
	
	  if (y&&(tmc.count())&&(checkPot==bitsetPot) )
	    {
	      if (whichVerb>5)
		std::cout<<"RPCCBits.first:"<<cc<<" .2nd.count():"<<(*TMCIter).getBS().count()<<std::endl;
	      if (whichVerb>5)
		std::cout<<"A:"<<cc.Arm()<<" S:"<<cc.Station()<<" RP:"<<cc.RomanPot()<<" dir:"<<cc.IsStripsCoordinateUDirection()<<std::endl;
	      (cc.IsStripsCoordinateUDirection() ? UTrig : VTrig)=true;
	      
	      if (trackPoints2Check)
		{
		  if (whichVerb>5)
		    {
		      std::cout<<"ref track points to checkPot";
		      if (cc.IsStripsCoordinateUDirection() ? trigUSrc : trigVSrc)
			std::cout<<" and source pot has exactly one trigSector on";
		      std::cout<<std::endl;
		      std::cout<<"target track trigger sector on: ";
		    }
		  if ((tmc&trgElBits).count())
		    {
		      for (int u=0;u<16;u++)
			{
			  if (tmc.test(u))
			    //		      if (cc.Direction())
			    {
			      if (whichVerb>5)
				std::cout<<" "<<u;
			      int a=(cc.IsStripsCoordinateUDirection() ? SrcSectorU : SrcSectorV);
			      //if ((u==a+1)||(u==a-1))
			      if ((u==a)||(u==a+1)||(u==a-1))
				{
				  (cc.IsStripsCoordinateUDirection() ? trigUtrg : trigVtrg)=true;
				  
				  ( cc.IsStripsCoordinateUDirection() ? trgU->Fill(u) : trgV->Fill(u));
				}
			      
			    }
			}
		    }
		  if (whichVerb>5)
		    std::cout<<std::endl;



		  if (cc.IsStripsCoordinateUDirection())
		    {
		      trgUexPos->Fill(vRef); //was filled for last detector, symbID = ###9, i.e. a V-plane, for which local-coord first 
		      //coordinate is V, second is U
		    }
		  else
		    {
		      trgVexPos->Fill(uRef);
		    }
		}
	      
	      
	    }
	}
    }

  int potsEmpty=0;

  for (int u=0;u<12;u++)
    {
      if ((HzAllPotU[u]==0)&&(HzAllPotV[u]==0))
	potsEmpty++;
      if ((HzAllPotU[u]==1)&&(HzAllPotV[u]==1))
	outputPerPot[u]->Fill(0.);
      if ((!HzAllPotU[u])&&(!HzAllPotV[u]))
	{
	  outputPerPot[u]->Fill(1.);
	  if (HzAllPotSTrk[u])
	    outputPerPot[u]->Fill(8.);
	}
      if (HzAllPotU[u]+HzAllPotV[u]==1)
	{
	  outputPerPot[u]->Fill(2.);
	  if (HzAllPotSTrk[u])
	    outputPerPot[u]->Fill(5.);
	  if (readMultiTrk&&HzAllPot1MTrk[u])
	    outputPerPot[u]->Fill(6.);
	  if (readMultiTrk&&HzPotM0S[u])
	    outputPerPot[u]->Fill(7.);
	}
      if ((HzAllPotU[u]==2)&&(HzAllPotV[u]==2))
	outputPerPot[u]->Fill(3.);
      if ((HzAllPotU[u]==4)&&(HzAllPotV[u]==4))
	outputPerPot[u]->Fill(4.);
    }

  if (potsEmpty==12)
    {
      AllPot00[falseBunch]++;

      if (!AllPot00.count(bunchNumSD))
	{
	  AllPot00[bunchNumSD]=1;
	}
      else
	{
	  AllPot00[bunchNumSD]++;
	}

    }
  else
    {
      AllPotUorV[falseBunch]++;
      if (!AllPotUorV.count(bunchNumSD))
	{
	  AllPotUorV[bunchNumSD]=1;
	}
      else
	{
	  AllPotUorV[bunchNumSD]++;
	}

    }
  //  bool allOn,allOff;

  //  allOn=UPotOn[0]&&VPotOn[0];
  //  for (int u=1;u<4;u++)
  //    allOn=allOn&&UPotOn[u]&&VPotOn[u];

  //  allOff=(!UPotOn[0])&&(!VPotOn[0]);
  //  for (int u=1;u<4;u++)
  //    allOff=allOff&&(!UPotOn[u])&&(!VPotOn[u]);

  //  if (trackPoints2Check&&elastU&&elastV&&xCutPass)
  //    {


#ifdef FOdebug    
  std::cout<<"DEBUG: Next form nondiagonal / all-8 / all-12 RP combinations"<<std::endl;
#endif



  refEvtsVale++;
  int numOriOn=0;
  int potsOn=0;

  for (int t=0;t<8;t++)
    {
      if ((allPotU[t]+allPotV[t])==1)
	{
	  num01->Fill(potsAll[t]);
	  AnyPot01++;
	  if (potHasMTrkCand[t]&&readMultiTrk)
	    {
	      MTrCand01AnyPot++;
	      mul01cand->Fill(potsAll[t]);
	    }
	  if (potHasSTrkCand[t])
	    {
	      TrCand01AnyPot++;
	      sing01cand->Fill(potsAll[t]);
	    }
	  if (AllPotSTrk[t])
	    {
	      sing01->Fill(potsAll[t]);
	      STr01AnyPot++;
	    }
	  if (AllPotMTrk[t])
	    {
	      mul01->Fill(potsAll[t]);
	      MTr01AnyPot++;
	    }
	}
    }
    
  int potsOff=0;
  int potMaxTMC=0;
  int potMediumT=0;
  int potHalf=0;
  unsigned int whichPotHalf[4]={0};
  unsigned int whichPotOff[4]={0};
  for (int u=0;u<4;u++)
    {
      if (UPotOn[u]==1)
	numOriOn++;
      if (VPotOn[u]==1)
	numOriOn++;
      if ((UPotOn[u]==1)&&(VPotOn[u]==1))
	potsOn++;
      if ((!UPotOn[u])&&(!VPotOn[u]))
	{
	  whichPotOff[potsOff]=potsD2[diagonal][u];
	  potsOff++;
	}
      if ((UPotOn[u]==2)&&(VPotOn[u]==2))
	potMediumT++;
      if ((UPotOn[u]==4)&&(VPotOn[u]==4))
	potMaxTMC++;
      //      if ((UPotOn[u]!=VPotOn[u])&&(UPotOn[u]<2)&&(VPotOn[u]<2))
      if (UPotOn[u]+VPotOn[u]==1)
	{
	  whichPotHalf[potHalf]=potsD2[diagonal][u];
	  potHalf++;
	}
    }

  bool armTZ[2]={false};
  int armTS[2]={0};
  if ((potHasSTrk1[0]!=potHasSTrk1[1])&&(potHasSTrk1[2]!=potHasSTrk1[3]))
    {
      for (int u=0;u<2;u++)
	{
	  if (potHasSTrk1[2*u]&&((UPotOff[2*u+1]+VPotOff[2*u+1])<2))
	    {
	      armTZ[u]=true;
	      armTS[u]=UPotOff[2*u+1]+VPotOff[2*u+1];
	    }
	  if (potHasSTrk1[2*u+1]&&((UPotOff[2*u]+VPotOff[2*u])<2))
	    {
	      armTZ[u]=true;
	      armTS[u]=UPotOff[2*u]+VPotOff[2*u];
	    }
	}
    }

  if (armTZ[0]&&armTZ[1]&&((armTS[0]+armTS[1])<2))
    {
      pot2T2Empty[falseBunch]++;
      if (!pot2T2Empty.count(bunchNumSD))
	{
	  pot2T2Empty[bunchNumSD]=1;
	}
      else
	{
	  pot2T2Empty[bunchNumSD]++;
	}
    }


  if (numOriOn==8)
    {
      potAll4On[falseBunch]++;

      if (!potAll4On.count(bunchNumSD))
	{
	  potAll4On[bunchNumSD]=1;
	}
      else
	{
	  potAll4On[bunchNumSD]++;
	}


      //      outputFulfillEvts->Fill(0);
    }
  //  if (numOriOn==7)
  if ((potsOn==3)&&(potHalf==1))
    {


     if (!pots3On1HalfAll.count(bunchNumSD))
	{
	  pots3On1HalfAll[bunchNumSD]=1;
	}
      else
	{
	  pots3On1HalfAll[bunchNumSD]++;
	}


      pots3On1HalfAll[falseBunch]++;
      oneOriOnly3p1->Fill(whichPotHalf[0]);
      for (int y=0;y<4;y++)
	{
	  if ((potsD2[diagonal][y]==whichPotHalf[0])&&lineInUOnly[y])
	    {
	      pot3On1HHu[falseBunch]++;

	      if (!pot3On1HHu.count(bunchNumSD))
		{
		  pot3On1HHu[bunchNumSD]=1;
		}
	      else
		{
		  pot3On1HHu[bunchNumSD]++;
		}


	    }

	  if ((potsD2[diagonal][y]==whichPotHalf[0])&&lineInVOnly[y])
	    {
	      pot3On1HHv[falseBunch]++;
	      if (!pot3On1HHv.count(bunchNumSD))
		{
		  pot3On1HHv[bunchNumSD]=1;
		}
	      else
		{
		  pot3On1HHv[bunchNumSD]++;
		}


	    }
	  
	  if ((potsD2[diagonal][y]==whichPotHalf[0])&&potHasSTrk1[y])
	    oneOriOnly3p011S->Fill(whichPotHalf[0]);

	  if ((potsD2[diagonal][y]==whichPotHalf[0])&&potHasMTrk1[y]&&readMultiTrk&&(!potHasSTrk1[y]))
	    {
	      pots3On1Half[falseBunch]++;


	      if (!pots3On1Half.count(bunchNumSD))
		{
		  pots3On1Half[bunchNumSD]=1;
		}
	      else
		{
		  pots3On1Half[bunchNumSD]++;
		}



	      mul01diag3p1->Fill(whichPotHalf[0]);
	    }
	}
    }
  if ((potsOn==3)&&(potsOff==1))
    {


     if (!pots3On1Off.count(bunchNumSD))
	{
	  pots3On1Off[bunchNumSD]=1;
	}
      else
	{
	  pots3On1Off[bunchNumSD]++;
	}


      pots3On1Off[falseBunch]++;





      for (int y=0;y<4;y++)
	{
	  if ((potsD2[diagonal][y]==whichPotOff[0])&&potHasSTrk1[y])
	    uv11Times3Plus001S->Fill(whichPotOff[0]);
	} 
    }
  if ((potsOn==3)&&(potMaxTMC==1))
    {
      pot3Plus1Max[falseBunch]++;

     if (!pot3Plus1Max.count(bunchNumSD))
	{
	  pot3Plus1Max[bunchNumSD]=1;
	}
      else
	{
	  pot3Plus1Max[bunchNumSD]++;
	}


    }
  if ((potsOn==3)&&(potMediumT==1))
    {
      pot3Plus1Medium[falseBunch]++;


    if (!pot3Plus1Medium.count(bunchNumSD))
	{
	  pot3Plus1Medium[bunchNumSD]=1;
	}
      else
	{
	  pot3Plus1Medium[bunchNumSD]++;
	}




    }

  if ((potsOn==2)&&(potMediumT==2))
    {
      pot2Plus2Medium[falseBunch]++;

    if (!pot2Plus2Medium.count(bunchNumSD))
	{
	  pot2Plus2Medium[bunchNumSD]=1;
	}
      else
	{
	  pot2Plus2Medium[bunchNumSD]++;
	}


    }
  //  if (numOriOn==0)
  if (potsOff==4)
    {


     if (!potAll4Off.count(bunchNumSD))
	{
	  potAll4Off[bunchNumSD]=1;
	}
      else
	{
	  potAll4Off[bunchNumSD]++;
	}


      potAll4Off[falseBunch]++;
      //      outputFulfillEvts->Fill(1);
      if (potsEmpty<12)
	{
	  OneDiag00SomepotOn[falseBunch]++;
	  if (!OneDiag00SomepotOn.count(bunchNumSD))
	    {
	      OneDiag00SomepotOn[bunchNumSD]=1;
	    }
	  else
	    {
	      OneDiag00SomepotOn[bunchNumSD]++;
	    }

	  
	}
    }
  //    }
  if ((potsOn==2)&&(potsOff==2))
    {
      TwoOnTwoOff[falseBunch]++;

     if (!TwoOnTwoOff.count(bunchNumSD))
	{
	  TwoOnTwoOff[bunchNumSD]=1;
	}
      else
	{
	  TwoOnTwoOff[bunchNumSD]++;
	}


 
    }
  if ((potsOn==2)&&(potHalf==2))
    {
      TwoOnTwoHalfAll[falseBunch]++;
      if (!TwoOnTwoHalfAll.count(bunchNumSD))
	{
	  TwoOnTwoHalfAll[bunchNumSD]=1;
	}
      else
	{
	  TwoOnTwoHalfAll[bunchNumSD]++;
	}




      oneOriOnly2p2->Fill(whichPotHalf[0]);
      oneOriOnly2p2->Fill(whichPotHalf[1]);
      for (int y=0;y<4;y++)
	{

	  //	  if ((potsD2[y]==whichPotHalf[0])&&potHasMTrk1[y])
	  //if ((potsD2[y]==whichPotHalf[1])&&potHasMTrk1[y])

	  if ((potsD2[diagonal][y]==whichPotHalf[0])&&potHasMTrk1[y]&&readMultiTrk)
	    {
	      for (int t=0;t<4;t++)
		{
		  if ((potsD2[diagonal][t]==whichPotHalf[1])&&potHasMTrk1[t])
		    {
		      mul01diag2p2->Fill(whichPotHalf[0]);
		      mul01diag2p2->Fill(whichPotHalf[1]);

		      TwoOnTwoHalf[falseBunch]++;

		      if (!TwoOnTwoHalf.count(bunchNumSD))
			{
			  TwoOnTwoHalf[bunchNumSD]=1;
			}
		      else
			{
			  TwoOnTwoHalf[bunchNumSD]++;
			}



		      if (t==y)
			std::cout<<"Found 2*(1,1)+2*(0,1), but the last two pots are the same pot!!"<<std::endl;
		    }
		}
	    }
	}
    }
  if ((potsOn==2)&&(potMaxTMC==2))
    {
      TwoOnTwoVMul[falseBunch]++;
      if (!TwoOnTwoVMul.count(bunchNumSD))
	{
	  TwoOnTwoVMul[bunchNumSD]=1;
	}
      else
	{
	  TwoOnTwoVMul[bunchNumSD]++;
	}

    }

  bool TMCtrk[4]={false};
  for (int u=0;u<4;u++)
    {
      if ((UPotOn[u]==1)&&(VPotOn[u]==1))
	TMCtrk[u]=true;
    }

  if ((TMCtrk[0]!=TMCtrk[1])&& (TMCtrk[2]!=TMCtrk[3]))
    {
      if ((potsOn==2)&&(potMediumT==2))
	{
	  pot2A2Plus2Medium[falseBunch]++;

	  if (!pot2A2Plus2Medium.count(bunchNumSD))
	    {
	      pot2A2Plus2Medium[bunchNumSD]=1;
	    }
	  else
	    {
	      pot2A2Plus2Medium[bunchNumSD]++;
	    }



	}
      if ((potsOn==2)&&(potsOff==2))
	{
	  TwoOn2ATwoOff[falseBunch]++;


	  if (!TwoOn2ATwoOff.count(bunchNumSD))
	    {
	      TwoOn2ATwoOff[bunchNumSD]=1;
	    }
	  else
	    {
	      TwoOn2ATwoOff[bunchNumSD]++;
	    }




	}


      if ((potsOn==2)&&(potHalf==2))
	{
	  TwoOn2ATwoHalfAll[falseBunch]++;


	  if (!TwoOn2ATwoHalfAll.count(bunchNumSD))
	    {
	      TwoOn2ATwoHalfAll[bunchNumSD]=1;
	    }
	  else
	    {
	      TwoOn2ATwoHalfAll[bunchNumSD]++;
	    }



	  for (int y=0;y<4;y++)
	    {
	      if ((potsD2[diagonal][y]==whichPotHalf[0])&&potHasMTrk1[y])
		{
		  for (int t=0;t<4;t++)
		    {
		      if ((potsD2[diagonal][t]==whichPotHalf[1])&&potHasMTrk1[t]&&readMultiTrk)
			{
			  TwoOn2ATwoHalf[falseBunch]++;

			  if (!TwoOn2ATwoHalf.count(bunchNumSD))
			    {
			      TwoOn2ATwoHalf[bunchNumSD]=1;
			    }
			  else
			    {
			      TwoOn2ATwoHalf[bunchNumSD]++;
			    }



			  if (t==y)
			    std::cout<<"Found 1*(1,1)+1*(0,1) in each arm, but the two half-on pots are the same pot!!"<<std::endl;
			}
		    }
		}
	    }
 
	}
      if ((potsOn==2)&&(potMaxTMC==2))
	{
	  TwoOn2ATwoVMul[falseBunch]++;

	  if (!TwoOn2ATwoVMul.count(bunchNumSD))
	    {
	      TwoOn2ATwoVMul[bunchNumSD]=1;
	    }
	  else
	    {
	      TwoOn2ATwoVMul[bunchNumSD]++;
	    }
	  

	}

    }

 
#ifdef FOdebug    
  std::cout<<"DEBUG: all-pot combinations done"<<std::endl;
#endif



  if ((whichVerb>5)&&elastU&&elastV)
    std::cout<<"only ask elastU & V in refPot TMC: U="<<SrcSectorU<<" and V="<<SrcSectorV<<std::endl;


  if (((nonPointTSU!=20)||(nonPointTSV!=20))&& (whichVerb>5) )
    {
      std::cout<<"Refpot has one sector on in U:"<<nonPointTSU<<" and V:"<<nonPointTSV<<" . Track points? tp2chk="<<trackPoints2Check;
      std::cout<<" . pointing sectors are U="<<SrcSectorU<<" and V="<<SrcSectorV<<" and elastic TMC is elU="<<elastU<<" and elV="<<elastV;
      std::cout<<"elBits="<<elBits<<"   and trgElBits="<<trgElBits;
      std::cout<<std::endl;
    }

  double y020,y021,y024,y025,y120,y121,y124,y125,spectParma;
  spectParma=-900;
  bool specElast=false;

 

  if (checkTrack&&trackPoints2Check&&trigUSrc&&trigVSrc)//&&elastU&&elastV)
    {
      if ((refPot!=120)&&(refPot!=124))
	{
	  if (elastU&&elastV)
	    srcX->Fill(xRefMid);
	}
      else
	{
	  if (elastU&&elastV&&(SrcSectorU+1==SrcSectorV))
	    srcX->Fill(xRefMid);
	}
      bool arm56=refPot/100;
      bool bottomPot=refPot%2;
      if (arm56)
	{
	  if (bottomPot)
	    {
	      //refPot=121 or 125
	      y125=((refPot==121) ? yCheck : yRefMid);
	      y121=((refPot==121) ? yRefMid : yCheck);
	      spectParma=cn125*y125+cn121*y121-delta56bt;
	      spectro56bt->Fill(spectParma);
	      if (fabs(spectParma)<lim56bt)
		specElast=true;
	    }
	  else
	    {
	      //add x<0.5mm krav
	      y124=((refPot==120) ? yCheck : yRefMid);
	      y120=((refPot==120) ? yRefMid : yCheck);
	      spectParma=cn124*y124+cn120*y120-delta56tp;
	      if (fabs(xRefMid)<0.5)
		{
		  spectro56tp->Fill(spectParma);
		  if (fabs(spectParma)<lim56tp)
		    specElast=true;
		}
	    }
	}
      else
	{
	  if (bottomPot)
	    {
	      y025=((refPot==21) ? yCheck : yRefMid);
	      y021=((refPot==21) ? yRefMid : yCheck);
	      spectParma=cn025*y025+cn021*y021-delta45bt;
	      spectro45bt->Fill(spectParma);
	      if (fabs(spectParma)<lim45bt)
		specElast=true;
	    }
	  else
	    {
	      y024=((refPot==20) ? yCheck : yRefMid);
	      y020=((refPot==20) ? yRefMid : yCheck);
	      spectParma=cn024*y024+cn020*y020-delta45tp;
	      spectro45tp->Fill(spectParma);
	      if (fabs(spectParma)<lim45tp)
		specElast=true;

	    }
	}
    }


  if (specElast)
    {
      if ((refPot!=120)&&(refPot!=124))
	{
	  if (elastU&&elastV)
	    spectroElastRefX->Fill(xRefMid);
	}
      else
	{
	  if (elastU&&elastV&&(SrcSectorU+1==SrcSectorV))
	    spectroElastRefX->Fill(xRefMid);

	}
      spectroElastUvsV->Fill(SrcSectorU,SrcSectorV);
    }

//  if (modeF==3)
//  if (modeF==4)

  //  if ((UTrig||VTrig)&&trackPoints2Check)
  //  if (UTrig&&VTrig&&trackPoints2Check)
  //    checkEvts++;
 

//  if ((UTrig||VTrig)&&trackPoints2Check)
  //    checkEvts++;

  
  bool elastTest2=(modeF ? (trigUtrg&&trigVtrg) : (trigUtrg || trigVtrg));

  if (modeF==3)
    elastTest2=trigUtrg;
  if (modeF==4)
    elastTest2=trigVtrg;


  //  if (elastTest2&&willCrossOneTS&&trackPoints2Check&&elastU&&elastV)
  //    checkEvts++;


  //  if (refTrack&&checkTrack&&trackPoints2Check&&elastU&&elastV&&elastTest2&&willCrossOneTS)



  //tracking efficiency
  if (refTrkEvt&&elastU&&elastV&&trackPoints2Check&&xCutPass)
    {
      refEvts++;

      //6.6.11
      //      if (checkTrack)
      //	{
      //	  refEvts++;
      if (elastTest2)
	checkEvts++;
      //	}

      (elastTest2 ? refEvtsOn++ : refEvtsOff++);
      if (checkTrack)
	{
	  trkMisPointXY->Fill(xRefEx-xCheck,yRefEx-yCheck);
	  bool ItHitsX=(fabs(xRefEx-xCheck)<fabs(errX));
	  bool ItHitsY=(fabs(yRefEx-yCheck)<fabs(errY));
	  if (ItHitsX&&ItHitsY)
	    {
	      (elastTest2 ? checkEvtsOn++ : checkEvtsOff++);
	      
	    }
	}
      else
	{
	  (elastTest2 ? noChkTrkOn++ : noChkTrkOff++);
	}
	
 
    }



    




  bool offDia=(SrcSectorU+1==SrcSectorV);
  if (refTrack&&checkTrack&&trackPoints2Check&&elastU&&elastV&&(((refPot==120)||(refPot==124)) ? offDia : true)) //&&elastTest2&&willCrossOneTS)
    {
      if (whichVerb>10)
	{
	  std::cout<<"PULL VALUES : sigma x: "<<sigmaX<<" , sigma y: "<<sigmaY<<" , xR:"<<xRefEx<<" , xC:"<<xCheck;
	  std::cout<<" , yR:"<<yRefEx<<" , yC:"<<yCheck<<std::endl;
	}
      //hardcoded number of sigmas in X & Y, will soon add param(?)
      double xSig=(xRefEx-xCheck)*2./(fabs(errX)+0.0000000001);
      double ySig=(yRefEx-yCheck)*3./(fabs(errY)+0.0000000001);
      if (errX)
	pullX->Fill(xSig);
	
      //	pullX->Fill((xRefEx-xCheck)/sigmaX);
      if (errY&&xCutPass)
	{
	  pullY->Fill(ySig);
	  pullYVsSpec->Fill(ySig,spectParma);
	}
      if (errY&&specElast&&xCutPass)
	{
	  pullYElast->Fill(ySig);
	  elasticRefVsCheckY->Fill(yRefMid,yCheck);
	}
      if (errY&&errX)
	{
	  pullXpY->Fill(xSig,ySig);
	  pullXXc->Fill(xSig,xCheck);
	  //	  pullYYc->Fill(ySig,yCheck);

	  //	  if (fabs(xCheck)<0.5)
	  pullYYc->Fill(ySig,yCheck);
	  //	  pullYvsXY->Fill(ySig,xCheck,yCheck);


	}
      //	pullY->Fill((yRefEx-yCheck)/sigmaY);
      //      if (xCutPass)
	//	{
	  //	  refEvts++;
	  //      //      if ((fabs(xSig)>2)||(fabs(ySig)>3))
	  //	  if (((fabs(xSig)>2)&&(fabs(xSig)<3)&&(fabs(ySig)<3)) ||  ((fabs(ySig)>3) && (fabs(ySig)<5)&&(fabs(xSig)<2))) //||((fabs(ySig)>3)
	  //	    checkEvts++;
	  //	}
    }


  //  if (UTrig&&VTrig)
  //  if (elastU&&elastV)
  //  if (elastTest2)
  if (elastU&&elastV&&(((refPot==120)||(refPot==124)) ? offDia : true))
    {
      if (whichVerb>5)
	std::cout<<"refPot has two TMC bits on in U=V=12"<<std::endl;
      //	std::cout<<"checkPot has at least one TMC bit on in U or V"<<std::endl;
//      if (trackPoints2Check&&(trigUtrg&&trigVtrg))
      if (trackPoints2Check&&elastTest2&&xCutPass) //&&willCrossOneTS)
	{
	  if (whichVerb>5)
	    std::cout<<", checkPot has TS +-1 on in both/either orientation(s), and refPot has a track pointing to checkPot too"<<std::endl;
	  //	  checkEvts++;
	}
    }

  

  // if (UTrig&&VTrig)
  //    {

  //      std::cout<<"checkPot has at least one TMC bit on in both U & V, i.e. should have triggered"<<std::endl;
  //      if (trackPoints2Check)
  //	{
  //	  std::cout<<", and reference pot has a track pointing to checkPot too"<<std::endl;
  //	  checkEvts++;
  //	}
  //    }

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, single-track class is processed (again)"<<std::endl;
#endif



  for (OneTrkIter=(*input1Trk).begin(); OneTrkIter!=(*input1Trk).end();OneTrkIter++)
    {

      if ((*OneTrkIter).first==checkPot)
	{
	  if (whichVerb>10)
	    std::cout<<"Single-trk found in checked pot: "<<(*OneTrkIter).first; //<<", mulReco found "<<OneTrkIter->second.size()<<" tracks";
	  if (trackPoints2Check)
	    {
	      //	      checkEvts++;
	      //	      checkEvts+=fracEvt; //checkPot AND refPot have tracks
	      if (whichVerb>10)
		std::cout<<", and reference pot has a track too"<<std::endl;
	
	      RPFittedTrack test=(*OneTrkIter).second;
	      if (test.IsValid())
		{


		  double detPos; //,errPX,errPY;
		  CLHEP::Hep3Vector extrapPos,localPos;
		  //	      std::vector<CLHEP::Hep3Vector> localPos;
		  //	      localPos.clear();

		  //		  fracEvt=0;
		  std::set<unsigned int> detects = Totem_RP_geometry->DetsInRP(checkPot);
		  for (std::set<unsigned int>::const_iterator st = detects.begin(); st != detects.end(); ++st)
		    {
		      TotRPDetId near = TotRPDetId(checkPot/100, (checkPot%100)/10, checkPot%10, ((*st)%10));
		      detPos=Totem_RP_geometry->GetDetTranslation(near).z();
		      test.Z0();
		      //		  std::cout<<" tot rp id is="<<near.DetectorDecId();
		      test.GetTxSigma();
		      test.GetTySigma();
		      TVector2 extrapol=test.GetTrackPoint(detPos);
		      extrapPos.setZ(detPos);
		      
		      extrapPos.setX(extrapol.X());
		      extrapPos.setY(extrapol.Y());
		      localPos=Totem_RP_geometry->GlobalToLocal(near,extrapPos);
		      if (whichVerb>5)
			{
			  std::cout<<" in checkPot: valid track with local coords";
			  if ((*st)%2)
			    {
			      std::cout<<" in det: "<<(*st)<<" diffU (last det):"<<localPos.y()-vRef<<" diffV:"<<localPos.x()-uRef<<" z:"<<localPos.z()<<" tx:"<<test.GetTx()<<" ty:"<<test.GetTy();
			    }
			  else
			    {
			      std::cout<<" in det: "<<(*st)<<" diffU (last det):"<<localPos.x()-vRef<<" diffV:"<<localPos.y()-uRef<<" z:"<<localPos.z()<<" tx:"<<test.GetTx()<<" ty:"<<test.GetTy();
			    }
			}
		    }

		  if (whichVerb>10)
		    std::cout<<std::endl;

		}

	    }
	  else
	    {
	      if (whichVerb>10)
		std::cout<<", but reference pot has no tracks pointing to checked pot"<<std::endl;
	    }

	}


      bool canHazTrks=false;
      if (readMultiTrk)
	{
	  for (TrkIter=(*inputTrk).begin(); TrkIter!=(*inputTrk).end();TrkIter++)
	    {
	      if ((*OneTrkIter).first==(*TrkIter).first)
		{
		  canHazTrks=true;
		  //single track was also found by multitrk
		  //	      std::cout<<"Mul-trk & single trk both found in pot: "<<(*TrkIter).first<<", mulReco found "<<TrkIter->second.size();
		  //	      std::cout<<" tracks"<<std::endl;
		  
		}
	    }
	  if ((!canHazTrks)&&(whichVerb>5))
	    std::cout<<"Only single trk found in pot: "<<(*OneTrkIter).first<<std::endl;
	}
    }

  int occupF[1280]={0};
  int totRPHits[2][4][8];
  for (int u=0;u<64;u++) 
    totRPHits[(u/32)%2][(u/8)%4][u%8]=0;

  
  unsigned int iArm=0;unsigned int iStation=0;unsigned int iRP=0;unsigned int iPlane=0;

  unsigned int planeEffHit[1280]={0};


  int partialHits[2][4][8];
  int planeHit[1280]={0};
  for (int u=0;u<64;u++) 
    partialHits[(u/32)%2][(u/8)%4][u%8]=0;

  int checkPotHits=0;
  int diagPotHits[4]={0};

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, process clusters/hits class"<<std::endl;
#endif




  if (readClusters)
    {
      for (; DSViter!=(*input).end();DSViter++)
	{
	  //      totalClu++;
	  int detidHits=DSViter->size();
	  //ifdef FOdebug      
	  if ((!(totalEvts%300000)))
	    std::cout<<"num hits within DetSet-RpDigClu: "<<detidHits<<std::endl;
	  //endif
	  DetId dP=DSViter->id;
	  TotRPDetId dW= (TotRPDetId) dP;
	  iArm=dW.Arm();
	  iStation=dW.Station();
	  iRP=dW.RomanPot();
	  iPlane=dW.Detector();
	  
	  unsigned int dsvPot=100*iArm+10*iStation+iRP;


	  for (int uu=0;uu<12;uu++)
	    {
	      if (potsAllHz[uu]==dsvPot)
		HzAllPotClu[uu]+=detidHits;
	    }

	  for (int uu=0;uu<4;uu++)
	    {
	      if (potsD2[diagonal][uu]==dsvPot)
		diagPotHits[uu]+=detidHits;
	    }

	  //      if (refTrack&&detidHits&&(checkPot==100*iArm+10*iStation+iRP))
	  //         checkPotHits+=detidHits;
	  
	  //      if (trackPoints2Check&&detidHits&&(checkPot==100*iArm+10*iStation+iRP))
	  //         checkPotHits+=detidHits;
	  
	  if (trackPoints2Check&&detidHits&&elastU&&elastV&&(checkPot==100*iArm+10*iStation+iRP))
	    checkPotHits+=detidHits;

	  partialHits[iArm][iStation][iRP]=1;
	  //ifdef FOdebug
	  if (!(totalEvts%90000)) 
	    {
	      std::cout<<"RP det id Arm:"<<iArm;
	      std::cout<<" Station:"<<iStation;
	      std::cout<<" RP:"<<iRP;
	      std::cout<<" Detector:"<<dW.Detector()<<std::endl;
	      //" station:"<<dW.Station()<<" pot:"dW.RomanPot()<<" det:"<<dW.Detector()<<std::endl;
	    }
	  //endif

	  planeHit[1000*iArm+100*iStation+10*iRP+iPlane]=1;

	  std::vector<RPDigCluster>::const_iterator cluN=DSViter->begin();
	  for (; cluN!=DSViter->end();cluN++)
	    {

	      //ifdef FOdebug      
	      if (!(totalEvts%30000)) 
		std::cout<<"cluSize: "<<(*cluN).GetNumberOfStrips()<<std::endl;
	      //endif
	      int tmp1=(*cluN).GetNumberOfStrips();
	      OneEvtCluSum++;
	      //	  std::cout<<"cluSize: "<<(*cluN).GetNumberOfStrips()<<" clu1stStrip:"<<(*cluN).StrBeg()<<" pot:"<<(100*iArm+10*iStation+iRP);
	      //	  std::cout<<" plane:"<<iPlane<<std::endl;


	      /*

//9.11.2011: commented out  this snippet, excluded an always-on strip found in last year's data
	      if ((tmp1==1)&&(100*iArm+10*iStation+iRP==25)&&(iPlane==6)&&((*cluN).StrBeg()==198))
		{
		  for (int uu=0;uu<12;uu++)
		    {
		      if (potsAllHz[uu]==25)
			HzAllPotClu[uu]--;
		    }
		  
		}


	      

	      //  if ((tmp1==1)&&(100*iArm+10*iStation+iRP==checkPot)&&(checkPot==25)&&(iPlane==6)&&((*cluN).StrBeg()==198)&&refTrack)
	      //	  if ((tmp1==1)&&(100*iArm+10*iStation+iRP==checkPot)&&(checkPot==25)&&(iPlane==6)&&((*cluN).StrBeg()==198)&&trackPoints2Check)
	      if ((tmp1==1)&&(100*iArm+10*iStation+iRP==checkPot)&&(checkPot==25)&&(iPlane==6)&&((*cluN).StrBeg()==198)&&trackPoints2Check&&elastU&&elastV)
		{
		  checkPotHits--;
		  if (whichVerb>10)
		    {
		      std::cout<<std::endl;
		      std::cout<<"Cluster in pot 025, plane 6, channel 198 disregarded"<<std::endl;
		    }
		}

*/

	      OnePotCluSum[100*iArm+10*iStation+iRP]++;
	      occupF[1000*iArm+100*iStation+10*iRP+iPlane]+=tmp1;
	      cluSize[iArm][iStation][iRP]+=tmp1;
	      planeEffHit[1000*iArm+100*iStation+10*iRP+iPlane]++;
	      
	      prfSClu->Fill(iArm*100+iStation*10+iRP,tmp1);
	      SOccup[iArm][iStation][iRP]+=( (float) (*cluN).GetNumberOfStrips())/5120.;
	    }
      
    
	  totalClu[iArm][iStation][iRP]+=detidHits;
	  totRPHits[iArm][iStation][iRP]+=detidHits;
	}
    }
//    if (refTrack)
//       numHitsInCheckPointing->Fill(checkPotHits);

//  if (trackPoints2Check&&elastU&&elastV)
    //    numHitsInCheckPointing->Fill(checkPotHits);
  //  if (trackPoints2Check&&(!(UTrig||VTrig)))
  //    numHitsInCheckPointing->Fill(checkPotHits);
  if (readClusters)
    {
      if (trackPoints2Check&&elastU&&elastV&&(!elastTest2))
	numHitsInCheckPointing->Fill(checkPotHits);

//      if (trackPoints2Check&&elastU&&elastV&&(checkPotHits==0))
    }
  if ((potS==2)&&(potNoS==2))
    {


     if (!mp2Trk.count(bunchNumSD))
	{
	  mp2Trk[bunchNumSD]=1;
	}
      else
	{
	  mp2Trk[bunchNumSD]++;
	}


      potDiag2STrk++;

      int arm1=vilkenNoSPot[0]/100;
      int arm2=vilkenNoSPot[1]/100;

       for (int u=0;u<4;u++)
	{
	  if ((vilkenNoSPot[0]==potsD2[diagonal][u])&&(!diagPotHits[u])&&readClusters)
	    for (int r=0;r<4;r++)
	      {
		if ((vilkenNoSPot[1]==potsD2[diagonal][r])&&(!diagPotHits[r]))
		  {

		    if (!potTrk2S2NoClu.count(bunchNumSD))
		      {
			potTrk2S2NoClu[bunchNumSD]=1;
		      }
		    else
		      {
			potTrk2S2NoClu[bunchNumSD]++;
		      }

		    potTrk2S2NoClu[falseBunch]++;
		    if (arm1!=arm2)
		      {
			pot2ATrk2S2NoClu[falseBunch]++;
			
			if (!pot2ATrk2S2NoClu.count(bunchNumSD))
			  {
			    pot2ATrk2S2NoClu[bunchNumSD]=1;
			  }
			else
			  {
			    pot2ATrk2S2NoClu[bunchNumSD]++;
			  }


		      }
		  }
	      }
	}
    }
      //      vilkenNoSPot,diagPotHits;
  if ((potS==1)&&(potNoS==3))
    {
      potDiag1STrk++;

     if (!mp1Trk.count(bunchNumSD))
	{
	  mp1Trk[bunchNumSD]=1;
	}
      else
	{
	  mp1Trk[bunchNumSD]++;
	}


    }
 
  if ((potS==0)&&(potNoS==4))
    {
      potDiag0STrk++;
      if (!mp0Trk.count(bunchNumSD))
	{
	  mp0Trk[bunchNumSD]=1;
	}
      else
	{
	  mp0Trk[bunchNumSD]++;
	}

    }
 

  if ((potS==3)&&(potNoS==1))
    {
     if (!mp3Trk.count(bunchNumSD))
	{
	  mp3Trk[bunchNumSD]=1;
	}
      else
	{
	  mp3Trk[bunchNumSD]++;
	}


      potDiag3STrk++;
      for (int u=0;u<4;u++)
	{
	  if ((vilkenNoSPot[0]==potsD2[diagonal][u])&&(!diagPotHits[u])&&readClusters)
	    {
	      potTrk3S1NoClu[falseBunch]++;

	      if (!potTrk3S1NoClu.count(bunchNumSD))
		{
		  potTrk3S1NoClu[bunchNumSD]=1;
		}
	      else
		{
		  potTrk3S1NoClu[bunchNumSD]++;
		}

	    }


	}
      //      vilkenNoSPot,diagPotHits;
    }

  if (refTrkEvt&&elastU&&elastV&&trackPoints2Check&&xCutPass&&(checkPotHits==0)&&readClusters)
    mptyEvts++;


#ifdef FOdebug    
  std::cout<<"DEBUG: Have read clusters class, next code pairing up near+far pot"<<std::endl;
#endif



  double TwoRpX0[2]={0};
  double TwoRpY0[2]={0};
  double TwoRpDZ=0;
  int pairPotsClu[6]={0};
  int pairPotsTmc[6]={0};
  int pairPotsTmcHalf[6]={0};
  bool pairPotsSTrk[6]={false};
  bool pairPotsBadTrk[6]={false};
  bool tempNearO,tempFarO,tempNearA,tempFarA;
  for (int r=0;r<2;r++)
    {
      tempNearO=(HzAllPotTmcU[r]||HzAllPotTmcV[r]);
      tempFarO=(HzAllPotTmcU[r+4]||HzAllPotTmcV[r+4]);
      tempNearA=(HzAllPotTmcU[r]&&HzAllPotTmcV[r]);
      tempFarA=(HzAllPotTmcU[r+4]&&HzAllPotTmcV[r+4]);
      if (tempFarO&&tempNearO)
	pairPotsTmc[r]=10;
      if (tempFarA&&tempNearA)
	pairPotsTmc[r]=11;
      
      if (tempFarO||tempNearO)
	pairPotsTmcHalf[r]=10;
      if (tempFarA||tempNearA)
	pairPotsTmcHalf[r]=11;

     tempNearA=(HzAllPotTmcU[r]>1)&&(HzAllPotTmcV[r]>1);
      tempFarA=(HzAllPotTmcU[r+4]>1)&&(HzAllPotTmcV[r+4]>1);

      if (tempFarA&&tempNearA)
	pairPotsTmc[r]=22;

      if (tempFarA||tempNearA)
	pairPotsTmcHalf[r]=22;

      pairPotsClu[r]=HzAllPotClu[r]+HzAllPotClu[r+4];
      pairPotsSTrk[r]=HzAllPotSTrk[r]&&HzAllPotSTrk[r+4];
      if (pairPotsSTrk[r])
	{
	  TwoRpX0[0]=x0rp[r];
	  TwoRpX0[1]=x0rp[r+4];
	  TwoRpY0[0]=y0rp[r];
	  TwoRpY0[1]=y0rp[r+4];
	  TwoRpDZ=fabs(z0rp[r]-z0rp[r+4]);

	  if (TwoRpDZ>10000)
	    throw cms::Exception("RPDataReduction") << "Distance between near and far pot ("<<TwoRpDZ<<" is too large, should be ~5m ~5000mm!";

	  if (TwoRpDZ<2500)
	    throw cms::Exception("RPDataReduction") << "Distance between near and far pot ("<<TwoRpDZ<<" is too small, should be ~5m ~5000mm!";

	}


      tempNearO=(HzAllPotTmcU[r+6]||HzAllPotTmcV[r+6]);
      tempFarO=(HzAllPotTmcU[r+10]||HzAllPotTmcV[r+10]);
      tempNearA=(HzAllPotTmcU[r+6]&&HzAllPotTmcV[r+6]);
      tempFarA=(HzAllPotTmcU[r+10]&&HzAllPotTmcV[r+10]);
      if (tempFarO&&tempNearO)
	pairPotsTmc[r+3]=10;
      if (tempFarA&&tempNearA)
	pairPotsTmc[r+3]=11;

      if (tempFarO||tempNearO)
	pairPotsTmcHalf[r+3]=10;
      if (tempFarA||tempNearA)
	pairPotsTmcHalf[r+3]=11;

      tempNearA=(HzAllPotTmcU[r+6]>1)&&(HzAllPotTmcV[r+6]>1);
      tempFarA=(HzAllPotTmcU[r+10]>1)&&(HzAllPotTmcV[r+10]>1);

      if (tempFarA&&tempNearA)
	pairPotsTmc[r+3]=22;
      if (tempFarA||tempNearA)
	pairPotsTmcHalf[r+3]=22;



      pairPotsClu[r+3]=HzAllPotClu[6+r]+HzAllPotClu[r+10];
      pairPotsSTrk[r+3]=HzAllPotSTrk[6+r]&&HzAllPotSTrk[r+10];
      if (pairPotsSTrk[r+3])
	{
	  TwoRpX0[0]=x0rp[r+6];
	  TwoRpX0[1]=x0rp[r+10];
	  TwoRpY0[0]=y0rp[r+6];
	  TwoRpY0[1]=y0rp[r+10];
	  TwoRpDZ=fabs(z0rp[r+6]-z0rp[r+10]);
	  if (TwoRpDZ>10000)
	    throw cms::Exception("RPDataReduction") << "Distance between near and far pot ("<<TwoRpDZ<<"(mm) is too large, should be ~5m ~5000mm!";

	  if (TwoRpDZ<2500)
	    throw cms::Exception("RPDataReduction") << "Distance between near and far pot ("<<TwoRpDZ<<"(mm) is too small, should be ~5m ~5000mm!";

	}


      tempNearO=(HzAllPotTmcU[6*r+2]||HzAllPotTmcV[6*r+2]);
      tempFarO=(HzAllPotTmcU[6*r+3]||HzAllPotTmcV[6*r+3]);
      tempNearA=(HzAllPotTmcU[6*r+2]&&HzAllPotTmcV[6*r+2]);
      tempFarA=(HzAllPotTmcU[6*r+3]&&HzAllPotTmcV[6*r+3]);
      if (tempFarO&&tempNearO)
	pairPotsTmc[3*r+2]=10;
      if (tempFarA&&tempNearA)
	pairPotsTmc[3*r+2]=11;

      if (tempFarO||tempNearO)
	pairPotsTmcHalf[3*r+2]=10;
      if (tempFarA||tempNearA)
	pairPotsTmcHalf[3*r+2]=11;

      tempNearA=(HzAllPotTmcU[6*r+2]>1)&&(HzAllPotTmcV[6*r+2]>1);
      tempFarA=(HzAllPotTmcU[6*r+3]>1)&&(HzAllPotTmcV[6*r+3]>1);

      if (tempFarA&&tempNearA)
	pairPotsTmc[3*r+2]=22;

      if (tempFarA||tempNearA)
	pairPotsTmcHalf[3*r+2]=22;



      pairPotsClu[3*r+2]=HzAllPotClu[6*r+2]+HzAllPotClu[6*r+3];
      pairPotsSTrk[3*r+2]=HzAllPotSTrk[6*r+2]&&HzAllPotSTrk[6*r+3];
      if (pairPotsSTrk[r*3+2])
	{
	  TwoRpX0[0]=x0rp[r*6+2];
	  TwoRpX0[1]=x0rp[r*6+3];
	  TwoRpY0[0]=y0rp[r*6+2];
	  TwoRpY0[1]=y0rp[r*6+3];
	  TwoRpDZ=fabs(z0rp[r*6+2]-z0rp[r*6+3]);
	  if (TwoRpDZ>10000)
	    throw cms::Exception("RPDataReduction") << "Distance between near and far pot ("<<TwoRpDZ<<" is too large, should be ~5m ~5000mm!";

	  if (TwoRpDZ<2500)
	    throw cms::Exception("RPDataReduction") << "Distance between near and far pot ("<<TwoRpDZ<<" is too small, should be ~5m ~5000mm!";

	}

      pairPotsBadTrk[r]=(HzAllPotSTrk[r]!=HzAllPotSTrk[r+4]);
      pairPotsBadTrk[r+3]=(HzAllPotSTrk[6+r]!=HzAllPotSTrk[r+10]);
      pairPotsBadTrk[3*r+2]=(HzAllPotSTrk[6*r+2]!=HzAllPotSTrk[6*r+3]);
     }


 
  bool badRp45Trk=false;
  bool badRp56Trk=false;
  
  bool all45LT10Clu=true;
  bool all56LT10Clu=true;

  int numTrk=0;
  int numBadTrk=0;
  int sumClusters=0;
  int whichPair=-1;
  int sumClu45=0;
  int sumClu56=0;
  int numTmc=0;
  int numTmcHalf=0;

  int maxU45=0;
  int maxV45=0;
  int maxU56=0;
  int maxV56=0;

  //SDCluNoTrk, rpTrgSDCluNoTrk

  for (int y=0;y<6;y++)
    {
      if (HzAllPotTmcU[y]>maxU45)
	maxU45=HzAllPotTmcU[y];
      if (HzAllPotTmcV[y]>maxV45)
	maxV45=HzAllPotTmcV[y];
      if (HzAllPotClu[y]>9)
	all45LT10Clu=false;
    }
  for (int y=6;y<12;y++)
    {
      if (HzAllPotTmcU[y]>maxU56)
	maxU56=HzAllPotTmcU[y];
      if (HzAllPotTmcV[y]>maxV56)
	maxV56=HzAllPotTmcV[y];
      if (HzAllPotClu[y]>9)
	all56LT10Clu=false;
 
    }


  for (int r=0;r<6;r++)
    {
      if (r<3)
	sumClu45+=pairPotsClu[r];
      if (r>2)
	sumClu56+=pairPotsClu[r];

      sumClusters+=pairPotsClu[r];

      if (pairPotsTmc[r]>0)
	{
	  numTmc++;
	  if (r<3)
	    rp45Tmc=true;
	  if (r>2)
	    rp56Tmc=true;
	}

      if (pairPotsTmcHalf[r]>0)
	{
	  numTmcHalf++;
	  //	  if (r<3)
	  //	    rp45Tmc=true;
	  //	  if (r>2)
	  //	    rp56Tmc=true;
	}


      if (pairPotsSTrk[r])
	{
	  if (r<3)
	    sumClu45-=pairPotsClu[r];
	  if (r>2)
	    sumClu56-=pairPotsClu[r];

	  whichPair=r;
	  sumClusters-=pairPotsClu[r];
	  numTrk++;
	  if (r<3)
	    rp45Trk=true;
	  else
	    rp56Trk=true;
	}
     if (pairPotsBadTrk[r])
	{
	  if (r<3)
	    sumClu45-=pairPotsClu[r];
	  if (r>2)
	    sumClu56-=pairPotsClu[r];

	  sumClusters-=pairPotsClu[r];
	  numBadTrk++;
	  if (r<3)
	    badRp45Trk=true;
	  else
	    badRp56Trk=true;
	}
    }


  if ((numTrk+numBadTrk)>1)
    {
      rp45Trk=false;
      rp56Trk=false;
      badRp45Trk=false;
      badRp56Trk=false;
    }

  if (rp45Trk||rp56Trk)
    TwoRPTrkNumProt->Fill(numProt);
	
  if ((numTmc>1)||(numTmcHalf>1))
    {
      rp45Tmc=false;
      rp56Tmc=false;
    }
	
  double rpEta=0;
  const double mProt=0.938272;
  if (rp45Trk)
    rpEta=log(2*optEb/mProt);
    
  if (rp56Trk)
    rpEta=-log(2*optEb/mProt);
    

  //T2
  bool bothArms=T2minusTrk&&T2plusTrk&&(rp45Trk||rp56Trk);
  bool sameArm=((T2minusTrk&&rp56Trk)||(T2plusTrk&&rp45Trk))&&(!bothArms);
  bool opposArm=((T2minusTrk&&rp45Trk)||(T2plusTrk&&rp56Trk))&&(!bothArms);
  bool rpOnly=(rp45Trk||rp56Trk)&&(!(T2minusTrk||T2plusTrk));
  bool bothArmsBad=T2minusTrk&&T2plusTrk&&(badRp45Trk||badRp56Trk);
  bool sameArmBad=((T2minusTrk&&badRp56Trk)||(T2plusTrk&&badRp45Trk))&&(!bothArmsBad);
  bool opposArmBad=((T2minusTrk&&badRp45Trk)||(T2plusTrk&&badRp56Trk))&&(!bothArmsBad);
  bool rpOnlyBad=(badRp45Trk||badRp56Trk)&&(!(T2minusTrk||T2plusTrk));

  bool bothArmsTmc=T2minusTrk&&T2plusTrk&&(rp45Tmc||rp56Tmc);
  bool sameArmTmc=((T2minusTrk&&rp56Tmc)||(T2plusTrk&&rp45Tmc))&&(!bothArmsTmc);
  bool opposArmTmc=((T2minusTrk&&rp45Tmc)||(T2plusTrk&&rp56Tmc))&&(!bothArmsTmc);
  bool rpOnlyTmc=(rp45Tmc||rp56Tmc)&&(!(T2minusTrk||T2plusTrk));

  bool elasticBox=false;
  bool sdBox=false;

  bool noT1Trk=(!T1plusTrk)&&(!T1minusTrk);
  //  bool diffTeleArm=((T1plusTrk&&T2minusTrk)||(T1minusTrk&&T2plusTrk));
  bool sameTeleArm=(((!T1plusTrk)&&T2minusTrk)||((!T1minusTrk)&&T2plusTrk))&&(!noT1Trk);
  bool bothT1Trk=T1plusTrk&&T1minusTrk;
  bool eitherT1Trk=((!T1plusTrk)&&T1minusTrk)||(T1plusTrk&&(!T1minusTrk));

if (numTrk+numBadTrk)
  {
    bothArmsTmc=false;
    sameArmTmc=false;
    opposArmTmc=false;
    rpOnlyTmc=false;
  }

  bool oppoSideFC=((all56LT10Clu&&rp45Trk)||(all45LT10Clu&&rp56Trk));

  //  bool noTrkClusters=(sumClusters==0);

  double angleX=(TwoRpX0[1]-TwoRpX0[0])/TwoRpDZ;
  double angleY=(TwoRpY0[1]-TwoRpY0[0])/TwoRpDZ;

  double ipAngleX,ipAngleY,calcTeeLoXi;
  ipAngleX=(optDVx*TwoRpX0[1]-optVx*angleX)/(optLx*optDVx-optDLx*optVx);
  ipAngleY=TwoRpY0[1]/optLy; //changed from near to far pot 5.12.11, because am using nominal optics at 220m = nearer to the far pot
  calcTeeLoXi=-optEb*optEb*(ipAngleX*ipAngleX+ipAngleY*ipAngleY);

  double crudeCollimX=TwoRpX0[0]-20000.*angleX;
  double crudeCollimY=TwoRpY0[0]-20000.*angleY;

#ifdef FOdebug    
  std::cout<<"DEBUG: Have formed SD p+T2 combinations, next code classifying event based on them"<<std::endl;
#endif


  if (readT2)
    { 

      if ((bothArms||rpOnly||sameArm||opposArm)&&(whichPair==-1))
        throw cms::Exception("RPDataReduction") << "2 RP tracks previously found, now missing, according to my IF-statements";

      int index=kT2NA;

      if (bothArmsTmc&&readLoNeg)
	{
	  rpTmcT2TrBoth[falseBunch]++;
	  if (!rpTmcT2TrBoth.count(bunchNumSD))
	    rpTmcT2TrBoth[bunchNumSD]=1;
	  else
	    rpTmcT2TrBoth[bunchNumSD]++;


	  if ((trigStatus/64)%2)
	    {
	      rpTmcT2TrBothT[falseBunch]++;
	      if (!rpTmcT2TrBothT.count(bunchNumSD))
		rpTmcT2TrBothT[bunchNumSD]=1;
	      else
		rpTmcT2TrBothT[bunchNumSD]++;

	      if (trigStatus%4)
		{
		  rpTmcT2TrBothRT[falseBunch]++;
		  if (!rpTmcT2TrBothRT.count(bunchNumSD))
		    rpTmcT2TrBothRT[bunchNumSD]=1;
		  else
		    rpTmcT2TrBothRT[bunchNumSD]++;
		  

		}

	    }
	}



     if (sameArmTmc&&readLoNeg)
	{
	  rpTmcT2TrSame[falseBunch]++;
	  if (!rpTmcT2TrSame.count(bunchNumSD))
	    rpTmcT2TrSame[bunchNumSD]=1;
	  else
	    rpTmcT2TrSame[bunchNumSD]++;


	  if ((trigStatus/64)%2) //T2
	    {
	      rpTmcT2TrSameT[falseBunch]++;
	      if (!rpTmcT2TrSameT.count(bunchNumSD))
		rpTmcT2TrSameT[bunchNumSD]=1;
	      else
		rpTmcT2TrSameT[bunchNumSD]++;

	      if (trigStatus%4) //RP V or RP H
		{
		  rpTmcT2TrSameRT[falseBunch]++;
		  if (!rpTmcT2TrSameRT.count(bunchNumSD))
		    rpTmcT2TrSameRT[bunchNumSD]=1;
		  else
		    rpTmcT2TrSameRT[bunchNumSD]++;
		  

		}

	    }
	}

    
     if (opposArmTmc&&readLoNeg)
	{
	  rpTmcT2TrOppo[falseBunch]++;
	  if (!rpTmcT2TrOppo.count(bunchNumSD))
	    rpTmcT2TrOppo[bunchNumSD]=1;
	  else
	    rpTmcT2TrOppo[bunchNumSD]++;


	  if ((trigStatus/64)%2)
	    {
	      rpTmcT2TrOppoT[falseBunch]++;
	      if (!rpTmcT2TrOppoT.count(bunchNumSD))
		rpTmcT2TrOppoT[bunchNumSD]=1;
	      else
		rpTmcT2TrOppoT[bunchNumSD]++;

	      if (trigStatus%4)
		{
		  rpTmcT2TrOppoRT[falseBunch]++;
		  if (!rpTmcT2TrOppoRT.count(bunchNumSD))
		    rpTmcT2TrOppoRT[bunchNumSD]=1;
		  else
		    rpTmcT2TrOppoRT[bunchNumSD]++;
		  

		}

	    }
	}





     if (rpOnlyTmc&&readLoNeg)
	{
	  rpTmcT2NoTr[falseBunch]++;
	  if (!rpTmcT2NoTr.count(bunchNumSD))
	    rpTmcT2NoTr[bunchNumSD]=1;
	  else
	    rpTmcT2NoTr[bunchNumSD]++;


	  if ((trigStatus/64)%2)
	    {
	      rpTmcT2NoTrT[falseBunch]++;
	      if (!rpTmcT2NoTrT.count(bunchNumSD))
		rpTmcT2NoTrT[bunchNumSD]=1;
	      else
		rpTmcT2NoTrT[bunchNumSD]++;

	      if (trigStatus%4)
		{
		  rpTmcT2NoTrRT[falseBunch]++;
		  if (!rpTmcT2NoTrRT.count(bunchNumSD))
		    rpTmcT2NoTrRT[bunchNumSD]=1;
		  else
		    rpTmcT2NoTrRT[bunchNumSD]++;
		  

		}

	    }
	}




      if (bothArms)
	{
	  
	  if (oppoSideFC)
	    {
	      rpT2TrBothFC[falseBunch]++;
	      if (!rpT2TrBothFC.count(bunchNumSD))
		rpT2TrBothFC[bunchNumSD]=1;
	      else
		rpT2TrBothFC[bunchNumSD]++;
	    }

	  if (T2TrigLN)
	    {
	      rpT2TrBothT2T[falseBunch]++;
	      if (!rpT2TrBothT2T.count(bunchNumSD))
		rpT2TrBothT2T[bunchNumSD]=1;
	      else
		rpT2TrBothT2T[bunchNumSD]++;
	    }

	  rpT2TrBoth[falseBunch]++;
	  if (!rpT2TrBoth.count(bunchNumSD))
	    rpT2TrBoth[bunchNumSD]=1;
	  else
	    rpT2TrBoth[bunchNumSD]++;
	    

	  if (!rp2T2b[whichPair].count(bunchNumSD))
	    rp2T2b[whichPair][bunchNumSD]=1;
	  else
	    rp2T2b[whichPair][bunchNumSD]++;


	  index=kT2both;

	  
	  if (rp56Trk)
	    {
	      SDT2Multi[2]->Fill(T2PlusMulti);
	      SDT2Multi[3]->Fill(T2MinusMulti);
	    }
	  else
	    {
	      SDT2Multi[2]->Fill(T2MinusMulti);
	      SDT2Multi[3]->Fill(T2PlusMulti);
	    }
	  

	}

      if (sameArm)
	{
	  rpT2TrSame[falseBunch]++;


	  if (!rpT2TrSame.count(bunchNumSD))
	    rpT2TrSame[bunchNumSD]=1;
	  else
	    rpT2TrSame[bunchNumSD]++;
	    
	  if (oppoSideFC)
	    {
	      rpT2TrSameFC[falseBunch]++;
	      if (!rpT2TrSameFC.count(bunchNumSD))
		rpT2TrSameFC[bunchNumSD]=1;
	      else
		rpT2TrSameFC[bunchNumSD]++;
	    }

	  if (T2TrigLN)
	    {
	      rpT2TrSameT2T[falseBunch]++;
	      if (!rpT2TrSameT2T.count(bunchNumSD))
		rpT2TrSameT2T[bunchNumSD]=1;
	      else
		rpT2TrSameT2T[bunchNumSD]++;
	    }


	  if (!rp2T2s[whichPair].count(bunchNumSD))
	    rp2T2s[whichPair][bunchNumSD]=1;
	  else
	    rp2T2s[whichPair][bunchNumSD]++;
	   


	  index=kT2same;
	  
	  if (rp45Trk)
	    {
	      SDT2Multi[0]->Fill(T2PlusMulti);
	    }
	  else
	    {
	      SDT2Multi[0]->Fill(T2MinusMulti);
	    }
	}
      
      if (opposArm)
	{
	  rpT2TrOppos[falseBunch]++;

	  if (!rpT2TrOppos.count(bunchNumSD))
	    rpT2TrOppos[bunchNumSD]=1;
	  else
	    rpT2TrOppos[bunchNumSD]++;
	  

	  if (oppoSideFC)
	    {
	      rpT2TrOpposFC[falseBunch]++;
	      if (!rpT2TrOpposFC.count(bunchNumSD))
		rpT2TrOpposFC[bunchNumSD]=1;
	      else
		rpT2TrOpposFC[bunchNumSD]++;
	    }

	  if (T2TrigLN)
	    {
	      rpT2TrOpposT2T[falseBunch]++;
	      if (!rpT2TrOpposT2T.count(bunchNumSD))
		rpT2TrOpposT2T[bunchNumSD]=1;
	      else
		rpT2TrOpposT2T[bunchNumSD]++;
	    }
  

	  if (!rp2T2o[whichPair].count(bunchNumSD))
	    rp2T2o[whichPair][bunchNumSD]=1;
	  else
	    rp2T2o[whichPair][bunchNumSD]++;


	  index=kT2oppo;

	  
	  if (rp56Trk)
	    {
	      SDT2Multi[1]->Fill(T2PlusMulti);
	    }
	  else
	    {
	      SDT2Multi[1]->Fill(T2MinusMulti);
	    }
	  
	}
      
      if (rpOnly)
	{
	  rpSTrT2NoTr[falseBunch]++;

	  if (!rpSTrT2NoTr.count(bunchNumSD))
	    rpSTrT2NoTr[bunchNumSD]=1;
	  else
	    rpSTrT2NoTr[bunchNumSD]++;
	  

	  if (T2TrigLN>0)
	    {
	      rpSTrT2NoTrT2T[falseBunch]++;
	      if (!rpSTrT2NoTrT2T.count(bunchNumSD))
		rpSTrT2NoTrT2T[bunchNumSD]=1;
	      else
		rpSTrT2NoTrT2T[bunchNumSD]++;

	    }


	  if (oppoSideFC)
	    {
	      rpSTrT2NoTrFC[falseBunch]++;
	      
	      if (!rpSTrT2NoTrFC.count(bunchNumSD))
		rpSTrT2NoTrFC[bunchNumSD]=1;
	      else
		rpSTrT2NoTrFC[bunchNumSD]++;
	    }
  

	  if (!rp2T2n[whichPair].count(bunchNumSD))
	    rp2T2n[whichPair][bunchNumSD]=1;
	  else
	    rp2T2n[whichPair][bunchNumSD]++;
	    

	  index=kT2none;
	}



#ifdef FOdebug    
  std::cout<<"DEBUG: Next, process evts with only 2/12 RP in 45 XOR 56"<<std::endl;
#endif




      if (rp45Trk||rp56Trk)
	{




	  int whichPot=12;
	  for (int y=0;y<6;y++)
	    {
	      if (pairPotsSTrk[y])
		whichPot=y;
	    }

	  if (index==kT2oppo)
	    {
	      if (readLoNeg) //&&noTrkClusters)
		{
		  lonegBunch[6]->Fill(bunchNumSD);

		  //

		  if (!sdBunchEv.count(bunchNumSD))
		    {
		      sdBunchEv[bunchNumSD]=1;
		    }
		  else
		    {
		      sdBunchEv[bunchNumSD]++;
		    }


		  if (trigStatus%2)
		    {
		      sdTrigs->Fill(0); // rp220v
		      if (!rpvTrig.count(bunchNumSD))
			{
			  rpvTrig[bunchNumSD]=1;
			}
		      else
			{
			  rpvTrig[bunchNumSD]++;
			}
		    }
		  if ((trigStatus/2)%2)
		    {
		      sdTrigs->Fill(1); // rp220h
		      if (!rphTrig.count(bunchNumSD))
			{
			  rphTrig[bunchNumSD]=1;
			}
		      else
			{
			  rphTrig[bunchNumSD]++;
			}

		    }
		  if ((trigStatus/4)%2)
		    sdTrigs->Fill(2); // rp220cross
		  if ((trigStatus/64)%2)
		    {
		      sdTrigs->Fill(3); // T2
		      if (!t2Trig.count(bunchNumSD))
			{
			  t2Trig[bunchNumSD]=1;
			}
		      else
			{
			  t2Trig[bunchNumSD]++;
			}


		    }
		  if ((trigStatus/256)%2)
		    sdTrigs->Fill(4); // T1
		  if ((trigStatus/512)%2)
		    {
		      sdTrigs->Fill(5); // Bx
		   
		      if (!bxTrig.count(bunchNumSD))
			{
			  bxTrig[bunchNumSD]=1;
			}
		      else
			{
			  bxTrig[bunchNumSD]++;
			}

		    }
		}
	    }

	  if (index==kT2NA)
	    throw cms::Exception("RPDataReduction") << "T2 neither has tracks nor doesn't according to my IF-statements";
	  if (whichPot==12)
	    throw cms::Exception("RPDataReduction") << "RP has no tracks, but also 1 track,  according to my IF-statements";

#ifdef FOdebug    
	  std::cout<<"DEBUG: Next fill sumClusters"<<std::endl;
#endif

	  if (readClusters)
	    TwoRPTrkSumCluRest[index]->Fill(sumClusters);

	
	  allTrkCollimatorXY[index]->Fill(crudeCollimX,crudeCollimY);



#ifdef FOdebug    
	  std::cout<<"DEBUG: Next do T1T2 rapidity gap & MADX proton studies"<<std::endl;
#endif

	  
	  if (readRecoProt&&readT1)
	    {
	      numRecoPInSD->Fill(Tvals.size());
	      double rgXi=12;
	      double rgTrk=50.;
	      double gapXi=-50.;
	      bool noInelastic=false;
	      
	      if ((T1T2maxTrk<-10)||(T1T2minTrk>10))
		noInelastic=true;
	      
	      if (!noInelastic)
		{
		  if (rpEta>0)
		    {
		      rgXi=-exp(-fabs(T1T2maxTrk-rpEta));
		      rgTrk=fabs(T1T2maxTrk-rpEta);
		    }
		  if (rpEta<0)
		    {
		      rgXi=-exp(-fabs(T1T2minTrk-rpEta));
		      rgTrk=fabs(T1T2minTrk-rpEta);
		    }
		  /*
		    if (rpEta>0)
		    rgXi=-exp(-fabs((T2maxTrk>-10 ? T2maxTrk : T2maxEta)-rpEta));
		    if (rpEta<0)
		    rgXi=-exp(-fabs((T2minTrk<10 ? T2minTrk : T2minEta)-rpEta));
		  */
		  
		}
	      double absXi;
	      double relErr;
	      bool kompatible=false;
	      bool Ycut=(fabs(TwoRpY0[0])>9)&&(fabs(TwoRpY0[0])<30);
	      DMapType::const_iterator pTee=Tvals.begin();
#ifdef FOdebug
          std::cout<<"DEBUG: Next, iterate over (t,xi) doublet"<<std::endl;
#endif
	      
	      if (!Tvals.empty())
		{

		  for (;pTee!=Tvals.end();pTee++)
		    {
		      if (rpEta!=0)
			{

			  if (!noInelastic)
			    {
			      xiRpVsXiRG[index]->Fill(-(pTee->second),-rgXi);
			      if ((index==kT2oppo)&&(readLoNeg&&((trigStatus/64)%2)))
				{
				  //				  if (diffTeleArm) //T1 also on p-side
				  if (bothT1Trk)
				    {
				      if (Ycut)
					xiRGT2oppo[kT1SameT2Oppo]->Fill(-(pTee->second),-rgXi);

				      if (-(pTee->second)>0)
					gapXi=fabs(log(-pTee->second));
				      RGTrkVsRGxiT1BothT2Oppo->Fill(rgTrk-gapXi); //rgTrk-gapXi;
				      if (Ycut)
					RGTrk2DRGxiT1BothT2Oppo->Fill(rgTrk,gapXi); //T1 & T2 Tracks RapGap vs ln(-xi)
				      if (fabs(rgTrk-gapXi)<2*sigmaRg)
					kompatible=true;
//				      if (fabs(angleX+0.0001)<0.00005)
//					{
//					  std::cout<<"EVENT AROUND ANGLE_X==-1E-4 :: angleX: "<<angleX<<" , x_nr: "<<TwoRpX0[0];
//					  std::cout<<" , xi: "<<(pTee->second);
//					  std::cout<<"  event_num: "<<dataT.event_num<<"  trigNum: "<<dataT.trigger_num<<"  orbit: ";
//					  std::cout<<dataT.orbit_num<<"  runNum: "<<dataT.run_num<<std::endl;
//					}

				      double sumT1Trk=T1PlusMulti+T1MinusMulti;
				      double sumT2Trk=T2PlusMulti+T2MinusMulti;

				      noYcutXiRGT2oppo[kT1SameT2Oppo][whichPair]->Fill(-(pTee->second),-rgXi);
				      if (Ycut)
					YcutXiRGT2oppo[kT1SameT2Oppo][whichPair]->Fill(-(pTee->second),-rgXi);

				      if (pTee->second>0.048)
					{
					  FiverAngleXvsX->Fill(angleX,TwoRpX0[0]);
					  FiverXvsY->Fill(TwoRpX0[0],TwoRpY0[0]);
					  
					  if (ThetaYs.size()==1)
					    FiverThetaYMineVsMadx->Fill(ipAngleY,ThetaYs[0]);
					  else
					    std::cout<<"Not exactly one recoProton for case 2/12 RP trk - size(Theta Y)="<<ThetaYs.size()<<std::endl;

					  if (ThetaXes.size()==1)
					    FiverThetaXMineVsMadx->Fill(ipAngleX,ThetaXes[0]);
					  else
					    std::cout<<"Not exactly one recoProton for case 2/12 RP trk - size(Theta X)="<<ThetaXes.size()<<std::endl;
					}
				      if (kompatible)
					{
					  okRgAngleXvsX->Fill(angleX,TwoRpX0[0]);
					  T1T2numTrk[0]->Fill(sumT1Trk);
					  T1T2numTrk[2]->Fill(sumT2Trk);
					  rgCollimatorXY[0]->Fill(crudeCollimX,crudeCollimY);
					  RgXiOkrpTee[0][0]->Fill(-pTee->first);
					}
				      else
					{
					  badRgAngleXvsX->Fill(angleX,TwoRpX0[0]);
					  T1T2numTrk[1]->Fill(sumT1Trk);
					  T1T2numTrk[3]->Fill(sumT2Trk);
					  rgCollimatorXY[1]->Fill(crudeCollimX,crudeCollimY);
					  RgXiOkrpTee[1][0]->Fill(-pTee->first);
					}
					
				   

				    }
#ifdef FOdebug
          std::cout<<"DEBUG: Next, same-arm & no-T1 conditions"<<std::endl;
#endif


				  if (sameTeleArm)
				    noYcutXiRGT2oppo[kT1T2Oppo][whichPair]->Fill(-(pTee->second),-rgXi);
					  
				  if (noT1Trk)
				    noYcutXiRGT2oppo[kNoT1][whichPair]->Fill(-(pTee->second),-rgXi);
				  
				  if (Ycut)
				    {
				      
				      if (sameTeleArm)
					YcutXiRGT2oppo[kT1T2Oppo][whichPair]->Fill(-(pTee->second),-rgXi);
					  
				      if (noT1Trk)
					YcutXiRGT2oppo[kNoT1][whichPair]->Fill(-(pTee->second),-rgXi);
				    }
				  
				    

				  if (sameTeleArm) //T1 & T2 both not on p-side
				    xiRGT2oppo[kT1T2Oppo]->Fill(-(pTee->second),-rgXi);
				  if (noT1Trk) 
				    xiRGT2oppo[kNoT1]->Fill(-(pTee->second),-rgXi);
				}
#ifdef FOdebug
          std::cout<<"DEBUG: Next, T2 both-arms condition"<<std::endl;
#endif


			      if (index==kT2both)
				{

#ifdef FOdebug
          std::cout<<"DEBUG: Next, fill gapXi"<<std::endl;
#endif


				  if (-(pTee->second)>0)
				    gapXi=fabs(log(-pTee->second));
				  if (fabs(rgTrk-gapXi)<2*sigmaRg)
				    kompatible=true;
				  else
				    kompatible=false;


				  if (bothT1Trk) 
				    {
#ifdef FOdebug
          std::cout<<"DEBUG: Next, fill RgXiOkrpTee"<<std::endl;
#endif

				       if (readLoNeg)
				          RgXiOkrpTee[kompatible ? 0 : 1][1]->Fill(-pTee->first);
#ifdef FOdebug
          std::cout<<"DEBUG: Next, fill xiRgT2Both & xiBadOkRGT2both[0]"<<std::endl;
#endif


				      xiRGT2both[kT1both]->Fill(-(pTee->second),-rgXi);
				      xiBadOkRGT2both[0][whichPair]->Fill(-(pTee->second),-rgXi);

				      if (Ycut)
					{
#ifdef FOdebug
          std::cout<<"DEBUG: Next, Ycut+fill xiRgT2Both & xiBadOkRGT2both[0]"<<std::endl;
#endif


					  xiBadOkRGT2both[1][whichPair]->Fill(-(pTee->second),-rgXi);
					  RGTrk2DRGxiT1T2Both->Fill(rgTrk,gapXi); 

					}
					
				    }
#ifdef FOdebug
          std::cout<<"DEBUG: Next, fill eitherT1 or noT1 histos"<<std::endl;
#endif


				  if (eitherT1Trk) //T1 trk only in one arm
				    xiRGT2both[kT1either]->Fill(-(pTee->second),-rgXi);
				  if (noT1Trk) 
				    xiRGT2both[kT1neither]->Fill(-(pTee->second),-rgXi);
				}


			    }
			  if (whichVerb>10)
			    std::cout<<"xi:"<<rgXi<<" , T2 min "<<T2minEta<<" , T2max:"<<T2maxEta<<" , rp eta="<<rpEta<<std::endl;
			}

#ifdef FOdebug
          std::cout<<"DEBUG: Next, all-pot xi & t-histograms filled"<<std::endl;
#endif


		      rpTee[index]->Fill((pTee->first)*(-1.));
		      absXi=fabs(pTee->second);
		      if (pTee->first)
			relErr=100.-100.*calcTeeLoXi/(pTee->first);
		      else
			relErr=0;
		      funcTeeRelErrVsXi[index]->Fill(pTee->second,relErr);
		      funcTeeVsElastTee[index]->Fill(-(pTee->first),-calcTeeLoXi);
		      
			  
			

		      
		      if (absXi<0.005)
			{
			  if (ThetaXes.size()==1)
			    {
			      if (ThetaXes[0])
				ipAngleXRelErr[index]->Fill(100.-100.*ipAngleX/ThetaXes[0]);
			    }
			  
			  if (ThetaYs.size()==1)
			    {
			      if (ThetaYs[0])
				ipAngleYRelErr[index]->Fill(100.-100.*ipAngleY/ThetaYs[0]);
			    }
			}
		      /*
			if (absXi<0.03)
			{
			if (whichVerb>10)
			std::cout<<"recoT: "<<(-pTee->first)<<"calcTee: "<<calcTeeLoXi<<std::endl;
			}
			else
			{
			if (whichVerb>10)
			std::cout<<"HI XI,recoT: "<<(-pTee->first)<<"calcTee: "<<calcTeeLoXi<<std::endl;
			
			}
		      */
		      if (readLoNeg&&((trigStatus/64)%2))
			{

			  if (ThetaXes.size()==1)
			    T2ipAngleX[index]->Fill(ThetaXes[0]);
			  if (ThetaYs.size()==1)
			    T2ipAngleY[index]->Fill(ThetaYs[0]);

			  t2rpTee[index]->Fill((pTee->first)*(-1.));
			  t2rpTeePair[index][whichPair]->Fill(-(pTee->first));
			  t2rpTeeVsXi[index]->Fill((pTee->first)*(-1.),pTee->second);
			}
		      
		    }
		}

	    }
	  if (readLoNeg&&readClusters)
	    {

	      bool rpVtrigged=(trigStatus%2);
	      int sameClu=(rp45Trk ? sumClu45 : sumClu56);
	      int oppClu=(rp45Trk ? sumClu56 : sumClu45);
	      int oppMaxU=(rp45Trk ? maxU56 : maxU45);
	      int oppMaxV=(rp45Trk ? maxV56 : maxV45);
	      int oppMaxUV= ((oppMaxV>oppMaxU) ? oppMaxV : oppMaxU);
	      //	      sumClu45,rp45Trk;
	      float same=(sameClu ? sameClu : 0.1);
	      float other=(oppClu ? -oppClu : -0.1);
	      SDCluNoTrk[index]->Fill(same);
	      SDCluNoTrk[index]->Fill(other);
	    
	      otherSideCluVsMaxTMC[index]->Fill(oppMaxUV,-other);
	      for (int a=0;a<6;a++)
		{
		  otherSideCluVsCluPerPot[index]->Fill(HzAllPotClu[rp56Trk ? a : a+6],-other);
		  if (bunchNumSD>0)
		    BNon0OtherSideClu[index]->Fill(HzAllPotClu[rp56Trk ? a : a+6],-other);
		}

	      if (rpVtrigged)
		{
		  rpTrgSDCluNoTrk[index]->Fill(same);
		  rpTrgSDCluNoTrk[index]->Fill(other);
		  RPTrigOtherSideCluVsMaxTMC[index]->Fill(oppMaxUV,-other);
		  for (int a=0;a<6;a++)
		    {
		      RPTrigOtherSideCluVsCluPerPot[index]->Fill(HzAllPotClu[rp56Trk ? a : a+6],-other);
		      if (bunchNumSD>0)
			BNon0RPTrigOtherSideClu[index]->Fill(HzAllPotClu[rp56Trk ? a : a+6],-other);
		    }

		}
	    }

	  for (int u=0;u<2;u++)
	    {
	      SDtrkX[index]->Fill(TwoRpX0[u]);
	      SDtrkY[index]->Fill(TwoRpY0[u]);
	      SDXvsY[index][rp45Trk ? u : (u+2)]->Fill(TwoRpX0[u],TwoRpY0[u]);
	    }
	  
	  for (int u=0;u<12;u++)
	    {
	      if (HzAllPotSTrk[u])
		SDTrkWhichPot[index]->Fill(potsAllHz[u]);
	    }

	  SDXNearVsFar[index][whichPot]->Fill(TwoRpX0[1],TwoRpX0[0]);
	  SDYNearVsFar[index][whichPot]->Fill(TwoRpY0[1],TwoRpY0[0]);
	  SDXNearMinusFar[index][whichPot]->Fill(TwoRpX0[0]-TwoRpX0[1]);
	  SDYNearMinusFar[index][whichPot]->Fill(TwoRpY0[0]-TwoRpY0[1]);
	  SD2pAngleXvsX[index]->Fill(angleX,TwoRpX0[0]);
	  SD2pAngleYvsY[index]->Fill(angleY,TwoRpY0[0]);
	  
	  if (readLoNeg)
	    {
	      if ((trigStatus/64)%2) 
		{ //T2 triggered
		  bool rpTrg=trigStatus%4; // rp-v = %2 & rp-h = (/2)%2

		  if ((fabs(angleX)<(2*cutAngleX))&&(fabs(TwoRpX0[0])<(2*cutX)))
		    elasticBox=true;
		  if (TwoRpX0[0]>(-2*cutX))
		    sdBox=true;

		  T2TrigSD2pAngleXvsX[index]->Fill(angleX,TwoRpX0[0]); 
		  if (bunchNumSD==bigBunch)
		    {
		      T2TrigSD2pAngleXvsXBB[index]->Fill(angleX,TwoRpX0[0]); 
		      T2SpectroQuarterBB[index][whichPair]->Fill(angleX,TwoRpX0[0]);
		  
		    }
		  else
		    {
		      T2TrigSD2pAngleXvsXOB[index]->Fill(angleX,TwoRpX0[0]); 
		      T2SpectroQuarterOB[index][whichPair]->Fill(angleX,TwoRpX0[0]);
		    }

		  T2SpectroQuarter[index][whichPair]->Fill(angleX,TwoRpX0[0]);

		  
		  if (!rpTrg)
		    T2noRPTrigSD2pAngleXvsX[index]->Fill(angleX,TwoRpX0[0]); 
		}
	    }


#ifdef FOdebug    
  std::cout<<"DEBUG: Next do box-cutted MADX proton studies"<<std::endl;
#endif



	  if (readRecoProt&&(!Tvals.empty()))
	    {
	      DMapType::const_iterator pTee=Tvals.begin();
	      for (;pTee!=Tvals.end();pTee++)
		{
		  if (readLoNeg&&((trigStatus/64)%2))
		    {
		      if (elasticBox)
			{
			  ElastT2rpTee[index]->Fill((pTee->first)*(-1.));
			  ElastT2rpTeeVsXi[index]->Fill((pTee->first)*(-1.),pTee->second);
			  if (ThetaXes.size()==1)
			    {
			      double absXi=fabs(pTee->second);
			      ElaT2ipAngleX[index]->Fill(ThetaXes[0]);
			      if (ThetaXes[0]&&(absXi<0.005))
				ElaIpAngleXRelErr[index]->Fill(100.-100.*ipAngleX/ThetaXes[0]);
			    }
			  if (ThetaYs.size()==1)
			    {
			      double absXi=fabs(pTee->second);

			      ElaT2ipAngleY[index]->Fill(ThetaYs[0]);
			      if (ThetaYs[0]&&(absXi<0.005))
				ElaIpAngleYRelErr[index]->Fill(100.-100.*ipAngleY/ThetaYs[0]);

			    }

			}
		      if (sdBox)
			{
			  SdT2rpTee[index]->Fill((pTee->first)*(-1.));
			  SDt2rpTeePair[index][whichPair]->Fill(-(pTee->first));
			  if (ThetaXes.size()==1)
			    SdT2ipAngleX[index]->Fill(ThetaXes[0]);
			  if (ThetaYs.size()==1)
			    SdT2ipAngleY[index]->Fill(ThetaYs[0]);

			}
		    }

		}
	    }

	}




#ifdef FOdebug    
  std::cout<<"DEBUG: Next, case with 1/12 pots w/ a track"<<std::endl;
#endif





      if (badRp45Trk||badRp56Trk)
	{
	  int indexB=kT2NA;
	  if (bothArmsBad)
	    {
	      badRpT2TrBoth[falseBunch]++;

	      if (!badRpT2TrBoth.count(bunchNumSD))
		{
		  badRpT2TrBoth[bunchNumSD]=1;
		}
	      else
		{
		  badRpT2TrBoth[bunchNumSD]++;
		}


	      indexB=kT2both;
	    }
      
	  if (sameArmBad)
	    {
	      badRpT2TrSame[falseBunch]++;


	      if (!badRpT2TrSame.count(bunchNumSD))
		{
		  badRpT2TrSame[bunchNumSD]=1;
		}
	      else
		{
		  badRpT2TrSame[bunchNumSD]++;
		}



	      indexB=kT2same;
	    }
      
	  if (opposArmBad)
	    {
	      badRpT2TrOppos[falseBunch]++;

	      if (!badRpT2TrOppos.count(bunchNumSD))
		{
		  badRpT2TrOppos[bunchNumSD]=1;
		}
	      else
		{
		  badRpT2TrOppos[bunchNumSD]++;
		}

	      indexB=kT2oppo;
	    }

	  if (rpOnlyBad)
	    {
	      rpBadTrT2NoTr[falseBunch]++;

	      if (!rpBadTrT2NoTr.count(bunchNumSD))
		{
		  rpBadTrT2NoTr[bunchNumSD]=1;
		}
	      else
		{
		  rpBadTrT2NoTr[bunchNumSD]++;
		}

	      indexB=kT2none;
	    }
	  for (int u=0;u<12;u++)
	    {
	      if (indexB==kT2NA)
		throw cms::Exception("RPDataReduction") << "T2 neither has tracks nor doesn't according to my IF-statements";

	      if (HzAllPotSTrk[u])
		SDBadTrkWhichPot[indexB]->Fill(potsAllHz[u]);
	    }

	}
    }

#ifdef FOdebug    
  std::cout<<"DEBUG: Next cluster hit quantities"<<std::endl;
#endif



  if (readClusters)
    {
      for (int u=0;u<64;u++) 
	partialEvts[(u/32)%2][(u/8)%4][u%8]+=partialHits[(u/32)%2][(u/8)%4][u%8];
      for (int u=0;u<1280;u++) 
	perPlaneEvts[u]+=planeHit[u];

      //    std::cout<<"num of RP clusters in evt: "<<totalClu<<std::endl;

      for (int u=0;u<64;u++) 
	{
	  if (totRPHits[(u/32)%2][(u/8)%4][u%8]) prfNClu->Fill(((u/32)%2)*100+((u/8)%4)*10+(u%8),totRPHits[(u/32)%2][(u/8)%4][u%8]);
	}
  
      for (int r=0;r<1280;r++)
	{
	  //      OneEvtCluSum+=occupF[r];
	  if (occupF[r])
	    prfOccup->Fill(r,occupF[r]/512.);
	}
      
      for (int u=0;u<128;u++)
	{
      
	  unsigned int UTr=0; unsigned int VTr=0;
	  for (int g=0;g<10;g++)
	    {
	      if (planeEffHit[10*u+g]>nMaxPrPl)
		planeEffHit[10*u+g]=0;
	      if (planeEffHit[10*u+g])
		{
		  if (g%2)
		    {		  // this is U-plane
		      UTr++;
		    }
		  else
		    {
		      VTr++;
		    }

		}
	
	    }

	  if ((UTr>=nPl)&&(VTr>=nPl))
	    {
	      evtsEff[u]++;
	      if (noTr[u])
		noTrkEff[u]++;
	    }
	  int occupPot=0;
	  for (int w=0;w<10;w++)
	    occupPot+=occupF[u*10+w];
	  if ((!multiT[u])&&(!siTr[u])&& ((u%100)<26)&&((u%100)>19)&&occupPot)
	    noTrkSomeClus[u]++;
	}


      prfMultipl->Fill(OneEvtCluSum,OneEvtTrkSum);

      for (int t=0;t<6;t++)
	{
	  prfMul[t]->Fill(OnePotCluSum[20+t],OnePotTrkSum[20+t]);
	  prfMul[t+6]->Fill(OnePotCluSum[120+t],OnePotTrkSum[120+t]);
	}
    }


}    


void RPDataReduction::beginRun(edm::Run const&, edm::EventSetup const& es)
{

  edm::ESHandle<BeamOpticsParams> BOParH;
  es.get<BeamOpticsParamsRcd>().get(BOParH);
  if(!BOParH.isValid())
    throw cms::Exception("RPDataReduction::beginRun") << " edm::ESHandle<BeamOpticsParams> is invalid";
  optObj=*BOParH;
 

}

// ------------ method called once each job just before starting event loop  ------------
void 
RPDataReduction::beginJob()
{

  //  const unsigned int potsD2[4]={20,24,121,125};

  for (int u=0;u<128;u++) 
    {
      singleTOnly[u]=0;
      noTrkEff[u]=0;
      evtsEff[u]=0;
      noTrk[u]=0;
      noTrkSomeClus[u]=0;
      multiTrk[u]=0;
    }

  for (int u=0;u<64;u++) 
    totalClu[(u/32)%2][(u/8)%4][u%8]=0;
  for (int u=0;u<64;u++) 
    partialEvts[(u/32)%2][(u/8)%4][u%8]=0;
  for (int u=0;u<64;u++) 
    cluSize[(u/32)%2][(u/8)%4][u%8]=0;
  for (int u=0;u<64;u++) 
    SOccup[(u/32)%2][(u/8)%4][u%8]=0.;
  for (int u=0;u<1280;u++)
    perPlaneEvts[u]=0;
  totalEvts=0; refEvtsOn=0; checkEvtsOn=0; refEvtsOff=0; checkEvtsOff=0;noChkTrkOn=0; noChkTrkOff=0;


  potAll4On[falseBunch]=0; potAll4Off[falseBunch]=0; pots3On1Off[falseBunch]=0; pots3On1Half[falseBunch]=0;
  potAll4Trk[falseBunch]=0;   potTrk3S1M[falseBunch]=0; potTrk2S2M[falseBunch]=0;   pots3On1HalfAll[falseBunch]=0;

  refEvts=0;
  pot3Plus1Medium[falseBunch]=0;
  pot3Plus1Max[falseBunch]=0; 

  pot2Plus2Medium[falseBunch]=0; 

  pot2T2Empty[falseBunch]=0; refEvtsVale=0;
  checkEvts=0; 


  TwoOnTwoHalfAll[falseBunch]=0 ; TwoOn2ATwoHalfAll[falseBunch]=0;
  TwoOnTwoOff[falseBunch]=0; TwoOnTwoHalf[falseBunch]=0; TwoOnTwoVMul[falseBunch]=0; TwoOn2ATwoOff[falseBunch]=0; TwoOn2ATwoHalf[falseBunch]=0; TwoOn2ATwoVMul[falseBunch]=0;
  MTrCand01AnyPot=0; TrCand01AnyPot=0; MTrC01_diag31=0; TrC01_diag31=0; MTrC01_diag22=0; TrC01_diag22=0;
  AnyPot01=0; MTr01AnyPot=0; STr01AnyPot=0; mptyEvts=0; 

  pot2A2Plus2Medium[falseBunch]=0; 

  // not used yet :  pot2On2HHv=0; pot2On2HHu=0;

  potTrk3S1NoClu[falseBunch]=0; potTrk2S2NoClu[falseBunch]=0;
  pot2ATrk2S2NoClu[falseBunch]=0; pot3On1HHv[falseBunch]=0; pot3On1HHu[falseBunch]=0; pot2ATrk2S2M[falseBunch]=0;
  //not used yet: TwoDiag00=0; TwoDiagUorV=0;
  AllPot00[falseBunch]=0; AllPotUorV[falseBunch]=0;  OneDiag00SomepotOn[falseBunch]=0;
  // not used yet:  rpSTrT2NoHit=0;
  
  rpSTrT2NoTr[falseBunch]=0;
  rpT2TrOppos[falseBunch]=0; rpT2TrSame[falseBunch]=0;
  rpT2TrBoth[falseBunch]=0;
  rpSTrT2NoTrFC[falseBunch]=0;
  rpT2TrOpposFC[falseBunch]=0; rpT2TrSameFC[falseBunch]=0;
  rpT2TrBothFC[falseBunch]=0;
  rpSTrT2NoTrT2T[falseBunch]=0;
  rpT2TrOpposT2T[falseBunch]=0; rpT2TrSameT2T[falseBunch]=0;
  rpT2TrBothT2T[falseBunch]=0;

  rpTmcT2NoTr[falseBunch]=0;
  rpTmcT2NoTrT[falseBunch]=0;
  rpTmcT2NoTrRT[falseBunch]=0;
  rpTmcT2TrSame[falseBunch]=0;
  rpTmcT2TrSameT[falseBunch]=0;
  rpTmcT2TrSameRT[falseBunch]=0;
  rpTmcT2TrOppo[falseBunch]=0;
  rpTmcT2TrOppoT[falseBunch]=0;
  rpTmcT2TrOppoRT[falseBunch]=0;
  rpTmcT2TrBoth[falseBunch]=0;
  rpTmcT2TrBothT[falseBunch]=0;
  rpTmcT2TrBothRT[falseBunch]=0;


  badRpT2TrSame[falseBunch]=0;
  badRpT2TrOppos[falseBunch]=0;
  badRpT2TrBoth[falseBunch]=0;
  rpBadTrT2NoTr[falseBunch]=0;
  trig220cross[falseBunch]=0;  trig220h[falseBunch]=0; trig220v[falseBunch]=0; trigT2[falseBunch]=0; trigT1[falseBunch]=0; trigBx[falseBunch]=0;
  potDiag3STrk=0;
  potDiag2STrk=0;
  potDiag1STrk=0;
  potDiag0STrk=0;

  for (int k=0;k<6;k++)
    {
      rp2T2n[k].clear();
      rp2T2s[k].clear();
      rp2T2o[k].clear();
      rp2T2b[k].clear();
    }

 //  for (int u=0;u<4;u++)
  //    {
  //      potUOn[u]=0;
  //      potUOff[u]=0;
  //    }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
RPDataReduction::endJob() 
{
  //  const unsigned int potsD2[4]={20,24,121,125};
  int iA=-1;int iS=-1;int iRP=-1;int totTot=0;


#ifdef FOdebug    
  std::cout<<"DEBUG: Begin endjob"<<std::endl;
#endif


  std::cout<<std::endl<<std::endl<<"Total events: "<<totalEvts<<std::endl<<std::endl;
  if (readClusters)
    {
      std::cout<<std::endl<<std::endl<<"CLUSTERING RESULTS:"<<std::endl<<"---------------------"<<std::endl;
      for (int u=0;u<2*4*8;u++) 
	{
	  iRP=u%8;
	  iS=(u/8)%4;
	  iA=(u/32)%2;
	  if (totalClu[iA][iS][iRP])
	    {
	      std::cout<<std::endl;
	      std::cout<<"Arm:"<<iA<<" Station:"<<iS<<" RP:"<<iRP;
	      std::cout<<" average number of clusters: "<<((float) totalClu[iA][iS][iRP])/totalEvts;
	      std::cout<<" avg clu size: "<<((float) cluSize[iA][iS][iRP])/totalClu[iA][iS][iRP]<<std::endl;
	      std::cout<<" avg occupancy (all evts): "<<( SOccup[iA][iS][iRP]*100.)/totalEvts<<" %"<<std::endl;
	      std::cout<<" avg occupancy (evts where this pot not empty): "<<( SOccup[iA][iS][iRP]*100.)/partialEvts[iA][iS][iRP]<<"%"<<std::endl;
	      
	    }
	  totTot+=totalClu[iA][iS][iRP];
	}
      std::cout<<"total clusters: "<<totTot<<std::endl;
    }
  std::cout<<std::endl<<std::endl<<"TRIGGER EFFICIENCY RESULTS:"<<std::endl<<"---------------------"<<std::endl;
  std::cout<<"Reference pot : "<<refPot<<" and checked pot: "<<checkPot<<" (opposite triggering pot: "<<oppositePot<<")"<<std::endl<<std::endl;
  double fracErrOn=1.;
  double fracErr=1.;
  double fracErrOff=1.;
  double fracErrMpty=1.;

  // Approximately, the binomial error
  //  if ((checkEvts>50)&&(refEvts>50))
  if ((checkEvtsOn>10)&&(refEvtsOn>10))
    {
      fracErrOn=sqrt(1-checkEvtsOn/refEvtsOn);
      fracErrOn=sqrt(checkEvtsOn)*fracErrOn/refEvtsOn;
    }
  if ((checkEvts>10)&&(refEvts>10))
    {
      fracErr=sqrt(1-checkEvts/refEvts);
      fracErr=sqrt(checkEvts)*fracErr/refEvts;
    }

  if ((mptyEvts>10)&&(refEvts>10))
    {
      fracErrMpty=sqrt(1-mptyEvts/refEvts);
      fracErrMpty=sqrt(mptyEvts)*fracErrMpty/refEvts;
    }
  //mptyEvts

  if ((checkEvtsOff>10)&&(refEvtsOff>10))
    {
      fracErrOff=sqrt(1-checkEvtsOff/refEvtsOff);
      fracErrOff=sqrt(checkEvtsOff)*fracErrOff/refEvtsOff;
    }

  /*
  else 
    {
      if ((checkEvts>10)&&(refEvts>10))
	{
	  fracErr=sqrt(refEvts*(refEvts+3+checkEvts)+checkEvts+2);
	  fracErr=sqrt(checkEvts)*fracErr/(refEvts*refEvts);
 	  
	}
    }
  */
  //  std::cout<<"Checked pot efficiency (%): "<<(100.*checkEvts)/refEvts<<"+-"<<100*fracErr<<" ; out of "<<refEvts<<" pointing tracks, "<<checkEvts<<" have >0 TMC bit(s) on in U and V in checkPot (arm 56 near pots have n+7.4mm cut, n="<<(elCut56-7.4)<<"), 45 nr tp:"<<elCut45t<<" and nr bt:"<<elCut45b<<std::endl;

  //  std::cout<<"Checked pot efficiency (%): "<<(100.*checkEvts)/refEvts<<"+-"<<100*fracErr<<" ; out of "<<refEvts<<" pointing elastic (tmc u=v=12, abs(extrapolCheck minus ref)>32 strips) tracks, "<<checkEvts<<" have >0 TMC bit(s) on in U and V (sectors 11-13 only considered) in checkPot (arm 56 near pots have n+7.4mm cut, n="<<(elCut56-7.4)<<"), 45 nr tp:"<<elCut45t<<" and nr bt:"<<elCut45b<<std::endl;

//  std::cout<<"Checked pot combined rate for trigger efficiency or being totally empty (%): "<<(100.*checkEvts)/refEvts<<"+-"<<100*fracErr<<" ; out of "<<refEvts<<" pointing tracks (excluding those in the cut region), "<<checkEvts<<" have >0 TMC bit(s) on in U or V in checkPot OR 0 clusters in checkPot (arm 56 near pots have n+7.4mm cut, n="<<(elCut56-7.4)<<"), 45 nr tp:"<<elCut45t<<" and nr bt:"<<elCut45b<<std::endl;
//  std::cout<<"Checked pot combined rate for trigger efficiency or being totally empty (%): "<<(100.*checkEvts)/refEvts<<"+-"<<100*fracErr<<" ; out of "<<refEvts<<" pointing tracks (excluding those in the cut region, and demanding refU=V=12), "<<checkEvts<<" have >0 TMC bit(s) on in U"<<((modeF==1) ? " and " : " or ")<<"V (#11-13) in checkPot OR 0 clusters in checkPot (arm 56 near pots have n+7.4mm cut, n="<<(elCut56-7.4)<<"), 45 nr tp:"<<elCut45t<<" and nr bt:"<<elCut45b<<std::endl;
  std::cout<<"Checked pot trigger efficiency (%): "<<(100.*checkEvts)/refEvts<<"+-"<<100*fracErr<<" ; out of "<<refEvts<<" pointing tracks (excluding those in the cut region, and demanding refU=V=12 ; 11...13/U+1=V for 56_tp), "<<checkEvts<<" have >0 TMC bit(s) on in U"<<((modeF==1) ? " and " : " or ")<<"V (#refU/V +- 1) in checkPot (arm 56 near pots have n+7.4mm cut, n="<<(elCut56-7.4)<<"), 45 nr tp:"<<elCut45t<<" and nr bt:"<<elCut45b<<std::endl;


  std::string whatCombi[9]={ " (45tp*56bt) "," (45bt*56tp) ", " (45tp*56tp) "," (45hr*56hr) ", " (45bt*56bt) ", " (45tp*56hr) "," (45bt*56hr) ", " (45hr*56tp) "," (45hr*56bt) "};


  std::cout<<std::endl<<std::endl<<"TRACKING EFFICIENCY RESULTS:"<<std::endl<<"---------------------"<<std::endl;

  std::cout<<"Checked pot tracking efficiency, TMC ON (%): "<<(100.*checkEvtsOn)/refEvtsOn<<"+-"<<100*fracErrOn<<" ; out of "<<refEvtsOn<<" pointing elastic tracks, "<<checkEvtsOn<<" have a track within 2 sigma in X & 3 sigma in Y of the extrapolated refTrack "<<"(56_tp nondiag, others check U,V=12+-1) (arm 56 near pots have n+7.4mm cut, n="<<(elCut56-7.4)<<"), 45 nr tp:"<<elCut45t<<" and nr bt:"<<elCut45b<<std::endl<<std::endl;
  std::cout<<"Checked pot tracking efficiency, TMC OFF (%): "<<(100.*checkEvtsOff)/refEvtsOff<<"+-"<<100*fracErrOff<<" ; out of "<<refEvtsOff<<" pointing elastic tracks, "<<checkEvtsOff<<" have a track within 2 sigma in X & 3 sigma in Y of the extrapolated refTrack "<<"(56_tp nondiag, others check U,V=12+-1 - found all FALSE) (arm 56 near pots have n+7.4mm cut, n="<<(elCut56-7.4)<<"), 45 nr tp:"<<elCut45t<<" and nr bt:"<<elCut45b<<std::endl<<std::endl;

  std::cout<<"Rate of no-track-in-checkPot for refEvents: "<<(100.*noChkTrkOn)/refEvtsOn<<" (% ; TMC ON) & "<<(100.*noChkTrkOff)/refEvtsOff<<" (% ; TMC OFF)"<<std::endl<<std::endl;

  if (readClusters)
    {
      std::cout<<std::endl<<std::endl<<"CLUSTERS RESULTS:"<<std::endl<<"---------------------"<<std::endl;
      std::cout<<"15jun11: Checked pot no clusters (%): "<<(100.*mptyEvts)/refEvts<<"+-"<<100*fracErrMpty<<" ; out of "<<refEvts<<" pointing tracks (excluding those in the cut region, and demanding refU=V=12 ; 11...13/U+1=V for 56_tp), "<<mptyEvts<<" have 0 clusters in checkPot (arm 56 near pots have n+7.4mm cut, n="<<(elCut56-7.4)<<"), 45 nr tp:"<<elCut45t<<" and nr bt:"<<elCut45b<<std::endl;
    }
  //  std::cout<<"Fraction of tracks outside the 3sigmaY*2sigmaX window, in one coordinate only (but still within a 5Y*3X window) (%): "<<(100.*checkEvts)/refEvts<<"+-"<<100*fracErr<<" ; out of "<<refEvts<<" events with track in both refPot and checkPot, pointing from ref 2 chk (excluding the cut region), "<<checkEvts<<" have a difference in extrapolated X vs measured checkPot X of 2...3 sigma (or 3...5 sigma for Y) (arm 56 near pots have n+7.4mm cut, n="<<(elCut56-7.4)<<"), 45 nr tp:"<<elCut45t<<" and nr bt:"<<elCut45b<<std::endl;
  //,,,pots3On1Half;

  std::cout<<std::endl<<std::endl<<"4-POT COMBINATION RESULTS:"<<std::endl<<"---------------------"<<std::endl;
  std::cout<<"number of reference events used for the below sums: "<<refEvtsVale<<std::endl;
  for (int t=0;t<64;t++)
    outputRefEvts->Fill(t,refEvtsVale);

  const char *lablsA[]={"2/12S+10/12(<10clu)!T2","2/12S+10/12(<10clu)T2,sam","2/12S+10/12(<10clu)T2,opp","2/12S+10/12(<10clu)T2,all","2/12S+!T2Trk(T2 trg)","2/12S+T2TrkSame,T2Trig","2/12S+T2TrkOpp,T2Trig","2/12S+T2TrkBoth,T2Trig","2/12Tmc!T2","2/12Tmc,T2same","2/12Tmc,T2oppo","2/12Tmc,T2both","TrigT2+Tmc+!T2","TrigT2+Tmc+T2same","TrigT2+Tmc+T2Oppo","TrigT2+Tmc+T2Both","TrigRPT2+Tmc+!T2","TrigRPT2+Tmc+T2same","TrigRPT2+Tmc+T2Oppo","TrigRPT2+Tmc+T2Both","","","","","","","","12*[00]","<12*[00]","","2/12S+T2!trk","2/12S+T2Tr,sam","2/12S+T2Tr,opp","2/12S+T2Tr,all","1/12S+T2Tr,sam","1/12S+T2Tr,opp","1/12S+T2Tr,all","1/12S+T2!trk","trig220X","trig220H","trig220V","trigT1","trigT2","trigBx","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""};
  const char *lablsD[]={"4(11)","4(00)","3(11)*1(00)","3(11)*1(01)","3(11)*1(01,m!s)","3(11)*1(2,2)","2(11)*2(2,2)","2(11)*2(2,2)/arm","3(11)*1(>4,>4)","2(11)*2(01+m)","2(11)*2(01)","2(11)*2(00)","2(11)*2(>4,>4)","2(11)*2(01+m)/arm","2(11)*2(01)/arm","2(11)*2(00)/arm","2(11)*2(44)/arm","3(11)*1(01=patt)","4*[00]+<8*[00]","4/4S","3S1M","2S2M","2S2(<1TS)/arm","3S1noClu","2S2noClu","2S2nClu/arm","2S2M/arm","","","","","","","","","","","","","","","","","","3/4S","2/4S","1/4S","0/4S","","","","","","","","","","","","","","","","","","","","","","","","","",""};
  int foundBin=0;

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, fill histos with many different combinations"<<std::endl;
#endif

  ofstream bunchNums(bunchText.c_str());
  if (bunchNums.fail())
    std::cout<<"RPDataReduction ERROR: Could not open output text file!"<<std::endl;

  if (!readMultiTrk)
    {
      lablsD[20]="";
      lablsD[21]="";
    }

  for (int r=0;r<64;r++)
    {
      foundBin=outputFulfillEvtsA->FindBin((float) r);
      outputFulfillEvtsA->GetXaxis()->SetBinLabel(foundBin,lablsA[r]);
      foundBin=outputFulfillEvtsD->FindBin((float) r);
      outputFulfillEvtsD->GetXaxis()->SetBinLabel(foundBin,lablsD[r]);
    }

  outputFulfillEvtsD->Fill(0.,potAll4On[falseBunch]);
  outputFulfillEvtsD->Fill(1.,potAll4Off[falseBunch]);
  outputFulfillEvtsD->Fill(2.,pots3On1Off[falseBunch]);
  outputFulfillEvtsD->Fill(3.,pots3On1HalfAll[falseBunch]);
  outputFulfillEvtsD->Fill(4.,pots3On1Half[falseBunch]);
  outputFulfillEvtsD->Fill(5.,pot3Plus1Medium[falseBunch]);
  outputFulfillEvtsD->Fill(6.,pot2Plus2Medium[falseBunch]);
  outputFulfillEvtsD->Fill(7.,pot2A2Plus2Medium[falseBunch]);
  outputFulfillEvtsD->Fill(8.,pot3Plus1Max[falseBunch]);
  outputFulfillEvtsD->Fill(9.,TwoOnTwoHalf[falseBunch]);
  outputFulfillEvtsD->Fill(10.,TwoOnTwoHalfAll[falseBunch]);
  outputFulfillEvtsD->Fill(11.,TwoOnTwoOff[falseBunch]);
  outputFulfillEvtsD->Fill(12.,TwoOnTwoVMul[falseBunch]);
  outputFulfillEvtsD->Fill(13.,TwoOn2ATwoHalf[falseBunch]);
  outputFulfillEvtsD->Fill(14.,TwoOn2ATwoHalfAll[falseBunch]);
  outputFulfillEvtsD->Fill(15.,TwoOn2ATwoOff[falseBunch]);
  outputFulfillEvtsD->Fill(16.,TwoOn2ATwoVMul[falseBunch]);
  outputFulfillEvtsD->Fill(17.,(pot3On1HHu[falseBunch]+pot3On1HHv[falseBunch]));

  outputFulfillEvtsD->Fill(19.,potAll4Trk[falseBunch]);

 
  outputFulfillEvtsD->Fill(44.,potDiag3STrk);
  outputFulfillEvtsD->Fill(45.,potDiag2STrk);
  outputFulfillEvtsD->Fill(46.,potDiag1STrk);
  outputFulfillEvtsD->Fill(47.,potDiag0STrk);

  if (readMultiTrk)
    {
      outputFulfillEvtsD->Fill(20.,potTrk3S1M[falseBunch]);
      outputFulfillEvtsD->Fill(21.,potTrk2S2M[falseBunch]);
      outputFulfillEvtsD->Fill(26.,pot2ATrk2S2M[falseBunch]);
    }
  outputFulfillEvtsD->Fill(22.,pot2T2Empty[falseBunch]);
  outputFulfillEvtsD->Fill(23.,potTrk3S1NoClu[falseBunch]);
  outputFulfillEvtsD->Fill(24.,potTrk2S2NoClu[falseBunch]);
  outputFulfillEvtsD->Fill(25.,pot2ATrk2S2NoClu[falseBunch]);

  outputFulfillEvtsA->Fill(27.,AllPot00[falseBunch]);
  outputFulfillEvtsA->Fill(28.,AllPotUorV[falseBunch]);

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, fill some T2-dependent combinations to the above-mentioned histograms"<<std::endl;
#endif

  if (readT2)
    {
      //      outputFulfillEvts->Fill(29.,rpSTrT2NoHit);
      outputFulfillEvtsA->Fill(30.,rpSTrT2NoTr[falseBunch]);
      outputFulfillEvtsA->Fill(31.,rpT2TrSame[falseBunch]);
      outputFulfillEvtsA->Fill(32.,rpT2TrOppos[falseBunch]);
      outputFulfillEvtsA->Fill(33.,rpT2TrBoth[falseBunch]);

      outputFulfillEvtsA->Fill(4.,rpSTrT2NoTrT2T[falseBunch]);
      outputFulfillEvtsA->Fill(5.,rpT2TrSameT2T[falseBunch]);
      outputFulfillEvtsA->Fill(6.,rpT2TrOpposT2T[falseBunch]);
      outputFulfillEvtsA->Fill(7.,rpT2TrBothT2T[falseBunch]);
      if (readClusters)
	{
	  outputFulfillEvtsA->Fill(0.,rpSTrT2NoTrFC[falseBunch]);
	  outputFulfillEvtsA->Fill(1.,rpT2TrSameFC[falseBunch]);
	  outputFulfillEvtsA->Fill(2.,rpT2TrOpposFC[falseBunch]);
	  outputFulfillEvtsA->Fill(3.,rpT2TrBothFC[falseBunch]);
	}

      outputFulfillEvtsA->Fill(34.,badRpT2TrSame[falseBunch]);
      outputFulfillEvtsA->Fill(35.,badRpT2TrOppos[falseBunch]);
      outputFulfillEvtsA->Fill(36.,badRpT2TrBoth[falseBunch]);
      outputFulfillEvtsA->Fill(37.,rpBadTrT2NoTr[falseBunch]);
    }

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, fill some loneg-dependent combinations to the above-mentioned histograms"<<std::endl;
#endif
  if (readLoNeg)
    {
      outputFulfillEvtsA->Fill(8.,rpTmcT2NoTr[falseBunch]);
      outputFulfillEvtsA->Fill(9.,rpTmcT2TrSame[falseBunch]);
      outputFulfillEvtsA->Fill(10.,rpTmcT2TrOppo[falseBunch]);
      outputFulfillEvtsA->Fill(11.,rpTmcT2TrBoth[falseBunch]);
      outputFulfillEvtsA->Fill(12.,rpTmcT2NoTrT[falseBunch]);
      outputFulfillEvtsA->Fill(13.,rpTmcT2TrSameT[falseBunch]);
      outputFulfillEvtsA->Fill(14.,rpTmcT2TrOppoT[falseBunch]);
      outputFulfillEvtsA->Fill(15.,rpTmcT2TrBothT[falseBunch]);
      outputFulfillEvtsA->Fill(16.,rpTmcT2NoTrRT[falseBunch]);
      outputFulfillEvtsA->Fill(17.,rpTmcT2TrSameRT[falseBunch]);
      outputFulfillEvtsA->Fill(18.,rpTmcT2TrOppoRT[falseBunch]);
      outputFulfillEvtsA->Fill(19.,rpTmcT2TrBothRT[falseBunch]);

      outputFulfillEvtsA->Fill(38.,trig220cross[falseBunch]);
      outputFulfillEvtsA->Fill(39.,trig220h[falseBunch]);
      outputFulfillEvtsA->Fill(40.,trig220v[falseBunch]);
      outputFulfillEvtsA->Fill(41.,trigT1[falseBunch]);
      outputFulfillEvtsA->Fill(42.,trigT2[falseBunch]);
      outputFulfillEvtsA->Fill(43.,trigBx[falseBunch]);
     }
  outputFulfillEvtsD->Fill(18.,OneDiag00SomepotOn[falseBunch]);
 
  bunchNums<<"ALL BUNCHES AGGREGATED ("<<refEvtsVale<<" events): 4-pot diagonal combinations, %-ages\n ------------------------\n\n";
  //  outputPercentages=outputFulfillEvts;
  double dPercent,aPercent;
  for (int r=1;r<65;r++)
    {
      dPercent=100.*outputFulfillEvtsD->GetBinContent(r)/outputRefEvts->GetBinContent(r);
      outputPercentagesD->SetBinContent(r,dPercent);
      outputPercentagesD->GetXaxis()->SetBinLabel(r,lablsD[r-1]);
      bunchNums<<lablsD[r-1]<<"\t"<<dPercent<<"\n";
    }

  bunchNums<<"\n\nALL BUNCHES AGGREGATED, 8/12-pot non-diagonal combinations, %-ages\n ------------------------\n\n";

  for (int r=1;r<65;r++)
    {
      aPercent=100.*outputFulfillEvtsA->GetBinContent(r)/outputRefEvts->GetBinContent(r);
      outputPercentagesA->SetBinContent(r,aPercent);
      outputPercentagesA->GetXaxis()->SetBinLabel(r,lablsA[r-1]);
      bunchNums<<lablsA[r-1]<<"\t"<<aPercent<<"\n";
   }

  if (readLoNeg)
    {
#ifdef FOdebug    
  std::cout<<"DEBUG: Next, fill some combinations to per-bunch histograms"<<std::endl;
#endif

      for (unsigned int x=0;x<10000;x++)
	{
	  outputPerBunchA[x]=std::auto_ptr<TH1D>(0);
	  outputPerBunchD[x]=std::auto_ptr<TH1D>(0);
	}
     for (unsigned int x=0;x<10000;x++)
       {
	 percentagesPerBunchA[x]=std::auto_ptr<TH1D>(0);
	 percentagesPerBunchD[x]=std::auto_ptr<TH1D>(0);
       }

      int r=0;
      std::stringstream bunch;
      std::string bunchS,name,title,nameP,nameA,titleA,nameAP,titleAP,titleP,title2;


      //      const char *lablsBunch[]={"evts","4/4S","3/4S","2/4S","1/4S","0/4S","","","","","","","","","","","","","","","4(11)","4(00)","3(11)*1(00)","3(11)*1(01)","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""};
      //  int foundBin=0;


      const char *lablsBunchA[]={"2/12S+10/12(<10clu)!T2","2/12S+10/12(<10clu)T2,sam","2/12S+10/12(<10clu)T2,opp","2/12S+10/12(<10clu)T2,all","2/12S+!T2,T2Trig","2/12S+T2Same,T2Trig","2/12S+T2Opp,T2Trig","2/12S+T2Both,T2Trig","2/12Tmc!T2","2/12Tmc,T2same","2/12Tmc,T2oppo","2/12Tmc,T2both","TrigT2+Tmc+!T2","TrigT2+Tmc+T2same","TrigT2+Tmc+T2Oppo","TrigT2+Tmc+T2Both","TrigRPT2+Tmc+!T2","TrigRPT2+Tmc+T2same","TrigRPT2+Tmc+T2Oppo","TrigRPT2+Tmc+T2Both","","","","","","","","12*[00]","<12*[00]","","2/12S+T2!trk","2/12S+T2Tr,sam","2/12S+T2Tr,opp","2/12S+T2Tr,all","1/12S+T2Tr,sam","1/12S+T2Tr,opp","1/12S+T2Tr,all","1/12S+T2!trk","trig220X","trig220H","trig220V","trigT1","trigT2","trigBx","","","","","bunEvts","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""};
      const char *lablsBunchD[]={"4(11)","4(00)","3(11)*1(00)","3(11)*1(01)","3(11)*1(01,m!s)","3(11)*1(2,2)","2(11)*2(2,2)","2(11)*2(2,2)/arm","3(11)*1(>4,>4)","2(11)*2(01+m)","2(11)*2(01)","2(11)*2(00)","2(11)*2(>4,>4)","2(11)*2(01+m)/arm","2(11)*2(01)/arm","2(11)*2(00)/arm","2(11)*2(44)/arm","3(11)*1(01=patt)","4*[00]+<8*[00]","4/4S","3S1M","2S2M","2S2(<1TS)/arm","3S1noClu","2S2noClu","2S2nClu/arm","2S2M/arm","","","","","","","","","","","","","","","","","","3/4S","2/4S","1/4S","0/4S","bunEvts","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""};
     
      for (MapType::const_iterator evBun=refEvtsBunches.begin();evBun!=refEvtsBunches.end();evBun++)
	{
	  const int rateBins=128;
	  name="outputEvtsInBunchDiag";
	  nameA="outputEvtsInBunchAll";
	  nameP="CombinationPercentagesForBunchDiag";
	  nameAP="CombinationPercentagesForBunchAll";
	  titleA="for nondiagonal combinations (reference class: all events), for bunch number ";
	  title="TMC & Track diagonal";
	  title2="combinations, for bunch number ";
 	  bunch<<(evBun->first);
	  bunchS=bunch.str();
	  name.append(bunchS);
	  nameP.append(bunchS);
	  nameA.append(bunchS);
	  nameAP.append(bunchS);
	  title2.append(bunchS);
	  titleA.append(bunchS);
	  titleAP=titleA;
	  titleAP.append(";combi;%");
	  titleA.append(";combi;evts");
	  title.append(whatCombi[diagonal]);
	  title.append(title2);
	  titleP=title;
	  titleP.append(";combi;%");
	  title.append(";combi;evts");
	  bunch.str("");
	  if ((r+1)>10000)
	    throw cms::Exception("RPDataReduction") << "Track & TMC combinations per bunch : found more than 10 000 bunches!";

	  //	  std::string

	  std::string potsP[]={"45tp,","45bt,","45hr,","56tp,","56bt,","56hr,"};
	  std::string T2P[]={"!T2","oT2","sT2","bT2"};
	  std::string rpP="2rp"; //{"2rp,","1rp,"};
	    //	  std::string temp="";
	  std::string Pbins[24];
	  int Pindex=0;

	  //	  for (int k=0;k<2;k++)
	  //	    {
	      for (int L=0;L<6;L++)
		{
		  for (int m=0;m<4;m++)
		    {
		      Pbins[Pindex].clear();
		      Pbins[Pindex].append(rpP);
		      Pbins[Pindex].append(potsP[L]);
		      Pbins[Pindex].append(T2P[m]);
		      Pindex++;
		      
		    }	      
		  
		}	      
	      //	    }

	  outputPerBunchA[r]=std::auto_ptr<TH1D>(new TH1D(nameA.c_str(),titleA.c_str(),rateBins,-0.5,rateBins-0.5));
	  percentagesPerBunchA[r]=std::auto_ptr<TH1D>(new TH1D(nameAP.c_str(),titleAP.c_str(),rateBins,-0.5,rateBins-0.5));
	  outputPerBunchD[r]=std::auto_ptr<TH1D>(new TH1D(name.c_str(),title.c_str(),rateBins,-0.5,rateBins-0.5));
	  percentagesPerBunchD[r]=std::auto_ptr<TH1D>(new TH1D(nameP.c_str(),titleP.c_str(),rateBins,-0.5,rateBins-0.5));

	  for (int rr=1;rr<(rateBins+1);rr++)
	    {
	       outputPerBunchA[r]->GetXaxis()->SetBinLabel(rr,lablsBunchA[rr-1]);
	       outputPerBunchD[r]->GetXaxis()->SetBinLabel(rr,lablsBunchD[rr-1]);
	       if (rr<49)
		 {
		   percentagesPerBunchA[r]->GetXaxis()->SetBinLabel(rr,lablsBunchA[rr-1]);
		   percentagesPerBunchD[r]->GetXaxis()->SetBinLabel(rr,lablsBunchD[rr-1]);
		 }
	    }

	  for (int rr=50;rr<(50+24);rr++)
	    {
	      outputPerBunchA[r]->GetXaxis()->SetBinLabel(rr,Pbins[rr-50].c_str());
	      percentagesPerBunchA[r]->GetXaxis()->SetBinLabel(rr,Pbins[rr-50].c_str());
	    }

	  outputPerBunchA[r]->Fill(48.,evBun->second);	  
	  outputPerBunchD[r]->Fill(48.,evBun->second);	  
	  
	  if (potAll4Trk.count(evBun->first))
	    outputPerBunchD[r]->Fill(19.,potAll4Trk[evBun->first]);	
  	    
	  
	  if (readMultiTrk)
	    {
	      if (potTrk3S1M.count(evBun->first))
		outputPerBunchD[r]->Fill(20.,potTrk3S1M[evBun->first]);
	      if (potTrk2S2M.count(evBun->first))
		outputPerBunchD[r]->Fill(21.,potTrk2S2M[evBun->first]);
	      if (pot2ATrk2S2M.count(evBun->first))
		outputPerBunchD[r]->Fill(26.,pot2ATrk2S2M[evBun->first]);
	    
	    }
	    




	  if (mp3Trk.count(evBun->first))
	    outputPerBunchD[r]->Fill(44.,mp3Trk[evBun->first]);	
  	    
	  if (mp2Trk.count(evBun->first))
	    outputPerBunchD[r]->Fill(45.,mp2Trk[evBun->first]);	
  	    
	  if (mp1Trk.count(evBun->first))
	    outputPerBunchD[r]->Fill(46.,mp1Trk[evBun->first]);	
  	    
	  if (mp0Trk.count(evBun->first))
	    outputPerBunchD[r]->Fill(47.,mp0Trk[evBun->first]);	
  	    
	  if (potAll4On.count(evBun->first))
	    outputPerBunchD[r]->Fill(0.,potAll4On[evBun->first]);	
  	    
	  if (potAll4Off.count(evBun->first))
	    outputPerBunchD[r]->Fill(1.,potAll4Off[evBun->first]);	
  	    
	  if (pots3On1Off.count(evBun->first))
	    outputPerBunchD[r]->Fill(2.,pots3On1Off[evBun->first]);	
  	    
	  if (pots3On1HalfAll.count(evBun->first))
	    outputPerBunchD[r]->Fill(3.,pots3On1HalfAll[evBun->first]);	
  	    
	  if (pots3On1Half.count(evBun->first))
	    outputPerBunchD[r]->Fill(4.,pots3On1Half[evBun->first]);	

	  if (pot3Plus1Medium.count(evBun->first))
	    outputPerBunchD[r]->Fill(5.,pot3Plus1Medium[evBun->first]);	

	  if (pot2Plus2Medium.count(evBun->first))
	    outputPerBunchD[r]->Fill(6.,pot2Plus2Medium[evBun->first]);	

	  if (pot2A2Plus2Medium.count(evBun->first))
	    outputPerBunchD[r]->Fill(7.,pot2A2Plus2Medium[evBun->first]);	


	  if (pot3Plus1Max.count(evBun->first))
	    outputPerBunchD[r]->Fill(8.,pot3Plus1Max[evBun->first]);	

 
	  if (TwoOnTwoHalf.count(evBun->first))
	    outputPerBunchD[r]->Fill(9.,TwoOnTwoHalf[evBun->first]);	

	  if (TwoOnTwoHalfAll.count(evBun->first))
	    outputPerBunchD[r]->Fill(10.,TwoOnTwoHalfAll[evBun->first]);	

	  if (TwoOnTwoOff.count(evBun->first))
	    outputPerBunchD[r]->Fill(11.,TwoOnTwoOff[evBun->first]);	

	  if (TwoOnTwoVMul.count(evBun->first))
	    outputPerBunchD[r]->Fill(12.,TwoOnTwoVMul[evBun->first]);	

	  if (TwoOn2ATwoHalf.count(evBun->first))
	    outputPerBunchD[r]->Fill(13.,TwoOn2ATwoHalf[evBun->first]);	

	  if (TwoOn2ATwoHalfAll.count(evBun->first))
	    outputPerBunchD[r]->Fill(14.,TwoOn2ATwoHalfAll[evBun->first]);	

	  if (TwoOn2ATwoOff.count(evBun->first))
	    outputPerBunchD[r]->Fill(15.,TwoOn2ATwoOff[evBun->first]);	

	  if (TwoOn2ATwoVMul.count(evBun->first))
	    outputPerBunchD[r]->Fill(16.,TwoOn2ATwoVMul[evBun->first]);	

 
	  if (pot3On1HHu.count(evBun->first))
	    outputPerBunchD[r]->Fill(17.,pot3On1HHu[evBun->first]);	

	  if (pot3On1HHv.count(evBun->first))
	    outputPerBunchD[r]->Fill(17.,pot3On1HHv[evBun->first]);	

	  if (OneDiag00SomepotOn.count(evBun->first))
	    outputPerBunchD[r]->Fill(18.,OneDiag00SomepotOn[evBun->first]);

	  if (pot2T2Empty.count(evBun->first))
	    outputPerBunchD[r]->Fill(22.,pot2T2Empty[evBun->first]);

	  if (potTrk3S1NoClu.count(evBun->first))
	    outputPerBunchD[r]->Fill(23.,potTrk3S1NoClu[evBun->first]);

	  if (potTrk2S2NoClu.count(evBun->first))
	    outputPerBunchD[r]->Fill(24.,potTrk2S2NoClu[evBun->first]);

	  if (pot2ATrk2S2NoClu.count(evBun->first))
	    outputPerBunchD[r]->Fill(25.,pot2ATrk2S2NoClu[evBun->first]);

	  if (AllPot00.count(evBun->first))
	    outputPerBunchA[r]->Fill(27.,AllPot00[evBun->first]);

	  if (AllPotUorV.count(evBun->first))
	    outputPerBunchA[r]->Fill(28.,AllPotUorV[evBun->first]);

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, fill some T2-dependent combinations to the above-mentioned histograms #2"<<std::endl;
#endif
 

	  if (readT2)
	    {


	      if (readLoNeg)
		{
		  if (rpTmcT2NoTr.count(evBun->first))
		    outputPerBunchA[r]->Fill(8.,rpTmcT2NoTr[evBun->first]);

		  if (rpTmcT2TrSame.count(evBun->first))
		    outputPerBunchA[r]->Fill(9.,rpTmcT2TrSame[evBun->first]);

		  if (rpTmcT2TrOppo.count(evBun->first))
		    outputPerBunchA[r]->Fill(10.,rpTmcT2TrOppo[evBun->first]);

		  if (rpTmcT2TrBoth.count(evBun->first))
		    outputPerBunchA[r]->Fill(11.,rpTmcT2TrBoth[evBun->first]);

		  if (rpTmcT2NoTrT.count(evBun->first))
		    outputPerBunchA[r]->Fill(12.,rpTmcT2NoTrT[evBun->first]);

		  if (rpTmcT2TrSameT.count(evBun->first))
		    outputPerBunchA[r]->Fill(13.,rpTmcT2TrSameT[evBun->first]);

		  if (rpTmcT2TrOppoT.count(evBun->first))
		    outputPerBunchA[r]->Fill(14.,rpTmcT2TrOppoT[evBun->first]);

		  if (rpTmcT2TrBothT.count(evBun->first))
		    outputPerBunchA[r]->Fill(15.,rpTmcT2TrBothT[evBun->first]);

		  if (rpTmcT2NoTrRT.count(evBun->first))
		    outputPerBunchA[r]->Fill(16.,rpTmcT2NoTrRT[evBun->first]);

		  if (rpTmcT2TrSameRT.count(evBun->first))
		    outputPerBunchA[r]->Fill(17.,rpTmcT2TrSameRT[evBun->first]);

		  if (rpTmcT2TrOppoRT.count(evBun->first))
		    outputPerBunchA[r]->Fill(18.,rpTmcT2TrOppoRT[evBun->first]);

		  if (rpTmcT2TrBothRT.count(evBun->first))
		    outputPerBunchA[r]->Fill(19.,rpTmcT2TrBothRT[evBun->first]);
		}

	      if (rpSTrT2NoTr.count(evBun->first))
		outputPerBunchA[r]->Fill(30.,rpSTrT2NoTr[evBun->first]);

	      if (rpT2TrSame.count(evBun->first))
		outputPerBunchA[r]->Fill(31.,rpT2TrSame[evBun->first]);

	      if (rpT2TrOppos.count(evBun->first))
		outputPerBunchA[r]->Fill(32.,rpT2TrOppos[evBun->first]);

	      if (rpT2TrBoth.count(evBun->first))
		outputPerBunchA[r]->Fill(33.,rpT2TrBoth[evBun->first]);


	      if (rpSTrT2NoTrT2T.count(evBun->first))
		outputPerBunchA[r]->Fill(4.,rpSTrT2NoTrT2T[evBun->first]);

	      if (rpT2TrSameT2T.count(evBun->first))
		outputPerBunchA[r]->Fill(5.,rpT2TrSameT2T[evBun->first]);

	      if (rpT2TrOpposT2T.count(evBun->first))
		outputPerBunchA[r]->Fill(6.,rpT2TrOpposT2T[evBun->first]);

	      if (rpT2TrBothT2T.count(evBun->first))
		outputPerBunchA[r]->Fill(7.,rpT2TrBothT2T[evBun->first]);

	      if (readClusters)
		{
		  if (rpSTrT2NoTrFC.count(evBun->first))
		    outputPerBunchA[r]->Fill(0.,rpSTrT2NoTrFC[evBun->first]);
		  
		  if (rpT2TrSameFC.count(evBun->first))
		    outputPerBunchA[r]->Fill(1.,rpT2TrSameFC[evBun->first]);

		  if (rpT2TrOpposFC.count(evBun->first))
		    outputPerBunchA[r]->Fill(2.,rpT2TrOpposFC[evBun->first]);

		  if (rpT2TrBothFC.count(evBun->first))
		    outputPerBunchA[r]->Fill(3.,rpT2TrBothFC[evBun->first]);

		}

	      if (badRpT2TrSame.count(evBun->first))
		outputPerBunchA[r]->Fill(34.,badRpT2TrSame[evBun->first]);

	      if (badRpT2TrOppos.count(evBun->first))
		outputPerBunchA[r]->Fill(35.,badRpT2TrOppos[evBun->first]);

	      if (badRpT2TrBoth.count(evBun->first))
		outputPerBunchA[r]->Fill(36.,badRpT2TrBoth[evBun->first]);

	      if (rpBadTrT2NoTr.count(evBun->first))
		outputPerBunchA[r]->Fill(37.,rpBadTrT2NoTr[evBun->first]);

	      for (int k=0;k<6;k++)
		{
		  if (rp2T2n[k].count(evBun->first))
		    outputPerBunchA[r]->Fill(49.+4.*k,rp2T2n[k][evBun->first]);

		  if (rp2T2o[k].count(evBun->first))
		    outputPerBunchA[r]->Fill(50.+4.*k,rp2T2o[k][evBun->first]);

		  if (rp2T2s[k].count(evBun->first))
		    outputPerBunchA[r]->Fill(51.+4.*k,rp2T2s[k][evBun->first]);

		  if (rp2T2b[k].count(evBun->first))
		    outputPerBunchA[r]->Fill(52.+4.*k,rp2T2b[k][evBun->first]);

		}

	    }

	  if (trig220cross.count(evBun->first))
	    outputPerBunchA[r]->Fill(38.,trig220cross[evBun->first]);
 
	  if (trig220h.count(evBun->first))
	    outputPerBunchA[r]->Fill(39.,trig220h[evBun->first]);
 
	  if (trig220v.count(evBun->first))
	    outputPerBunchA[r]->Fill(40.,trig220v[evBun->first]);
 
	  if (trigT1.count(evBun->first))
	    outputPerBunchA[r]->Fill(41.,trigT1[evBun->first]);
 
	  if (trigT2.count(evBun->first))
	    outputPerBunchA[r]->Fill(42.,trigT2[evBun->first]);
 

	  if (trigBx.count(evBun->first))
	    outputPerBunchA[r]->Fill(43.,trigBx[evBun->first]);
 
	  int foundBin2=outputPerBunchD[r]->FindBin((float) 48.);
 	  int foundBin2A=outputPerBunchA[r]->FindBin((float) 48.);

	  bunchNums<<"\n\n BUNCH NUMBER "<<bunchS<<", with "<<outputPerBunchA[r]->GetBinContent(foundBin2A)<<" events\n --------------------------";
	  
	  double aBPercent,dBPercent;
	  bunchNums<<"\n 8/12-pot non-diagonal combinations, %-ages\n ------------------------\n\n";

	  for (int rr=1;rr<(rateBins+1);rr++)
	    {

	      if (rr!=49)
		{
		  std::string hisLabl=outputPerBunchA[r]->GetXaxis()->GetBinLabel(rr);
		  aBPercent=100.*outputPerBunchA[r]->GetBinContent(rr)/outputPerBunchA[r]->GetBinContent(foundBin2A);
		  percentagesPerBunchA[r]->SetBinContent(rr,aBPercent);
		  bunchNums<<hisLabl<<"\t"<<aBPercent<<"\n";
		}
	    }

	  bunchNums<<"\n 4-pot diagonal combinations, %-ages\n ------------------------\n\n";

	  for (int rr=1;rr<(rateBins+1);rr++)
	    {
	      if (rr!=49)
		{
		
		  std::string hisLabl=outputPerBunchD[r]->GetXaxis()->GetBinLabel(rr);
		  dBPercent=100.*outputPerBunchD[r]->GetBinContent(rr)/outputPerBunchD[r]->GetBinContent(foundBin2);
		  percentagesPerBunchD[r]->SetBinContent(rr,dBPercent);
		  bunchNums<<hisLabl<<"\t"<<dBPercent<<"\n";
		}
	    }
 
  	    
	  r++;
	}
    }


  bunchNums.close();
#ifdef FOdebug    
  std::cout<<"DEBUG: Next, fill some T2-dependent LoNeg numbers"<<std::endl;
#endif

  if (readT2&&readLoNeg)
    {
      rpvBunchTrigs=std::auto_ptr<TH1D>(new TH1D("rpvBunchTrigs","LoNeg RP-V trigger bit for SD events with the given bunch colliding (2/12 RP Trk (nr+fr)+Opposite side T2 Trk)",rpvTrig.size(),-0.5,rpvTrig.size()-0.5));
      
      rphBunchTrigs=std::auto_ptr<TH1D>(new TH1D("rphBunchTrigs","LoNeg RP-H trigger bit for SD events with the given bunch colliding (2/12 RP Trk (nr+fr)+Opposite side T2 Trk)",rphTrig.size(),-0.5,rphTrig.size()-0.5));
      t2BunchTrigs=std::auto_ptr<TH1D>(new TH1D("t2BunchTrigs","LoNeg T2 trigger bit for SD events with the given bunch colliding (2/12 RP Trk (nr+fr)+Opposite side T2 Trk)",t2Trig.size(),-0.5,t2Trig.size()-0.5));
      bxBunchTrigs=std::auto_ptr<TH1D>(new TH1D("bxBunchTrigs","LoNeg BX trigger bit for SD events with the given bunch colliding (2/12 RP Trk (nr+fr)+Opposite side T2 Trk)",bxTrig.size(),-0.5,bxTrig.size()-0.5));
      sdEvBunches=std::auto_ptr<TH1D>(new TH1D("sdEvBunches","Number of SD events with the given bunch colliding (2/12 RP Trk (nr+fr)+Opposite side T2 Trk)",sdBunchEv.size(),-0.5,sdBunchEv.size()-0.5));
      
      


      int foundBinSd=0;
      
      float r=0;
      
      MapType::const_iterator rpv=rpvTrig.begin();
      MapType::const_iterator rph=rphTrig.begin();
      MapType::const_iterator t2=t2Trig.begin();
      MapType::const_iterator bx=bxTrig.begin();
      MapType::const_iterator sdEv=sdBunchEv.begin();
      std::stringstream bunch;
      std::string bunchS;
      
      for (;rpv!=rpvTrig.end();rpv++)
	{
	  
	  bunch<<(rpv->first);
	  bunchS=bunch.str();

	  foundBinSd=rpvBunchTrigs->FindBin(r);
	  rpvBunchTrigs->Fill(r,(rpv->second));
	  rpvBunchTrigs->GetXaxis()->SetBinLabel(foundBinSd,bunchS.c_str());
	  bunch.str("");
	  r++;
	}

      r=0.;
      for (;rph!=rphTrig.end();rph++)
	{
	  
	  bunch<<(rph->first);
	  bunchS=bunch.str();
	  
	  foundBinSd=rphBunchTrigs->FindBin(r);
	  rphBunchTrigs->Fill(r,(rph->second));
	  rphBunchTrigs->GetXaxis()->SetBinLabel(foundBinSd,bunchS.c_str());
	  bunch.str("");
	  r++;
	}
    
      r=0.;
      for (;t2!=t2Trig.end();t2++)
	{
	  
	  bunch<<(t2->first);
	  bunchS=bunch.str();
	  
	  foundBinSd=t2BunchTrigs->FindBin(r);
	  t2BunchTrigs->Fill(r,(t2->second));
	  t2BunchTrigs->GetXaxis()->SetBinLabel(foundBinSd,bunchS.c_str());
	  bunch.str("");
	  r++;
	}
    
      r=0.;
      for (;bx!=bxTrig.end();bx++)
	{
	  
	  bunch<<(bx->first);
	  bunchS=bunch.str();
	  
	  foundBinSd=bxBunchTrigs->FindBin(r);
	  bxBunchTrigs->Fill(r,(bx->second));
	  bxBunchTrigs->GetXaxis()->SetBinLabel(foundBinSd,bunchS.c_str());
	  bunch.str("");
	  r++;
	}
 

      r=0.;
      for (;sdEv!=sdBunchEv.end();sdEv++)
	{
	  
	  bunch<<(sdEv->first);
	  bunchS=bunch.str();
	  
	  foundBinSd=sdEvBunches->FindBin(r);
	  sdEvBunches->Fill(r,(sdEv->second));
	  sdEvBunches->GetXaxis()->SetBinLabel(foundBinSd,bunchS.c_str());
	  bunch.str("");
	  r++;
	}
 






    }


  //  outputPercentages->Scale(100.);

  //  for (int r=0;r<16;r++)
  //    {
  //      foundBin=outputPercentages->FindBin((float) r);
  //      outputPercentages->GetXaxis()->SetBinLabel(foundBin,labls[r]);
  //    }


#ifdef FOdebug    
  std::cout<<"DEBUG: Next, print numerical values"<<std::endl;
#endif


  std::cout<<std::endl<<std::endl<<"TMC RESULTS:"<<std::endl<<"---------------------"<<std::endl;
  //  std::string vilkenDiag=(diagonal ? " 45_bt*56_tp " :  " 45_tp*56_bt ");
  std::cout<<"Percentage of events where all 4 pots in combi"<<whatCombi[diagonal]<<"are U=V=on (1 TS on in each) :"<<(100.*potAll4On[falseBunch]/refEvtsVale)<<std::endl<<std::endl;
  std::cout<<"Percentage of events where all 4 pots in combi"<<whatCombi[diagonal]<<"are U=V=off (0 TS on in each) :"<<( 100.*potAll4Off[falseBunch]/refEvtsVale)<<std::endl<<std::endl;
  std::cout<<"Percentage of events where 3 pots in combi"<<whatCombi[diagonal]<<"are U=V=on (1 TS on in each) & 1 pot is off U=V=0 TS on:"<<( pots3On1Off[falseBunch]*100./refEvtsVale )<<std::endl<<std::endl;

  std::cout<<"Percentage of events where 3 pots in combi"<<whatCombi[diagonal]<<"are U=V=on (1 TS on in each) & 1 pot is on U,V=(0,1)/(1,0):"<<( 100.*pots3On1HalfAll[falseBunch]/refEvtsVale )<<std::endl<<std::endl;

  if (readMultiTrk)  
    std::cout<<"Percentage of events where 3 pots in combi"<<whatCombi[diagonal]<<"are U=V=on (1 TS on in each) & 1 pot is on U,V=(0,1)/(1,0).AND. it has a multitrack :"<<( 100.*pots3On1Half[falseBunch]/refEvtsVale )<<std::endl<<std::endl;

  std::cout<<"Percentage of events where 3 pots in combi"<<whatCombi[diagonal]<<"are U=V=on (1 TS on in each) & 1 pot is crowded; U&V belong [2,4] :"<<( 100.*pot3Plus1Medium[falseBunch]/refEvtsVale )<<std::endl<<std::endl;
  std::cout<<"Percentage of events where 2 pots in combi"<<whatCombi[diagonal]<<"are U=V=on (1 TS on in each) & 2 pots are crowded; U&V belong [2,4] :"<<( 100.*pot2Plus2Medium[falseBunch]/refEvtsVale )<<std::endl<<std::endl;
  std::cout<<"Percentage of events where 2 pots in combi"<<whatCombi[diagonal]<<"are U=V=on (1 TS on in each), 1 per arm & the 2 other pots are crowded; U&V belong [2,4] :"<<( 100.*pot2A2Plus2Medium[falseBunch]/refEvtsVale )<<std::endl<<std::endl;
  std::cout<<"Percentage of events where 3 pots in combi"<<whatCombi[diagonal]<<"are U=V=on (1 TS on in each) & 1 pot is very crowded; U&V>4 :"<<( 100.*pot3Plus1Max[falseBunch]/refEvtsVale )<<std::endl<<std::endl;


  std::cout<<"Percentage of events where 4 pots in combi"<<whatCombi[diagonal]<<"are U=V=off (0 TS on in each) & at least 1 other pot is on (U or V >0, includes horizontals):"<<( 100.*OneDiag00SomepotOn[falseBunch]/refEvtsVale )<<std::endl<<std::endl;

  std::cout<<"Percentage of events where all 12 pots are U=V=off (0 TS on in each) (includes horizontals):"<<( 100.*AllPot00[falseBunch]/refEvtsVale )<<std::endl<<std::endl;

  std::cout<<"Percentage of events where at least one out of 12 pots is (U or V=on; >0 TS on in each) (includes horizontals):"<<( 100.*AllPotUorV[falseBunch]/refEvtsVale )<<std::endl<<std::endl;

  if (readMultiTrk)
    std::cout<<"Percentage of events with (U,V)=(1,1) in 2 pots & (U,V)=(0,1) or (1,0).AND.a multitrack in both the other 2 pots, in combi"<<whatCombi[diagonal]<<": "<<((100.*TwoOnTwoHalf[falseBunch])/refEvtsVale)<<std::endl<<std::endl;

  std::cout<<"Percentage of events with (U,V)=(1,1) in 2 pots & (U,V)=(0,1) or (1,0) in the other 2 pots, in combi"<<whatCombi[diagonal]<<": "<<((100.*TwoOnTwoHalfAll[falseBunch])/refEvtsVale)<<std::endl<<std::endl;
  std::cout<<"Percentage of events with (U,V)=(1,1) in 2 pots & (U,V)=(0,0) in the other 2 pots, in combi"<<whatCombi[diagonal]<<": "<<((100.*TwoOnTwoOff[falseBunch])/refEvtsVale)<<std::endl<<std::endl;
  std::cout<<"Percentage of events with (U,V)=(1,1) in 2 pots & (U,V)=(>4,>4) in the other 2 pots, in combi"<<whatCombi[diagonal]<<": "<<((100.*TwoOnTwoVMul[falseBunch])/refEvtsVale)<<std::endl<<std::endl<<std::endl<<std::endl;
  std::cout<<"same as the 3 numbers above, but demanding the two pots with (U,V)=(1,1) be in different arms"<<std::endl<<"-----------------"<<std::endl;

  if (readMultiTrk)
    std::cout<<"Percentage of events with (U,V)=(1,1) in 2 pots, 1 per arm & (U,V)=(0,1) or (1,0).AND.a multitrack in both the other 2 pots, in combi"<<whatCombi[diagonal]<<": "<<((100.*TwoOn2ATwoHalf[falseBunch])/refEvtsVale)<<std::endl<<std::endl;

  std::cout<<"Percentage of events with (U,V)=(1,1) in 2 pots, 1 per arm & (U,V)=(0,1) or (1,0) in the other 2 pots, in combi"<<whatCombi[diagonal]<<": "<<((100.*TwoOn2ATwoHalfAll[falseBunch])/refEvtsVale)<<std::endl<<std::endl;
  std::cout<<"Percentage of events with (U,V)=(1,1) in 2 pots, 1 per arm & (U,V)=(0,0) in the other 2 pots, in combi"<<whatCombi[diagonal]<<": "<<((100.*TwoOn2ATwoOff[falseBunch])/refEvtsVale)<<std::endl<<std::endl;
  std::cout<<"Percentage of events with (U,V)=(1,1) in 2 pots, 1 per arm & (U,V)=(>4,>4) in the other 2 pots, in combi"<<whatCombi[diagonal]<<": "<<((100.*TwoOn2ATwoVMul[falseBunch])/refEvtsVale)<<std::endl<<std::endl;
  std::cout<<"Percentage of events with (U,V)=(1,1) in 3 pots & (U,V)=(01/10).AND.recognized (nonPara singleTrk) pattern in U or V only in the other pot, in combi"<<whatCombi[diagonal]<<": "<<((100.*(pot3On1HHu[falseBunch]+pot3On1HHv[falseBunch]))/refEvtsVale)<<std::endl<<std::endl;
  //  std::cout<<"Percentage of events with (U,V)=(1,1) in 3 pots & (U,V)=(01/10).AND.recognized pattern in V only in the other pot, in combi"<<whatCombi[diagonal]<<": "<<((100.*pot3On1HHv[falseBunch])/refEvtsVale)<<std::endl;



  std::cout<<std::endl<<std::endl<<"TRACK RESULTS:"<<std::endl<<"---------------------"<<std::endl;

  std::cout<<"Percentage of events with SingleTrack in all 4 pots in combi"<<whatCombi[diagonal]<<": "<<((100.*potAll4Trk[falseBunch])/refEvtsVale)<<std::endl<<std::endl;

  std::cout<<"Percentage of events with SingleTrack in 3 pots in combi"<<whatCombi[diagonal]<<": "<<((100.*potDiag3STrk)/refEvtsVale)<<std::endl<<std::endl;

  std::cout<<"Percentage of events with SingleTrack in 2 pots in combi"<<whatCombi[diagonal]<<": "<<((100.*potDiag2STrk)/refEvtsVale)<<std::endl<<std::endl;

  std::cout<<"Percentage of events with SingleTrack in 1 pot in combi"<<whatCombi[diagonal]<<": "<<((100.*potDiag1STrk)/refEvtsVale)<<std::endl<<std::endl;

  std::cout<<"Percentage of events with no SingleTrack in any pot in combi"<<whatCombi[diagonal]<<": "<<((100.*potDiag0STrk)/refEvtsVale)<<std::endl<<std::endl;
  
  if (readMultiTrk)
    std::cout<<"Percentage of events with SingleTrack in 3 pots & only MultiTrk in the other 1, in combi"<<whatCombi[diagonal]<<": "<<((100.*potTrk3S1M[falseBunch])/refEvtsVale)<<std::endl<<std::endl;

  if (readMultiTrk)
    std::cout<<"Percentage of events with SingleTrack in 2 pots & only MultiTrk in the other 2, in combi"<<whatCombi[diagonal]<<": "<<((100.*potTrk2S2M[falseBunch])/refEvtsVale)<<std::endl<<std::endl;

  if (readMultiTrk)
    std::cout<<"Percentage of events with SingleTrack in 2 pots & only MultiTrk in the other 2 (1 per arm), in combi"<<whatCombi[diagonal]<<": "<<((100.*pot2ATrk2S2M[falseBunch])/refEvtsVale)<<std::endl<<std::endl;
  std::cout<<"Percentage of events with SingleTrack in 2 pots, 1 per arm & max 1 TS on in total in the other 2 pots, in combi"<<whatCombi[diagonal]<<": "<<((100.*pot2T2Empty[falseBunch])/refEvtsVale)<<std::endl<<std::endl;

  if (readClusters)
    {
      std::cout<<std::endl<<std::endl<<"TRACK+CLUSTERS RESULTS:"<<std::endl<<"---------------------"<<std::endl;
      std::cout<<"Percentage of events with SingleTrack in 3 pots & no clusters in the other one, in combi"<<whatCombi[diagonal]<<": "<<((100.*potTrk3S1NoClu[falseBunch])/refEvtsVale)<<std::endl<<std::endl;
      std::cout<<"Percentage of events with SingleTrack in 2 pots & no clusters in the other 2, in combi"<<whatCombi[diagonal]<<": "<<((100.*potTrk2S2NoClu[falseBunch])/refEvtsVale)<<std::endl<<std::endl;
      std::cout<<"Percentage of events with SingleTrack in 2 pots, 1 per arm & no clusters in the other 2, in combi"<<whatCombi[diagonal]<<": "<<((100.*pot2ATrk2S2NoClu[falseBunch])/refEvtsVale)<<std::endl<<std::endl;
    }
  //  pot3Plus1Max,pot3Plus1Medium,STr01AnyPot,MTr01AnyPot

  //  =0; =0; TwoOnTwoVMul=0; TwoOn2ATwoOff=0; TwoOn2ATwoHalf=0; TwoOn2ATwoVMul=0;


  std::cout<<std::endl<<std::endl<<"SINGLE POT RESULTS:"<<std::endl<<"---------------------"<<std::endl;
  std::cout<<"TRACKING FOR POTS with U=0 & V=1 TMC sector on, or vice versa : total number of pots with (0,1) TMC: "<<AnyPot01<<std::endl;


  if (readMultiTrk)
    std::cout<<"Percent of these with a MultiTrkCand of size >0 in same pot :"<<(100*MTrCand01AnyPot/AnyPot01)<<std::endl;

  std::cout<<"Percent of these with a (Single)TrkCand that is Fittable() in same pot :"<<(100*TrCand01AnyPot/AnyPot01)<<std::endl;

  if (readMultiTrk)
    std::cout<<"Percent of these with a MultiTrack of size >0 in same pot :"<<(100*MTr01AnyPot/AnyPot01)<<std::endl;
  std::cout<<"Percent of these with a (NonParall) SingleTrk in same pot :"<<(100*STr01AnyPot/AnyPot01)<<std::endl;

  if (readLoNeg)
    {
      std::cout<<std::endl<<std::endl<<"LONEG TRIGGER RATES:"<<std::endl<<"---------------------"<<std::endl;
      std::cout<<"RP 220 cross-triggers: "<<trig220cross[falseBunch]<<" events."<<std::endl;
      std::cout<<"RP 220 vertical triggers: "<<trig220v[falseBunch]<<" events."<<std::endl;
      std::cout<<"RP 220 horizontal triggers: "<<trig220h[falseBunch]<<" events."<<std::endl;
      std::cout<<"T2 triggers: "<<trigT2[falseBunch]<<" events."<<std::endl;
      std::cout<<"T1 triggers: "<<trigT1[falseBunch]<<" events."<<std::endl;
      std::cout<<"Bunch Crossing (bx) triggers: "<<trigBx[falseBunch]<<" events."<<std::endl;
    }


  std::cout<<std::endl<<std::endl<<"SETTINGS USED:"<<std::endl<<"---------------------"<<std::endl;
  std::cout<<"RefTrack X0()-cut WAS used: xLow= "<<exxLow<<" and xHigh= "<<exxHi<<std::endl;
  if (modeF==0)
    std::cout<<"CheckPot TMC considered on if U.OR.V have >=1 bit on in refTS+-1"<<std::endl;
  if (modeF==1)
    std::cout<<"CheckPot TMC considered on if U.AND.V have >=1 bit on in refTS+-1"<<std::endl;
  if (modeF==3)
    std::cout<<"Was U or V with V disconnected, i.e. actually only ask 'U=on?'"<<std::endl;
  
  if (modeF==4)
    std::cout<<"Was U or V with U disconnected, i.e. actually only ask 'V=on?'"<<std::endl;


  std::string outputFile=outFile;

  std::cout<<"Number of events fulfilling all these criteria (refU=refV=12, pointing refTrk, refTrk xCut,spectrometer cut, checkTrk):"<<pullYElast->GetEntries()<<std::endl;

  std::auto_ptr<TH1D> hSF=std::auto_ptr<TH1D>(new TH1D("hSF","Single trk fraction per pot (%)",128,-0.5,127.5));
  std::auto_ptr<TH1D> hZF=std::auto_ptr<TH1D>(new TH1D("hZF","No-trk fraction per pot (%)",128,-0.5,127.5));
  std::auto_ptr<TH1D> hZCF;
  std::auto_ptr<TH1D> hTrkIneff;
  if (readClusters)
    {
      hZCF=std::auto_ptr<TH1D>(new TH1D("hZCF","No-trk+some-clusters, fraction per pot (%)",128,-0.5,127.5));
      hTrkIneff=std::auto_ptr<TH1D>(new TH1D("hTrkIneff","No-trk+1toMax-clusters-perPl in Min 3 U and V pl, fraction per pot (%)",128,-0.5,127.5));
    }
  std::auto_ptr<TH1D> hMF=std::auto_ptr<TH1D>(new TH1D("hMF","Multi-trk fraction per pot (%)",128,-0.5,127.5));

  for (int hh=0;hh<128;hh++)
    {
      hSF->Fill(hh,singleTOnly[hh]*100./totalEvts);
      hZF->Fill(hh,noTrk[hh]*100./totalEvts);
      if (readClusters)
	{ 
	  hZCF->Fill(hh,noTrkSomeClus[hh]*100./totalEvts);
	  if (evtsEff[hh])
	    hTrkIneff->Fill(hh,noTrkEff[hh]*100./evtsEff[hh]);
	}
      hMF->Fill(hh,multiTrk[hh]*100./totalEvts);
    }
  std::auto_ptr<TH1D> hNClu,hSClu,hOccuAll,hOccuNE,hZPo,hZPl;

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, cluster/hit occupancy"<<std::endl;
#endif


  if (readClusters)
    {
      hNClu=std::auto_ptr<TH1D>(new TH1D("hNClu","Avg num of clusters",128,-0.5,127.5));
      hSClu=std::auto_ptr<TH1D>(new TH1D("hSClu","Avg cluster size",128,-0.5,127.5));
      hOccuAll=std::auto_ptr<TH1D>(new TH1D("hOccuAll","Strip occupancy (%)",128,-0.5,127.5));
      hOccuNE=std::auto_ptr<TH1D>(new TH1D("hOccuNE","Strip occupancy (NClu.gt.0 in this pot, %)",128,-0.5,127.5));
      hZPo=std::auto_ptr<TH1D>(new TH1D("hZPo","Fraction of events where pot has no clusters (%)",128,-0.5,127.5));
      hZPl=std::auto_ptr<TH1D>(new TH1D("hZPl","Fraction of events where plane has no clusters (%)",1280,-0.5,1279.5));
      
      for (int u=0;u<2*4*8;u++) 
	{
	  iRP=u%8;
	  iS=(u/8)%4;
	  iA=(u/32)%2;
	  if ((iS==2)&&(iRP<6)&&(iRP>=0))
	    {
	      hNClu->Fill(100*iA+10*iS+iRP,((float) totalClu[iA][iS][iRP])/partialEvts[iA][iS][iRP]);
	      hSClu->Fill(100*iA+10*iS+iRP,((float) cluSize[iA][iS][iRP])/totalClu[iA][iS][iRP]);
	      hOccuAll->Fill(100*iA+10*iS+iRP, (SOccup[iA][iS][iRP]*100.)/totalEvts);
	      hOccuNE->Fill(100*iA+10*iS+iRP,(SOccup[iA][iS][iRP]*100.)/partialEvts[iA][iS][iRP]);
	      hZPo->Fill(100*iA+10*iS+iRP,100.- (float) partialEvts[iA][iS][iRP]*100./totalEvts);
	    }
	}

      for (int u=0;u<1280;u++) 
	{
	  if (((u%1000)>199)&&((u%1000)<260))
	    hZPl->Fill(u,100.- (float) perPlaneEvts[u]*100./totalEvts);
	
	}
    }

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, open ROOT file for editing"<<std::endl;
#endif

  TFile *of = TFile::Open(outputFile.c_str(), "recreate");
  TDirectory *dir;

  if(!of || !of->IsWritable()||of->IsZombie())
    {
      throw cms::Exception("RPDataReduction") << "Output file not opened correctly, nothing was written!";
    }
  else
    {
#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write cluPlots"<<std::endl;
#endif

      if (readClusters)
	{
	  dir = of->mkdir("cluPlots", "Plots about clusters in RP");
	  dir->cd();
	  hNClu->Write();
	  hSClu->Write();
	  hZPo->Write();
	  hZPl->Write();
	  hTrkIneff->Write();
	  hOccuAll->Write();
	  hOccuNE->Write();
	  prfNClu->Write();
	  numHitsInCheckPointing->Write();
	  prfSClu->Write();
	  prfOccup->Write();
	  prfMultipl->Write();
	  for (int t=0;t<12;t++)
	    prfMul[t]->Write();
	  hZCF->Write();
	}

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write uvPlots"<<std::endl;
#endif

      dir = of->mkdir("uvPlots", "Plots about RP coordinates (U,V)");
      dir->cd();
      trgU->Write();
      trgV->Write();
      srcU->Write();
      srcV->Write();
      trgUexPos->Write();
      trgVexPos->Write();
      srcUpos->Write();
      srcVpos->Write();
      UvsTMCu->Write();
      UvsTMCv->Write();
      VvsTMCu->Write();
      VvsTMCv->Write();
      refUVsCheckU->Write();
      refVVsCheckV->Write();
      refTSvsExtrapU->Write();
      refTSvsCheckU->Write();
      refTSvsCheckV->Write();
      refTSvsExtrap->Write();

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write T1T2eta Plots"<<std::endl;
#endif


      dir = of->mkdir("etaPlots", "Plots about T1T2Track eta : track fit vs assumed origin at IP");
      dir->cd();
      if (readT2)
	{
	  T2FitEtaVs000Eta->Write();
	  T2Eta000Dispersion->Write();
	}
      if (readT1&&readT2)
	{
	  T1FitEtaVs000Eta->Write();
	  T1Eta000Dispersion->Write();
	}

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write xyPlots"<<std::endl;
#endif

      dir = of->mkdir("xyPlots", "Plots about RP Tracks in coordinates (x,y)");
      dir->cd();
      trkTx->Write();
      trkTy->Write();
      pullX->Write();
      pullXpY->Write();
      pullY->Write();
      pullXXc->Write();
      pullYYc->Write();
      srcY->Write();
      chkY->Write();
      srcYHit->Write();
      chkYHit->Write();
      refEvtXY->Write();
      chkEvtXY->Write();
      srcX->Write();
      pullYElast->Write();
      trkMisPointXY->Write();
#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write potPlots"<<std::endl;
#endif
 
      dir = of->mkdir("potPlots", "Plots about RP Tracks");
      dir->cd();
      hSF->Write();
      hZF->Write();
      hMF->Write();
      multiPot->Write();


#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write tmcPlots"<<std::endl;
#endif

      dir = of->mkdir("tmcPlots", "Plots about RP TMC bits");
      dir->cd();
      trksWOneTmcBitOn->Write();
      NotrksWOneTmcBitOn->Write();
      uv11Times3Plus001S->Write();
      elasticRefVsCheckY->Write();
      mul01->Write();
      num01->Write();
      sing01->Write();
      mul01cand->Write();
      sing01cand->Write();
      mul01diag3p1->Write();
      mul01diag2p2->Write();
      oneOriOnly3p1->Write();
      oneOriOnly2p2->Write();
      oneOriOnly3p011S->Write();

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write spectrometrPlots"<<std::endl;
#endif

      dir = of->mkdir("spectromPlots", "Plots about Hubert spectrometer cut");
      dir->cd();
      spectro45tp->Write();
      spectro45bt->Write();
      pullYVsSpec->Write();
      spectro56tp->Write();
      spectro56bt->Write();
      spectroElastUvsV->Write();
      spectroElastRefX->Write();
      if (readT2)
	{
	  for (int u=0;u<4;u++)
	    allTrkCollimatorXY[u]->Write();
	}

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write all various combinations, for all & per bunch"<<std::endl;
#endif

      dir = of->mkdir("combiPlotsRPT2LoNeg", "Rates for RP track & TMC diagonal combinations, nondiagonal T2+RP SD rates & LoNeg rates");
      dir->cd();
      outputPercentagesA->Write();
      outputFulfillEvtsA->Write();
      outputPercentagesD->Write();
      outputFulfillEvtsD->Write();
      outputRefEvts->Write();
      for (int t=0;t<12;t++)
	outputPerPot[t]->Write();

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write perBunchPlots"<<std::endl;
#endif

      if (readLoNeg)
	{
	  dir = of->mkdir("LoNegPlots", "LoNeg params, bunch number etc");
	  dir->cd();
	  for (int j=0;j<7;j++)
	    lonegBunch[j]->Write();

	  dir = of->mkdir("diagBunchPlots", "Per-bunch diagonal combinations");
	  dir->cd();

	  for (unsigned int x=0;x<refEvtsBunches.size();x++)
	    {
	      outputPerBunchD[x]->Write();
	      percentagesPerBunchD[x]->Write();
	    }

	  dir = of->mkdir("nonDiagBunchPlots", "Per-bunch all-events combinations");
	  dir->cd();
	  for (unsigned int x=0;x<refEvtsBunches.size();x++)
	    {
	      outputPerBunchA[x]->Write();
	      percentagesPerBunchA[x]->Write();
	    }
	    

	}

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write recoProtPlots"<<std::endl;
#endif

      if (readRecoProt)
	{
	  dir = of->mkdir("protonPlots", "Plots about reconstructed protons");
	  dir->cd();
	  recoProtXiVsN->Write();
	  recoProtXiVsErr->Write();
	  TwoRPTrkNumProt->Write();
	  numRecoPInSD->Write();
	  if (readT2)
	    {
	      for (int u=0;u<4;u++)
		{
		  rpTee[u]->Write();
		  funcTeeVsElastTee[u]->Write();
		  funcTeeRelErrVsXi[u]->Write();
		  ipAngleXRelErr[u]->Write();
		  ipAngleYRelErr[u]->Write();
		  ElaIpAngleXRelErr[u]->Write();
		  ElaIpAngleYRelErr[u]->Write();
		  //		  numRecoPInSD[u]->Write();
		}

	      if (readLoNeg)
		{
		
#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write xiPlots"<<std::endl;
#endif
		  if (readT1)
		    {
		      dir = of->mkdir("xiPlots", "Plots about reconstructed protons vs T1 & T2 tracks (rap.gap)");
		      dir->cd();
	 

		      for (int u=0;u<4;u++)
			RgXiOkrpTee[u%2][u/2]->Write();
		      okRgAngleXvsX->Write();
		      badRgAngleXvsX->Write();
		      FiverAngleXvsX->Write();
		      FiverXvsY->Write();
		      for (int u=0;u<4;u++)
			T1T2numTrk[u]->Write();
		      rgCollimatorXY[0]->Write();
		      rgCollimatorXY[1]->Write();
		      RGTrkVsRGxiT1BothT2Oppo->Write();
		      RGTrk2DRGxiT1BothT2Oppo->Write();
		      RGTrk2DRGxiT1T2Both->Write();
		      for (int u=0;u<3;u++)
			xiRGT2oppo[u]->Write();
		      for (int u=0;u<3;u++)
			xiRGT2both[u]->Write();
		      for (int u=0;u<4;u++)
			xiRpVsXiRG[u]->Write();
		      
#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write xiPerPlots"<<std::endl;
#endif

		      dir = of->mkdir("xiPlotsPerPot", "Plots about reconstructed protons vs T1 & T2 tracks, divided by RP quarters (rap.gap)");
		      dir->cd();

		      
		      for (int x=0;x<3;x++)
			{
			  for (int y=0;y<6;y++)
			    {
			      YcutXiRGT2oppo[x][y]->Write();
			      noYcutXiRGT2oppo[x][y]->Write();
			      if (x<2)
				xiBadOkRGT2both[x][y]->Write();
			    }
			}

		    }

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write protSDPlots"<<std::endl;
#endif

		  dir = of->mkdir("protonSDPlots", "Plots about reconstructed protons vs RP tracks, T2 SD-cases");
		  dir->cd();

		  FiverThetaYMineVsMadx->Write();
		  FiverThetaXMineVsMadx->Write();

		  for (int u=0;u<4;u++)
		    T2ipAngleX[u]->Write();
		  for (int u=0;u<4;u++)
		    T2ipAngleY[u]->Write();
		  for (int u=0;u<4;u++)
		    ElaT2ipAngleX[u]->Write();
		  for (int u=0;u<4;u++)
		    ElaT2ipAngleY[u]->Write();
		  for (int u=0;u<4;u++)
		    SdT2ipAngleX[u]->Write();
		  for (int u=0;u<4;u++)
		    SdT2ipAngleY[u]->Write();
		  for (int u=0;u<4;u++)
		    t2rpTee[u]->Write();
		  for (int u=0;u<4;u++)
		    t2rpTeeVsXi[u]->Write();
		  for (int u=0;u<4;u++)
		    ElastT2rpTee[u]->Write();
		  for (int u=0;u<4;u++)
		    SdT2rpTee[u]->Write();
		  for (int u=0;u<4;u++)
		    ElastT2rpTeeVsXi[u]->Write();
		  
#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write protSDPerPlotPlots"<<std::endl;
#endif
		    

		  dir = of->mkdir("protonSDPerPotPlots", "Plots about reconstructed protons vs RP tracks, T2 SD-cases & RP quarters");
		  dir->cd();

		  for (int u=0;u<4;u++)
		    {
		      for (int v=0;v<6;v++)
			{
			  t2rpTeePair[u][v]->Write();
			  SDt2rpTeePair[u][v]->Write();
			  T2SpectroQuarter[u][v]->Write();
			  T2SpectroQuarterBB[u][v]->Write();
			  T2SpectroQuarterOB[u][v]->Write();
			}
		    }

		}
	    }
	}

#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write sdPlots"<<std::endl;
#endif


      if (readT2)
	{
	  dir = of->mkdir("sdPlots", "Plots about SD: RP & T2 combination");
	  dir->cd();
	  
	  if (readLoNeg)
	    {
	      sdTrigs->Write();
	      rpvBunchTrigs->Write();
	      rphBunchTrigs->Write();
	      t2BunchTrigs->Write();
	      bxBunchTrigs->Write();
	      sdEvBunches->Write();
		
	    }

	  for (int i=0;i<4;i++)
	    {
	      if (readClusters)
		{
		  TwoRPTrkSumCluRest[i]->Write();
		  if (readLoNeg)
		    {
		      SDCluNoTrk[i]->Write();
		      rpTrgSDCluNoTrk[i]->Write();
		      RPTrigOtherSideCluVsMaxTMC[i]->Write();
		      otherSideCluVsMaxTMC[i]->Write();
		      RPTrigOtherSideCluVsCluPerPot[i]->Write();
		      otherSideCluVsCluPerPot[i]->Write();
		      BNon0RPTrigOtherSideClu[i]->Write();
		      BNon0OtherSideClu[i]->Write();
		    }
		}

	      SDtrkX[i]->Write();
	      SDtrkY[i]->Write();
	      SDT2Multi[i]->Write();

	      for (int j=0;j<4;j++)
		SDXvsY[i][j]->Write();
	      if (readLoNeg)
		{
		  T2TrigSD2pAngleXvsX[i]->Write();
		  T2TrigSD2pAngleXvsXBB[i]->Write();
		  T2TrigSD2pAngleXvsXOB[i]->Write();
		  T2noRPTrigSD2pAngleXvsX[i]->Write();
		}
	      SD2pAngleXvsX[i]->Write();
	      SD2pAngleYvsY[i]->Write();
	      SDTrkWhichPot[i]->Write();
	      SDBadTrkWhichPot[i]->Write();
	    }
	    
#ifdef FOdebug    
  std::cout<<"DEBUG: Next, write sdPerPotPlots"<<std::endl;
#endif

	  dir = of->mkdir("sdPerPotPlots", "Plots about SD: RP & T2 combination, divided by RP pot");
	  dir->cd();


	  for (int i=0;i<4;i++)
	    {
	      for (int j=0;j<6;j++)
		{
		  SDXNearVsFar[i][j]->Write();
		  SDYNearVsFar[i][j]->Write();
		  SDXNearMinusFar[i][j]->Write();
		  SDYNearMinusFar[i][j]->Write();

		}
	    }
	}

	


    }

  delete of;

}

//define this as a plug-in
DEFINE_FWK_MODULE(RPDataReduction);
