//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>

//STANDARD C++ INCLUDES
#include <iostream>
#include <sstream>
using namespace std;

//OUR OWN CLASSES TO READ THE TREE
#include "../../UADataFormat/src/MassParticles.h"
#include "../../UADataFormat/src/MyBaseJet.h"
#include "../../UADataFormat/src/MyBeamSpot.h"
#include "../../UADataFormat/src/MyCaloJet.h"
#include "../../UADataFormat/src/MyCastorDigi.h"
#include "../../UADataFormat/src/MyCastorJet.h"
#include "../../UADataFormat/src/MyCastorRecHit.h"
#include "../../UADataFormat/src/MyDiJet.h"
#include "../../UADataFormat/src/MyElectron.h"
#include "../../UADataFormat/src/MyEvtId.h"
#include "../../UADataFormat/src/MyFwdGap.h"
#include "../../UADataFormat/src/MyGenJet.h"
#include "../../UADataFormat/src/MyGenKin.h"
#include "../../UADataFormat/src/MyGenMet.h"
#include "../../UADataFormat/src/MyGenPart.h"
#include "../../UADataFormat/src/MyHLTrig.h"
#include "../../UADataFormat/src/MyJet.h"
#include "../../UADataFormat/src/MyL1Trig.h"
#include "../../UADataFormat/src/MyL1TrigOld.h"
#include "../../UADataFormat/src/MyMITEvtSel.h"
#include "../../UADataFormat/src/MyMet.h"
#include "../../UADataFormat/src/MyMuon.h"
#include "../../UADataFormat/src/MyPFCand.h"
#include "../../UADataFormat/src/MyPFJet.h"
#include "../../UADataFormat/src/MyPUSumInfo.h"
#include "../../UADataFormat/src/MyPart.h"
#include "../../UADataFormat/src/MySimVertex.h"
#include "../../UADataFormat/src/MyTracks.h"
#include "../../UADataFormat/src/MyVertex.h"

bool isMC  = true;
bool debug = false;

void template_UATreeReader(const Int_t nevt_max = 100){
  

  //Declaration of my th1s
  TH1F* pt_gen  = new TH1F("pt_gen" , "pt_gen;pt;nTracks" , 120 , 0 , 6);
  TH1F* pt_reco = new TH1F("pt_reco", "pt_reco;pt;nTracks" , 120 , 0 , 6);
  
  //getting the list of files
  //vector<TString>* vfiles = getListOfFiles("/mydir/toto/UABaseTree*.root"); //Needs the fileManager.C from UAmulti
  vector<TString>* vfiles = new vector<TString>(1,"UABaseTree.root"); 
  //vfiles->push_back("UABaseTree2.root"); 
  
  //Declaration of tree and its branches variables
  TTree* tree = new TTree("evt","");
  MyEvtId*           evtId        = NULL ;
  vector<MyGenPart>* genPart      = NULL;
  vector<MyTracks>*  tracks       = NULL;
  vector<MyVertex>*  vertex       = NULL;
  
  // etc ....
  //Put All Classes you want to read here !!
  
  int i_tot = 0 , nevt_tot = 0;
  //starting Loop over files, stops at end of list of files or when reached nevt_max
  for(vector<TString>::iterator itfiles = vfiles->begin() ; itfiles != vfiles->end() && i_tot < nevt_max ; ++itfiles){
  
    TFile* file = TFile::Open(*itfiles,"READ");
    
    //getting the tree form the current file
    tree = (TTree*) file->Get("evt");
    
    //adding branches to the tree ----------------------------------------------------------------------
    tree->SetBranchAddress("evtId",&evtId);
    tree->SetBranchAddress("generalTracks",&tracks); 
    tree->SetBranchAddress("offlinePrimaryVertices",&vertex);
    
    if(isMC) tree->SetBranchAddress("genPart",&genPart);
  

    //Getting number of events
    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;
  
    //starting loop over events, stops when reached end of file or nevt_max
    for(int i = 0; i < nev && i_tot < nevt_max; ++i , ++i_tot){
    
      //printing the % of events done every 10k evts
      if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;
    
      //Filling the variables defined setting branches
      tree->GetEntry(i);


      //AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !


      //-------------------------------------------------------------------------------------------------
      //filling pt distribution for the generated particles
      //ie those from pythia generator, without reconstruction
      if(isMC){
        for(vector<MyGenPart>::iterator p=genPart->begin() ; p!=genPart->end() ; p++ )
          pt_gen->Fill(p->Pt());
      }
      
      

      //-------------------------------------------------------------------------------------------------
      //filling pt distribution for the observed/reconstructed tracks in the detector
      for(vector<MyTracks>::iterator it_tr = tracks->begin() ; it_tr != tracks->end() ; ++it_tr)
        pt_reco->Fill(it_tr->Pt());



    }//end of loop over events
    
    //Closing current files
    file->Close();
    
  }//end of loop over files
  
  
  //output file
  TFile* output = new TFile("output.root","RECREATE");
  output->cd();
  
  pt_gen->Write();
  pt_reco->Write();

  output->Close();
}
