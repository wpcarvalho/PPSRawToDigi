#include "TotemRPValidation/RPGeant4Validation/interface/RPOscarSimReader.h"

#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include "ctype.h"

RPOscarSimReader::RPOscarSimReader(const std::string &file_name)
{
  verbosity_ = 0;
  debug_info_ = true;

  TChain *chain = new TChain("ntuple");
  chain->Add(file_name.c_str());
  Init(chain);
  entries_no = fChain->GetEntries();
  entries_counter = 0;
}


bool RPOscarSimReader::RetrieveNewEvent()
{
  the_event_hits.clear();
  the_event_debug_info.clear();
  
  if(entries_counter >= entries_no)
    return false;

  if(fChain == 0)
    return false;

  Long64_t ientry = LoadTree(entries_counter);
  if (ientry < 0)
    return false;
  fChain->GetEntry(entries_counter);

  Long64_t curr_event = (Long64_t) Event;
  while(curr_event == Event)
  {
    double theta = TMath::ATan2(sqrt(p_x*p_x+p_y*p_y),p_z);
    double phi = TMath::ATan2(p_y, p_x);
    int unit_id = (int)UnitID;
    if(unit_id>=0)
    {
      the_event_hits.push_back( PSimHit( 
          Local3DPoint(lx, ly, lz), Local3DPoint(lx_ex, ly_ex, lz_ex), 
          PABS, 0, ELoss, (int)Ptype, (int)UnitID, (int)/*TrackID*/ParentID, theta, phi, (unsigned short)prim_ver_id) );
      //track id contains ParentID, process_type = prim_ver_id - the id for the place of the creation of secondaries
    }
    else if(debug_info_)
    {
      if(unit_id==-1)  //get primary protons
      {
        //std::cout<<"get primary protons"<<std::endl;
        the_event_debug_info.push_back( RPPSimHitDebugInfo( 
            Local3DPoint(0, 0, 0), Local3DPoint(0, 0, 0), 
            PABS, 0, 0, (int)Ptype, (int)UnitID, (int)/*TrackID*/ParentID, theta, phi, 0,
            Local3DPoint(vx, vy, vz), Local3DPoint(vx, vy, vz), 
            -1, (int)ParentID, -1) );
      }
      else if(unit_id==-2 || unit_id==-6 || unit_id==-7 || unit_id==-8 || unit_id==-10)  //get particles leaving the station or entering
      {
        //std::cout<<"get particles leaving the station"<<std::endl;
        the_event_debug_info.push_back( RPPSimHitDebugInfo( 
            Local3DPoint(0, 0, 0), Local3DPoint(0, 0, 0), 
            PABS, 0, 0, (int)Ptype, (int)UnitID, (int)/*TrackID*/ParentID, theta, phi, 0,
            Local3DPoint(x, y, z), Local3DPoint(vx, vy, vz), 
            -1, (int)ParentID, -1) );
      }
      else if(unit_id==-3 || unit_id==-4)
      {
        the_event_debug_info.push_back( RPPSimHitDebugInfo( 
            Local3DPoint(0, 0, 0), Local3DPoint(0, 0, 0), 
            PABS, 0, 0, (int)Ptype, (int)UnitID, (int)/*TrackID*/ParentID, theta, phi, 0,
            Local3DPoint(x, y, z), Local3DPoint(vx, vy, vz), 
            -1, (int)ParentID, (int)ELoss*10) );
      }
      else if(unit_id==-5)  //primary proton track stopped inside RP station
      {
        //std::cout<<"gprimary proton track stopped inside RP station"<<std::endl;
        the_event_debug_info.push_back( RPPSimHitDebugInfo( 
            Local3DPoint(0, 0, 0), Local3DPoint(0, 0, 0), 
            PABS, 0, 0, (int)Ptype, (int)UnitID, (int)/*TrackID*/ParentID, theta, phi, 0,
            Local3DPoint(x, y, z), Local3DPoint(vx, vy, vz), 
            (int)prim_ver_id, (int)ParentID, (int)ELoss*10, (int)ParentID));
      }
      else if(unit_id==-9)  //particle leaves the front wall of the rp
      {
        //std::cout<<"particle leaves the front wall of the rp"<<std::endl;
        the_event_debug_info.push_back( RPPSimHitDebugInfo( 
            Local3DPoint(0, 0, 0), Local3DPoint(0, 0, 0), 
            PABS, 0, 0, (int)Ptype, (int)UnitID, (int)/*TrackID*/ParentID, theta, phi, 0,
            Local3DPoint(x, y, z), Local3DPoint(vx, vy, vz), 
            -1, (int)ParentID, (int)ELoss*10) );
      }
    }
    
    ++entries_counter;
    
    if(entries_counter < entries_no)
    {
      Long64_t ientry = LoadTree(entries_counter);
      if (ientry < 0) break;
      fChain->GetEntry(entries_counter);
    }
    else
      break;
    if(verbosity_)
    {
      std::cout<<"Assembling hits, event: "<<Event<<std::endl;
    }
  }
  //std::cout<<"debughits vector size="<<the_event_debug_info.size()<<std::endl;
  
  if(verbosity_)
  {
    std::cout<<"PSimHitVector size: "<<the_event_hits.size()<<std::endl;
  }
  return true;
}


//TChain *RPOscarSimReader::ChainedFiles(const char *list_file_name, int max_files_no)
//{
//  fstream in;
//  in.open(list_file_name, ios::in);
//  
//  if(verbosity_)
//  {
//    std::cout<<"file with files: "<<list_file_name<<std::endl;
//  }
//  
//  if(in.bad())
//    return NULL;
//    
//  TChain *chain = new TChain("ntuple");
//  int file_no=0;
//  std::string name;
//  
//  while(!in.eof() && in.good() && file_no<max_files_no)
//  {
//    name = "";
//    in>>name;
//    while( !in.eof() && isspace(in.peek()) )
//      in.get();
//      
//    if(name.size())
//    {
//      chain->Add(name.c_str());
//      file_no++;
//      if(verbosity_)
//      {
//        std::cout<<file_no<<" "<<name<<std::endl;
//      }
//    }
//    else
//      break;
//  }
//  in.close();
//  return chain;
//}


RPOscarSimReader::~RPOscarSimReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RPOscarSimReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RPOscarSimReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void RPOscarSimReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if (tree == 0) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event",&Event);
   fChain->SetBranchAddress("UnitID",&UnitID);
   fChain->SetBranchAddress("Ptype",&Ptype);
   fChain->SetBranchAddress("TrackID",&TrackID);
   fChain->SetBranchAddress("ParentID",&ParentID);
   fChain->SetBranchAddress("ELoss",&ELoss);
   fChain->SetBranchAddress("PABS",&PABS);
   fChain->SetBranchAddress("vx",&vx);
   fChain->SetBranchAddress("vy",&vy);
   fChain->SetBranchAddress("vz",&vz);
   fChain->SetBranchAddress("x",&x);
   fChain->SetBranchAddress("y",&y);
   fChain->SetBranchAddress("z",&z);
   fChain->SetBranchAddress("lx",&lx);
   fChain->SetBranchAddress("ly",&ly);
   fChain->SetBranchAddress("lz",&lz);
   fChain->SetBranchAddress("x_ex",&x_ex);
   fChain->SetBranchAddress("y_ex",&y_ex);
   fChain->SetBranchAddress("z_ex",&z_ex);
   fChain->SetBranchAddress("lx_ex",&lx_ex);
   fChain->SetBranchAddress("ly_ex",&ly_ex);
   fChain->SetBranchAddress("lz_ex",&lz_ex);
   fChain->SetBranchAddress("p_x",&p_x);
   fChain->SetBranchAddress("p_y",&p_y);
   fChain->SetBranchAddress("p_z",&p_z);
   fChain->SetBranchAddress("prim_ver_id",&prim_ver_id);
   Notify();
}

Bool_t RPOscarSimReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.

   // Get branch pointers
   b_Event = fChain->GetBranch("Event");
   b_UnitID = fChain->GetBranch("UnitID");
   b_Ptype = fChain->GetBranch("Ptype");
   b_TrackID = fChain->GetBranch("TrackID");
   b_ParentID = fChain->GetBranch("ParentID");
   b_ELoss = fChain->GetBranch("ELoss");
   b_PABS = fChain->GetBranch("PABS");
   b_vx = fChain->GetBranch("vx");
   b_vy = fChain->GetBranch("vy");
   b_vz = fChain->GetBranch("vz");
   b_x = fChain->GetBranch("x");
   b_y = fChain->GetBranch("y");
   b_z = fChain->GetBranch("z");
   b_lx = fChain->GetBranch("lx");
   b_ly = fChain->GetBranch("ly");
   b_lz = fChain->GetBranch("lz");
   b_x_ex = fChain->GetBranch("x_ex");
   b_y_ex = fChain->GetBranch("y_ex");
   b_z_ex = fChain->GetBranch("z_ex");
   b_lx_ex = fChain->GetBranch("lx_ex");
   b_ly_ex = fChain->GetBranch("ly_ex");
   b_lz_ex = fChain->GetBranch("lz_ex");
   b_p_x = fChain->GetBranch("p_x");
   b_p_y = fChain->GetBranch("p_y");
   b_p_z = fChain->GetBranch("p_z");
   b_prim_ver_id = fChain->GetBranch("prim_ver_id");

   return kTRUE;
}

void RPOscarSimReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RPOscarSimReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


void RPOscarSimReader::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L RPOscarSimReader.C
//      Root > RPOscarSimReader t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
