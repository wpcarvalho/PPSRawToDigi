//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 28 21:50:22 2006 by ROOT version 5.08/00
// from TTree ntuple/Ntuple
// found on file: Memory Directory
//////////////////////////////////////////////////////////

#ifndef TotemRPValidation_RPGeant4Validation_RPOscarSimReader_h
#define TotemRPValidation_RPGeant4Validation_RPOscarSimReader_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

#include "TMath.h"
#include "fstream"
#include "TChain.h"
#include <vector>
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPPSimHitDebugInfo.h"


class RPOscarSimReader
{
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Double_t        Event;
   Double_t        UnitID;
   Double_t        Ptype;
   Double_t        TrackID;
   Double_t        ParentID;
   Double_t        ELoss;
   Double_t        PABS;
   Double_t        vx;
   Double_t        vy;
   Double_t        vz;
   Double_t        x;
   Double_t        y;
   Double_t        z;
   Double_t        lx;
   Double_t        ly;
   Double_t        lz;
   Double_t        x_ex;
   Double_t        y_ex;
   Double_t        z_ex;
   Double_t        lx_ex;
   Double_t        ly_ex;
   Double_t        lz_ex;
   Double_t        p_x;
   Double_t        p_y;
   Double_t        p_z;
   Double_t        prim_ver_id;

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_UnitID;   //!
   TBranch        *b_Ptype;   //!
   TBranch        *b_TrackID;   //!
   TBranch        *b_ParentID;   //!
   TBranch        *b_ELoss;   //!
   TBranch        *b_PABS;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_lx;   //!
   TBranch        *b_ly;   //!
   TBranch        *b_lz;   //!
   TBranch        *b_x_ex;   //!
   TBranch        *b_y_ex;   //!
   TBranch        *b_z_ex;   //!
   TBranch        *b_lx_ex;   //!
   TBranch        *b_ly_ex;   //!
   TBranch        *b_lz_ex;   //!
   TBranch        *b_p_x;   //!
   TBranch        *b_p_y;   //!
   TBranch        *b_p_z;   //!
   TBranch        *b_prim_ver_id;   //!

   RPOscarSimReader(const std::string &file_name);
   virtual ~RPOscarSimReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   inline void ToBegin() {entries_counter=0;}
   inline int GetRetrievedEventHitsNo() {return the_event_hits.size();}
   bool RetrieveNewEvent();
   const std::vector<PSimHit>& GetRetrievedEvent()  {return the_event_hits;}
   const std::vector<RPPSimHitDebugInfo>& GetRetrievedEventDebugInfo()  {return the_event_debug_info;}
   
private:
//   TChain *ChainedFiles(const char *list_file_name, int max_files_no=1000);
   Long64_t entries_counter;
   Long64_t entries_no;
   std::vector<PSimHit> the_event_hits;
   std::vector<RPPSimHitDebugInfo> the_event_debug_info;
   int verbosity_;
   bool debug_info_;
};

#endif
