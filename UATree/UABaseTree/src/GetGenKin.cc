// Description: Function to retrieve Generated Kinematic et al.

// DataFormats
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"


// UABaseTree Analysis class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"


bool GenKinDebug = false;

void UABaseTree::GetGenKin(const edm::Event& iEvent)
{
   genKin.Reset();

   // MC Process Id and PtHat for RECO
   /*Handle<HepMCProduct> hepMCHandle;
   iEvent.getByLabel(hepmc_, hepMCHandle ) ;
   
   if( hepMCHandle.isValid() ){
      const HepMC::GenEvent* GenEvt = hepMCHandle->GetEvent() ;
      genKin.MCProcId = GenEvt->signal_process_id();
      genKin.PtHat = GenEvt->event_scale();

      //-- PDF Info

      const HepMC::PdfInfo* pdfInfo = GenEvt->pdf_info();

      genKin.x1      = pdfInfo->x1();
      genKin.x2      = pdfInfo->x2();
      genKin.Q       = pdfInfo->scalePDF();
      genKin.Part1Id = pdfInfo->id1();
      genKin.Part2Id = pdfInfo->id2();
   }*/
   Handle<GenEventInfoProduct> genEvtInfoH;
   iEvent.getByLabel(hepmc_, genEvtInfoH ) ;
   
   if( genEvtInfoH.isValid() ){
      genKin.MCProcId = genEvtInfoH->signalProcessID();
      const std::vector<double>& binningValues = genEvtInfoH->binningValues();
      if( binningValues.size() ) genKin.PtHat = binningValues.at(0);
      genKin.genWeight = genEvtInfoH->weight();

      //-- PDF Info
      gen::PdfInfo const* pdfInfo = genEvtInfoH->pdf();
      if(pdfInfo){
	 genKin.x1      = pdfInfo->x.first;
	 genKin.x2      = pdfInfo->x.second;
	 genKin.Q       = pdfInfo->scalePDF;
	 genKin.Part1Id = pdfInfo->id.first;
	 genKin.Part2Id = pdfInfo->id.second;
      }
   }

   //K Factor For Signal
   edm::Handle<double> KFactor; 
   iEvent.getByLabel("KFactorProducer",KFactor);
   if( KFactor.isValid() ){
    genKin.kfactor = *KFactor;
   }
   
   if(GenKinDebug) genKin.Print();
}

