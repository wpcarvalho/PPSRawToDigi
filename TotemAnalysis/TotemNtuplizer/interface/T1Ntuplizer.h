
#include "FWCore/Utilities/interface/InputTag.h"
#include "TotemAnalysis/TotemNtuplizer/interface/Ntuplizer.h"
#include "TotemAnalysis/TotemNtuplizer/interface/T1Event.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfatCollection.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include "DataFormats/T1Cluster/interface/T1Cluster.h"
#include "DataFormats/T1Cluster/interface/T1ClusterCollection.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"



class T1Ntuplizer : public Ntuplizer
{
public:
  T1Ntuplizer(const edm::ParameterSet&);
  virtual ~T1Ntuplizer() {}
  
  virtual void CreateBranches(const edm::EventSetup&, TTree *);
  virtual void FillEvent(const edm::Event&, const edm::EventSetup&);

  T1Event t1obj;

private:
 
  std::auto_ptr<T1Geometry> layer;
  std::string RawDataName;
  std::string outputFileName, HepMCProductLabel, CluLabel, HitLabel, RoadLabel,trackLabel;
  edm::InputTag t1DigiWireCollectionLabel;
  edm::InputTag t1DigiVfatCollectionLabel;
  edm::InputTag t1RecHit2DCollectionLabel;

};


