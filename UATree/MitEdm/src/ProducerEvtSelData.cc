// $Id: ProducerEvtSelData.cc,v 1.1 2012/04/25 15:18:02 rougny Exp $

#include "UATree/MitEdm/interface/ProducerEvtSelData.h"
#include "UATree/MitEdm/interface/EvtSelData.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include <TMath.h>

using namespace std;
using namespace edm;
using namespace mitedm;

//--------------------------------------------------------------------------------------------------
ProducerEvtSelData::ProducerEvtSelData(const ParameterSet& cfg)
  : srcHF_(cfg.getUntrackedParameter<string>("hfRecHits","hfreco")),
    srcHBHE_(cfg.getUntrackedParameter<string>("hbheRecHits","hbhereco")),
    srcCastor_(cfg.getUntrackedParameter<string>("castorRecHits","castorreco")),
    srcZDC_(cfg.getUntrackedParameter<string>("zdcRecHits","zdcreco")),
    srcPixels_(cfg.getUntrackedParameter<string>("pixelRecHits","siPixelRecHits")),
    srcVertex_(cfg.getUntrackedParameter<string>("vertex","")),
    srcTowers_(cfg.getUntrackedParameter<string>("caloTowers","towerMaker")),
    hfEthresh_(cfg.getUntrackedParameter<double>("hfEThreshold",3.)),
    hfETowerh_(cfg.getUntrackedParameter<double>("hfETowerThreshold",3.)),
    srcTrk_(cfg.getUntrackedParameter<string>("tracks","generalTracks"))
{
  // Constructor.

  produces<EvtSelData>();
}

//--------------------------------------------------------------------------------------------------
ProducerEvtSelData::~ProducerEvtSelData()
{
  // Destructor.
}

//--------------------------------------------------------------------------------------------------
void ProducerEvtSelData::produce(Event &evt, const EventSetup &setup)
{
  // Produce event selection data for this event.

  double eHfNeg       = 0;
  double eHfPos       = 0;
  double eHfPosTime   = 0;
  double eHfNegTime   = 0;
  int    eHfPcounts   = 0;
  int    eHfNcounts   = 0;
  double eHcalNeg     = 0;
  double eHcalPos     = 0;
  double eCaNeg       = 0;
  double eCaPos       = 0;
  double eCaPosTime   = 0;
  double eCaNegTime   = 0;
  double eZdNeg       = 0;
  double eZdPos       = 0;
  double eZdPosTime   = 0;
  double eZdNegTime   = 0;
  int    ePxbHits     = 0;
  int    ePxHits      = 0;
  double eClusVtxQual = 0;
  double eClusVtxDiff = 0;
  int    nHfTowersP      = 0;
  int    nHfTowersN      = 0;
  double sumEsubEpPos    = 0;
  double sumEaddEpPos    = 0;
  double sumHfEsubEpPos  = 0;
  double sumHfEaddEpPos  = 0;
  double sumEsubEpNeg    = 0;
  double sumEaddEpNeg    = 0;
  double sumHfEsubEpNeg  = 0;
  double sumHfEaddEpNeg  = 0;
  double eHPTrkFrac      = 0;

  // HF
  Handle<HFRecHitCollection> hfhits;
  try {
    evt.getByLabel(edm::InputTag(srcHF_),hfhits);
  } catch (...) {}  

  if (hfhits.isValid()) {
    for(size_t ihit = 0; ihit<hfhits->size(); ++ihit){
      const HFRecHit h = (*hfhits)[ihit];
      double energy = h.energy();
      double time = h.time();
      const HcalDetId id(h.id());
      if(energy<0) continue; 
      if (id.zside()<0) {
        eHfNeg     += energy;
        eHfNegTime += energy * time;
        if (energy>hfEthresh_)
          ++eHfNcounts;
      } else {
        eHfPos     += energy;
        eHfPosTime += energy * time;
        if (energy>hfEthresh_)
          ++eHfPcounts;
      }
    }
  }

  // HCAL (barrel + endcap)
  Handle<HBHERecHitCollection> hbhehits;
  try {
    evt.getByLabel(edm::InputTag(srcHBHE_),hbhehits);
  } catch (...) {}  

  if (hbhehits.isValid()) {
    for(size_t ihit = 0; ihit<hbhehits->size(); ++ihit){
      const HBHERecHit h = (*hbhehits)[ihit];
      double energy = h.energy();
      const HcalDetId id(h.id()); 
      if(energy<0) continue;
      if (id.zside()<0) {
        eHcalNeg   += energy;
      } else {
        eHcalPos   += energy;
      }
    }
  }

  // CASTOR
  Handle<CastorRecHitCollection> castorhits;
  try {
    evt.getByLabel(edm::InputTag(srcCastor_),castorhits);
  } catch (...) {}  

  if (castorhits.isValid()) {
    for(size_t ihit = 0; ihit<castorhits->size(); ++ihit){
      const CastorRecHit h = (*castorhits)[ihit];
      double energy = h.energy();
      double time = h.time();
      const HcalCastorDetId id(h.id()); 
      if(energy<0) continue;
      if (id.zside()<0) {
        eCaNeg     += energy;
        eCaNegTime += energy * time;
      } else {
        eCaPos     += energy;
        eCaPosTime += energy * time;
      }
    }
  }

  // ZDC
  Handle<ZDCRecHitCollection> zdchits;
  try {
    evt.getByLabel(edm::InputTag(srcZDC_),zdchits);
  } catch (...) {}  

  if (zdchits.isValid()) {
    for(size_t ihit = 0; ihit<zdchits->size(); ++ihit){
      const ZDCRecHit h = (*zdchits)[ihit];
      double energy = h.energy();
      double time = h.time();
      const HcalZDCDetId id(h.id()); 
      if(energy<0) continue;
      if (id.zside()<0) {
        eZdNeg     += energy;
        eZdNegTime += energy * time;
      } else {
        eZdPos     += energy;
        eZdPosTime += energy * time;
      }
    }
  }

  // time values are energy-weighted averages
  if(eHfPos>0)
    eHfPosTime /= eHfPos;
  if(eHfNeg>0)
    eHfNegTime /= eHfNeg;
  if(eCaPos>0)
    eCaPosTime /= eCaPos;
  if(eCaNeg>0)
    eCaNegTime /= eCaNeg;
  if(eZdPos>0)
    eZdPosTime /= eZdPos;
  if(eZdNeg>0)
    eZdNegTime /= eZdNeg;
  

  // pixel cluster shape compatibility with vertex
  Handle<SiPixelRecHitCollection> hRecHits;
  try {
    evt.getByLabel(edm::InputTag(srcPixels_),hRecHits);
  } catch (...) {}

  if (hRecHits.isValid()) {
    ESHandle<TrackerGeometry> trackerHandle;
    setup.get<TrackerDigiGeometryRecord>().get(trackerHandle);
    const TrackerGeometry *tgeo = trackerHandle.product();
    const SiPixelRecHitCollection *hits = hRecHits.product();

    vector<VertexHit> vhits;
    for(SiPixelRecHitCollection::DataContainer::const_iterator hit = hits->data().begin(), 
          end = hits->data().end(); hit != end; ++hit) {
      if (!hit->isValid())
        continue;
      ++ePxHits;
      DetId id(hit->geographicalId());
      if(id.subdetId() != int(PixelSubdetector::PixelBarrel))
        continue;
      ++ePxbHits;
      const PixelGeomDetUnit *pgdu = static_cast<const PixelGeomDetUnit*>(tgeo->idToDet(id));
      if (1) {
        const RectangularPixelTopology *pixTopo = 
          static_cast<const RectangularPixelTopology*>(&(pgdu->specificTopology()));
        vector<SiPixelCluster::Pixel> pixels(hit->cluster()->pixels());
        bool pixelOnEdge = false;
        for(std::vector<SiPixelCluster::Pixel>::const_iterator pixel = pixels.begin(); 
            pixel != pixels.end(); ++pixel) {
          int pixelX = pixel->x;
          int pixelY = pixel->y;
          if(pixTopo->isItEdgePixelInX(pixelX) || pixTopo->isItEdgePixelInY(pixelY)) {
            pixelOnEdge = true;
            break;
          }
        }
        if (pixelOnEdge)
          continue;
      }
      
      LocalPoint lpos = LocalPoint(hit->localPosition().x(),
                                   hit->localPosition().y(),
                                   hit->localPosition().z());
      GlobalPoint gpos = pgdu->toGlobal(lpos);
      VertexHit vh;
      vh.z = gpos.z(); 
      vh.r = gpos.perp(); 
      vh.w = hit->cluster()->sizeY();
      vhits.push_back(vh);
    }

    double zest = 0.0;
    if (srcVertex_.length()!=0) { // use specified vertex
      edm::Handle<reco::VertexCollection> vertexCol;
      evt.getByLabel(srcVertex_,vertexCol);
      const reco::VertexCollection *vertices = vertexCol.product();
      unsigned int vtracks = 0;
      for(reco::VertexCollection::const_iterator
            vertex = vertices->begin(); vertex!= vertices->end(); ++vertex) {
        if(vertex->tracksSize()>vtracks) {
          vtracks = vertex->tracksSize();
          zest = vertex->z();
        }
      }
    } else { // or use pixel cluster-shape vertex
      int nhits = 0, nhits_max = 0;
      double chi = 0, chi_max = 1e+9;
      for(double z0 = -15.9; z0 <= 15.95; z0 += 0.1) {
        nhits = getContainedHits(vhits, z0, chi);
        if(nhits > 0) {
          if(nhits >  nhits_max) { 
            chi_max = 1e+9; 
            nhits_max = nhits; 
          }
          if(nhits >= nhits_max) {
            if(chi < chi_max) { 
              chi_max = chi; 
              zest = z0; 
            }
          }
        }
      }
    }

    double chi = 0;
    int nbest=0, nminus=0, nplus=0;
    nbest = getContainedHits(vhits,zest,chi);
    nminus = getContainedHits(vhits,zest-10.,chi);
    nplus = getContainedHits(vhits,zest+10.,chi);

    eClusVtxDiff = nbest - (nminus+nplus)/2.;
    if ((nminus+nplus)> 0)
      eClusVtxQual = (2.0*nbest)/(nminus+nplus);  // A/B
    else if (nbest>0)
      eClusVtxQual = 1000.0;                      // A/0 (set to arbitrarily large number)
    else
      eClusVtxQual = 0;                           // 0/0 (already the default)

  }

  // calo based variables
  edm::Handle<CaloTowerCollection> towers; 
  try {
    evt.getByLabel(edm::InputTag(srcTowers_),towers);
  } catch (...) {}  

  if (towers.isValid()) {
    for(CaloTowerCollection::const_iterator cal = towers->begin(); cal != towers->end(); ++cal) {
      if (cal->energy()<hfETowerh_)
        continue;

      double tmp = TMath::Cos(cal->theta());
      double Esub = cal->energy()*(1-tmp);
      double Eadd = cal->energy()*(1+tmp);

      if (cal->eta()>0) {
        sumEsubEpPos += Esub;
        sumEaddEpPos += Eadd;
      } else if (cal->eta()<0) {
        sumEsubEpNeg += Esub;
        sumEaddEpNeg += Eadd;
      }
      if (cal->eta()>3) {
        sumHfEsubEpPos += Esub;
        sumHfEaddEpPos += Eadd;
      } else if (cal->eta()<-3) {
        sumHfEsubEpNeg += Esub;
        sumHfEaddEpNeg += Eadd;
      }

      for(unsigned int i = 0; i < cal->constituentsSize(); ++i) {
        const DetId id = cal->constituent(i);
        if(id.det() != DetId::Hcal)
          continue;
        HcalSubdetector subdet=(HcalSubdetector(id.subdetId()));
        if(subdet != HcalForward)
          continue;
        if (cal->eta()<-3)
          ++nHfTowersN;
        if (cal->eta()>+3)
          ++nHfTowersP;
      }    
    }
  }

  // high-purity track fraction
  Handle<reco::TrackCollection> tkRef;
  try {
    evt.getByLabel(edm::InputTag(srcTrk_),tkRef);
  } catch (...) {}  

  if (tkRef.isValid()) {

    const reco::TrackCollection* tkColl = tkRef.product();
    int numhighpurity=0;
    reco::TrackBase::TrackQuality trackQuality(reco::TrackBase::qualityByName("highPurity"));
    if(tkColl->size()>0){ 
      reco::TrackCollection::const_iterator itk = tkColl->begin();
      reco::TrackCollection::const_iterator itk_e = tkColl->end();
      for(;itk!=itk_e;++itk){
        if(itk->quality(trackQuality)) numhighpurity++;
      }
      eHPTrkFrac = (double)numhighpurity/(double)tkColl->size();
    }
  }

  // fill EvtSelData object
  std::auto_ptr<EvtSelData> output(new EvtSelData(eHcalNeg, eHcalPos,
                                                  eHfNeg,eHfPos,eHfNegTime,eHfPosTime,
                                                  eCaNeg,eCaPos,eCaNegTime,eCaPosTime,
                                                  eZdNeg,eZdPos,eZdNegTime,eZdPosTime,
						  ePxbHits,ePxHits,
                                                  eClusVtxQual,eClusVtxDiff,
                                                  eHfPcounts, eHfNcounts,
                                                  nHfTowersP, nHfTowersN,
                                                  sumEsubEpPos,  sumEaddEpPos,
                                                  sumHfEsubEpPos,sumHfEaddEpPos,
                                                  sumEsubEpNeg,  sumEaddEpNeg,
                                                  sumHfEsubEpNeg,sumHfEaddEpNeg,
                                                  eHPTrkFrac));
  evt.put(output);
}

//--------------------------------------------------------------------------------------------------
int ProducerEvtSelData::getContainedHits(const std::vector<VertexHit> &hits, double z0, double &chi) 
{
  // Calculate number of hits contained in v-shaped window in cluster y-width vs. z-position.

  int n = 0;
  chi   = 0.;

  for(vector<VertexHit>::const_iterator hit = hits.begin(); hit!= hits.end(); hit++) {
    double p = 2 * fabs(hit->z - z0)/hit->r + 0.5; // FIXME
    if(TMath::Abs(p - hit->w) <= 1.) { 
      chi += fabs(p - hit->w);
      n++;
    }
  }
  return n;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ProducerEvtSelData);
