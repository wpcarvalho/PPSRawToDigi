#ifndef DataFormats_T2Hit_Classes_H
#define DataFormats_T2Hit_Classes_H

#include <boost/cstdint.hpp>
#include <vector>
#include <string>

#include <DataFormats/Common/interface/Wrapper.h>
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/hit_entry.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2Hit/interface/T2HitCollectionWrapper.h"
#include "DataFormats/T2Hit/interface/T2HitCollectionMapping.h"
#include "DataFormats/T2Hit/interface/T2HitCollectionMappingWrapper.h"
#include "DataFormats/T2Hit/interface/T2PadStripAssociator.h"
#include "DataFormats/T2Hit/interface/T2PadStripAssociatorWrapper.h"
#include "DataFormats/T2Hit/interface/T2Hit_to_Track_Map.h"
#include "DataFormats/T2Hit/interface/T2Hit_to_Track_MapWrapper.h"



  namespace {
    T2Hit t2hit_;
    hit_entry t2hent;
    std::vector<hit_entry> t2hentv;	
    T2HitCollection t2hitColl_;
    T2HitCollectionWrapper t2hitCollW_; 
    T2HitCollectionMapping t2hitCollMap_;	
    T2HitCollectionMappingWrapper t2hitCollMapW_;
    T2PadStripAssociator t2padstrAss_;	
    T2PadStripAssociatorWrapper t2padstrAssW_;
 
    T2Hit_to_Track_Map t2hittotrk_;
    T2Hit_to_Track_MapWrapper t2hittotrkW_;

    std::vector<T2Hit> dumt2hitColl_;
    edm::Wrapper<std::vector<T2Hit> > dumt2hitCollW_;
     
    std::map<std::pair<long int,long int>,unsigned long int > dummap2_;
    edm::Wrapper<std::map<std::pair<long int,long int>,unsigned long int > > dummap1W_; 

    std::map<long int, vector<long int> > dumMapCluIdVsStrip_;
    edm::Wrapper<std::map<long int, vector<long int> > > dumMapCluIdVsStripWr_;

    std::map<long int, int > dumMapHittoTrk_; 
    edm::Wrapper<std::map<long int, int > > dumMapHittoTrkWr_;
  }

#endif // DataFormats_T2Hit_Classes_H
