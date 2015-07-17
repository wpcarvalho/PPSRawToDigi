#include "UATree/UABaseTree/interface/UABaseTree.h"

Bool_t FilterEventsDebug = true;

Bool_t UABaseTree::FilterEvents(){
  
  
  //L1 CUT
 /* if(storeL1Trig_){
    if(    L1Trig.GetTechDecisionBefore(36)
        || L1Trig.GetTechDecisionBefore(37)
        || L1Trig.GetTechDecisionBefore(38)
        || L1Trig.GetTechDecisionBefore(39)
      ) return 0;
  }*/
  
  //HLT CUT
  /*if(HLTrig.HLTmap.size() > 0){
    if(! (   HLTrig.HLTmap["HLT_L1_BscMinBiasOR_BptxPlusORMinus"] 
          || HLTrig.HLTmap["HLT_PixelTracks_Multiplicity40"]
          || HLTrig.HLTmap["HLT_PixelTracks_Multiplicity70"]
          || HLTrig.HLTmap["HLT_PixelTracks_Multiplicity85"]
	  || HLTrig.HLTmap["HLT_MinBiasPixel_SingleTrack"]   )
      )  return 0;
  }*/
  //2.76TeV Menu
  if(HLTrig.HLTmap.size() > 0){
    if(! (   HLTrig.HLTmap["HLT_Jet20_v1"]
          || HLTrig.HLTmap["HLT_Jet40_v1"] 
          || HLTrig.HLTmap["HLT_Jet60_v1"] 
          || HLTrig.HLTmap["HLT_L1BscMinBiasORBptxPlusANDMinus_v1"]
          || HLTrig.HLTmap["HLT_PixelTracks_Multiplicity50_Loose"]
          || HLTrig.HLTmap["HLT_PixelTracks_Multiplicity60_Loose"]
          || HLTrig.HLTmap["HLT_PixelTracks_Multiplicity70_Loose"]
          || HLTrig.HLTmap["HLT_ZeroBiasPixel_SingleTrack_v1"]
	  || HLTrig.HLTmap["HLT_ZeroBias_v1"]   )
      )  return 0;
  }
  
  //Vtx cut 
  Double_t vtxz_cut = 30.;
  Bool_t   goodVtx  = 0 ;
  for(map<string,vector<MyVertex> >::iterator vtxcoll = allVertices.begin() ; vtxcoll != allVertices.end() ; ++ vtxcoll){
    for(vector<MyVertex>::iterator vtx = vtxcoll->second.begin() ; vtx != vtxcoll->second.end() ; ++vtx){
      if(vtx->validity && fabs(vtx->z) < vtxz_cut && vtx->ntracks > 0 ){
        goodVtx = true;
	break;
      }
    }
    if(goodVtx) break;
  }
  if(!goodVtx) return 0;
  
  
  //If passed everything, keep event !
  return 1;
}
