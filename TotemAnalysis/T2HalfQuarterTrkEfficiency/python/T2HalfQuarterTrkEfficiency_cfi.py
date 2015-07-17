import FWCore.ParameterSet.Config as cms

T2HalfQuarterTrkEfficiency = cms.EDAnalyzer("T2HalfQuarterTrkEfficiency",

	 
         verbosity=cms.bool(False),
         UnassHitThr=cms.double(0.2),
         RefTrkHitMult=cms.uint32(4),                                    
         RefTrkMultiplicity=cms.uint32(3),                            
         maxdphihit=cms.double(3.0), 
	 maxdrhit=cms.double(1.0), 
	 MaxDphi=cms.double(12.0),
         MaxPad=cms.uint32(10),
         MaxStrip=cms.uint32(10),                   

         xmlfilenameFull = cms.string(''),
         xmlfilenameNotDead = cms.string(''),                   

         Chi2ProbRefTrk=cms.double(0.01),
	 UseRestrictedOverlap=cms.bool(False),
         OldReco=cms.bool(False),
                                     
	 Trk_RSeparation=cms.double(1.0),
         Trk_PhiSeparation=cms.double(3.0),
         Trk_EtaSeparation=cms.double(0.15),                   
         MaxPadCluInOverlapUporDown =cms.uint32(23),                           
         LookToRawEvent=cms.bool(False),
         UseUncorrupetdEventMap=cms.bool(False),
         OnlyShadowAnalysis=cms.bool(False),
        
         TrkEtamin= cms.double (5.1),  # 4.8 for 10 GeV
         TrkEtaMAX= cms.double (6.9),  # 7.0 for 10 GeV
         RefCluLabel = cms.string('T2MCl'),    
         CluLabel = cms.string('T2MCl'),
         HitLabel = cms.string('T2Hits'),
         RoadLabel = cms.string('T2RoadPadFinderA'),#T2RoadColl2 apparently don't used
         TrackLabel = cms.string('T2TrackColl3A'),#T2TrackColl2 apparently don't used
	 HitLabelTestedQ = cms.string('T2Hits'),# apparently don't used


         RoadLabelH0 = cms.string('T2RoadPadFinderA'),
         TrackLabelH0 = cms.string('T2TrackColl3A'),
         RoadLabelH1 = cms.string('T2RoadPadFinderB'),
         TrackLabelH1 = cms.string('T2TrackColl3B'),

         RoadLabelH2 = cms.string('T2RoadPadFinderA'),#T2RoadColl2H2
         TrackLabelH2 = cms.string('T2TrackColl3A'),#T2TrackColl2H2
         RoadLabelH3 = cms.string('T2RoadPadFinderB'),#T2RoadColl2H3
         TrackLabelH3 = cms.string('T2TrackColl3B'),#T2TrackColl2H3
                                            
         TrackLabelFirstPlanes= cms.string('T2TrackColl3A'),
         TrackLabelLastPlanes= cms.string('T2TrackColl3B'),
         PadCluLabelFirstPlanes= cms.string('T2MClSTDA'),
         PadCluLabelLastPlanes= cms.string('T2MClSTDB'),                              
	 Refstring = cms.string('A'), # choose First->A or Last group  of planes in the half as a Ref
         AnType = cms.string('INCLUSIVE'),#INCLUSIVE DOUBLE SINGLE
                                            
                                            
         OutputFile = cms.untracked.string('valPythia90T2PlotsReco.root'),        
         AllowedDRTrackDistance=cms.double(1.0),
         ShadowTan=cms.double(0.7),
         MaxTrkInProcess=cms.uint32(10),                            
         ReferenceQuarters=cms.vuint32(0),
	 T2PadDigiCollectionLabel = cms.InputTag("T2Digis", "T2PadDigi"),
	 T2StripDigiCollectionLabel = cms.InputTag("T2Digis", "T2StripDigi"),
	 T2VfatInformationLabel = cms.InputTag("Raw2DigiProducer")
    )
