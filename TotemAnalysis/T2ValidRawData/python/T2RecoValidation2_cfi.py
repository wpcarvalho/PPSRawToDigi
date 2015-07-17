import FWCore.ParameterSet.Config as cms

T2ValidRaw = cms.EDAnalyzer("T2AnalyzerRaw",
        FitdzEnergy = cms.double (10.0),   # choose 10. or 50. GeV
        Chicut= cms.double (3.0),
        DZScale= cms.double (1.0),
        TrkEtamin= cms.double (5.1),  # 4.8 for 10 GeV
        TrkEtaMAX= cms.double (6.9),  # 7.0 for 10 GeV
        DispVtx=cms.bool(False),
	# DetEffRWind =cms.double(10.0),
	# DetEffPhiWind =cms.double(10.0),
	 Numhittrkeff =cms.int32(4),	
	 chiRredCut=cms.double(5.0),
         chiPhiredCut=cms.double(5.0),
	 AllowedDRTrackDistance=cms.double(1.0),
         Testedcamera = cms.vuint32(0,1,2,3,4),
         PhiChiProbCut=cms.double(0.01),
         RChiProbCut=cms.double(0.01),

	 EffMaxPad=cms.uint32(12),
         EffMaxStrip= cms.uint32(12),
 	 Effmaxdphihit=cms.double(3.0),
         Effmaxdrhit=cms.double(3.0),
         NoiseDphiMAX=cms.double(30.0),
         NoiseDrMAX=cms.double(20.0),
	 Effgoodhitnumber=cms.uint32(3),
	
	 useRZforResol=cms.uint32(0),
	 simufile=cms.bool(False),
         verbosity=cms.bool(False),
                            
         MaxPadAllowedInQuarter=cms.uint32(40),                 
	 MaxPad=cms.uint32(10),
         MaxStrip=cms.uint32(10),
	 MaxEvents=cms.uint32(2400),
	 maxdphihit=cms.double(3.0), 
	 maxdrhit=cms.double(1.0), 
	 NumHitGood=cms.uint32(3), 
	 MaxDphi=cms.double(12.0),
         OnlycorruptionAnalysis=cms.bool(False),                    
         OnlyClusterAnalysis=cms.bool(False),
         VFATMonitoring=cms.bool(False),
                            
	 FitgravcenterZ=cms.double(0.),
	 UseJointProb=cms.uint32(1),
	 #DoALign=cms.bool(True),
	 HitNumb4Align=cms.uint32(8),
	 MeasuredXYResol=cms.double(0.15),
	 SHIFTprescale=cms.double(1.0),
         Idreferencedet=cms.uint32(13),
         MaxStepalignstep=cms.uint32(100),
         AlignmentHitRMax=cms.double(110.0),

	 DetForNoiseStudies=cms.uint32(0), 
	 PhiMinForNoiseStudies=cms.double(270), 
	 PhiMaxForNoiseStudies=cms.double(359),

         ExcludeNoisyplane=cms.bool(False),
         CommonNoiseClSize=cms.uint32(12),

         SelectedHalf=cms.uint32(0),       
         MaxTrkInQuarter=cms.uint32(1),                        
         MinTrkInQuarter=cms.uint32(1),
         produceVfatEffiFile=cms.bool(False),
              
         skipSelectedEvents=cms.bool(False),
         skipEventFileName= cms.string('corruptedEvents.txt'),
                            
         CluLabel = cms.string('T2MCl'),
         HitLabel = cms.string('T2Hits'),
         RoadLabel = cms.string('T2RoadColl'),
         TrackLabel = cms.string('T2TrackColl2'),
         OutputFile = cms.untracked.string('valPythia90T2PlotsReco.root'),
         xmlfilenameFull = cms.string(''),
         xmlfilenameNotDead = cms.string(''),
         DeadSectFileName = cms.string(''),
         LookToRawEvent=cms.bool(False),
         UseUncorrupetdEventMap=cms.bool(False),
         requiregoodChi=cms.bool(False),                   
         bunchesToAnalyse=cms.vuint32()
    )
