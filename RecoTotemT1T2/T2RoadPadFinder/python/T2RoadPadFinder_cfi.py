import FWCore.ParameterSet.Config as cms

T2RoadPadFinder = cms.EDProducer("T2RoadPadFinder",
 
    VplaneToExclude= cms.vint32(),
    QuartedSelected = cms.vint32(0,1,2,3),
    MinCluSize_considered_asBlobs = cms.int32(5),
    #Tube size                             
    TolleranceThetaX = cms.double(0.002), 
    TolleranceThetaY = cms.double(0.002), 
    TolleranceDX = cms.double(1.0), 
    TolleranceDY = cms.double(1.0),                                  
    InefficiencyMaxJump= cms.int32(2),
    AllowsPadReAssociation = cms.bool(False),
    AllowsConcurrentBranches= cms.bool(False),
    Nmin_padsFinal = cms.int32(4),
    TwoPointsTubesAngularCollinearity= cms.double(0.05),                            MinimumNumCl1Hit= cms.int32(3),
    chi2XYProb_Thr= cms.double(0.01),  
    HitLabel = cms.string('T2Hits'),
    CluLabel = cms.string('T2MCl'),                        
    T2RoadCollProdName = cms.string('T2RoadColl'),                         
    verbosity  = cms.int32(0),
    ResolveOverlapDoubleCount= cms.bool(True),
    BiggestTubeAngleConsidered= cms.double(100.),
    NumSigma= cms.double(2.0),
    NumPadCluOccupancyAlert= cms.double(50.),
    useStraightPadTowers= cms.bool(False),

    OverlapDoubleCountDR = cms.double(2.0), #Depend on your alignment Resol  
    OverlapDoubleCountDPhi =cms.double(3.5),
    OverlapDoubleCountDTheta =  cms.double(0.001)  #This is DTheta of the averages                    
)


