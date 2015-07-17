import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('RPDataReduction',
    nPl=cms.uint32(3), #minimum number of planes to call a mulTrack
    nMaxPrPl=cms.uint32(10), #max hits per plane for multiTrack
    verb=cms.uint32(4),   #verbosity
    sigmaCut=cms.double(3.0), #nominal number of sigmas in y-edge-cut
    andNotOr=cms.uint32(0), #trigger mode: 0==u.or.v, 1==u.and.v, 3==u, 4==v
    elCut56=cms.double(8.4), #cut applied on track Y0() if reference Pot is 56_nr_* (120 or 121 == 56_nr_tp or 56_nr_bt)
                             #(mm from beam center)
   #arm 5-6, far pots are 1mm further from the beam than the near pots
    elCut45t=cms.double(8.5), #cut applied for reference pot 45_nr_tp (020)
    elCut45b=cms.double(0.0), #cut applied for reference pot 45_nr_bt (021)
    errX=cms.double(2.0), # total distance from edge demanded, in x (mm)
    errY=cms.double(0.3), # total distance from edge demanded, in y (mm)
    xLow=cms.untracked.double(-0.6), # x-cut on refPot track to raise real elastic purity of sample
    xHigh=cms.untracked.double(0.0), #same 
    whichDiag=cms.uint32(0), #4-pot calculated quantities: chooses which 4 pots to check
                      #0 == 020,024,121,125 (45tp*56bt - elastic)
                      #1 == 021,025,120,124 (45bt*56tp - elastic)
                      #2 == 020,024,120,124 (45tp*56tp)
                      #3 == 022,023,122,123 (45hr*56hr)
                      #4 == 021,025,121,125 (45bt*56bt)
                      #5 == 020,024,122,123 (45tp*56hr)
                      #6 == 021,025,122,123 (45bt*56hr)
                      #7 == 022,023,120,124 (45hr*56tp)
                      #8 == 022,023,121,125 (45hr*56bt)
    whichBunchAll=cms.int32(-1), #4-pot calculated quantities: all-bunches histogram is filled using the same
                      # "int to int" map as the per-bunch histograms, only using this fake bunch number
    bigBunch=cms.int32(0),
    dLxDs=cms.double(-0.536),                  #nominal betaStar=90m optics November 2011
    Ly=cms.double(263000.),                       #L_x=v_y=0 for this optics, Ly=263m, given here in mm
    Lx=cms.double(0.),                       #L_x=0...3m for this optics, given here in mm
    Vx=cms.double(-1.87),                       #v_y=0 for this optics, v is unitless
    dVxDs=cms.double(0.000056),                       #v_y=0 for this optics, dv/ds=0.056/m, given here in inverse mm
    Eb=cms.double(3500.),                       #Beam energy
    ElaSigmaX=cms.double(0.53),                  
    ElaSigmaAngleX=cms.double(0.000046),                  
    sigmaXi=cms.double(0.0087),
    sigmaRg=cms.double(0.94),
    readMultiTrk=cms.bool(True),
    readT2=cms.bool(False),
    readT1=cms.bool(False),
    readLoNeg=cms.bool(False),
    readClusters=cms.bool(False),
    readRecoProt=cms.bool(False),
    refPot=cms.uint32(24),
    oppositePot=cms.uint32(125),
    checkPot=cms.uint32(20),
    tmcModule = cms.string('RPCC'),
    tmcProd = cms.string('RPCCSimuBits'),
    t2trModule = cms.string('T2TrackColl3'),
    t2trProd = cms.string('T2TrackColl'),
    t1trModule = cms.string('t1tracks'),
    t1trProd = cms.string('T1TrackColl'),
    t2padModule = cms.string('T2MCl'),
    t2padProd = cms.string('T2PadClusters'),
    fileName = cms.string('RPDataReduction.root'),
    bunchText = cms.string('bunchCombinNums-rXXXX.txt'),

    RawEventLabel = cms.InputTag("source"),
    DetSetVectorRPDigClusterLabel = cms.InputTag("RPClustProd"),

    RPRecognizedPatternsCollectionLabel = cms.InputTag("NonParallelTrackFinder"),
    RPTrackCandidateCollectionLabel = cms.InputTag("RPSinglTrackCandFind"),
    RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit"),
    RPMulTrackCandidateCollectionLabel = cms.InputTag("RPMulTrackCandFind"),
    RPMulFittedTrackCollectionLabel = cms.InputTag("RPMulTrackCandCollFit"),

    RPReconstructedProtonCollectionLabel = cms.InputTag("RPCC")
)
