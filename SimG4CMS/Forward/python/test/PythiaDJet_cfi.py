import FWCore.ParameterSet.Config as cms

from GeneratorInterface.Pythia6Interface.pythiaDefault_cff import *
generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    comEnergy = cms.double(10000.0),

    PythiaParameters = cms.PSet(
        # Default (mostly empty - to keep PYTHIA default) card file
        # Name of the set is "pythiaDefault"
        pythiaDefaultBlock,
        # User cards - name is "myParameters"
        # Pythia's random generator initialization 
        myParameters = cms.vstring(),
        # This is a vector of ParameterSet names to be read, in this order
        # The first two are in the include files below
        # The last one are simply my additional parameters
        parameterSets = cms.vstring('pythiaDefault', 
            'pythiaMinBias', 
            'myParameters'),

        pythiaMinBias = cms.vstring('MSEL=0         ! User defined processes (full user control): select events involving only quarks', 
            'MSUB(11)=1     ! qq->qq, qq"->qq" ', 
            'MSUB(12)=1     ! qq(bar)->qq(bar), qq(bar)->q"q"(bar), qq"(bar)->qq"(bar) ', 
            'MSTJ(11)=3     ! Choice of the fragmentation function', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10.   ! for which ctau  10 mm', 
            'MSTP(2)=2      ! which order running alphaS', 
            'MSTP(33)=3     ! no K factors in hard cross sections', 
            'MSTP(51)=7     ! choice of proton parton-distribution set (D=7 and means CTEQ 5L)',
            'MSTP(52)=1     ! choice of proton pdf library (D=1 and means internal pythia one, according to MSTP(51) above',  
            'MSTP(61)=0     ! master switch for initial-state QCD and QED radiation 1=on=default', 
            'MSTP(81)=0     ! multiple parton interactions 1 (= on) is Pythia default, 0 = off', 
            'MSTP(82)=4     ! Defines the multi-parton model', 
            'MSTU(21)=1     ! Check on possible errors during program execution', 
            'PARP(82)=1.9409   ! pt cutoff for multiparton interactions', 
            'PARP(89)=1960. ! sqrts for which PARP82 is set', 
            'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter', 
            'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter', 
            'PARP(90)=0.16  ! Multiple interactions: rescaling power', 
            'PARP(67)=2.5   ! amount of initial-state radiation', 
            'PARP(85)=1.0   ! gluon prod. mechanism in MI', 
            'PARP(86)=1.0   ! gluon prod. mechanism in MI', 
            'PARP(62)=1.25  ! ', 
            'PARP(64)=0.2   ! ', 
            'MSTP(91)=1     !', 
            'PARP(91)=2.1   ! kt distribution', 
            'PARP(93)=15.0  ! ', 
            'CKIN(3)=3.     ! Outgoing parton min Pt (c.m. Frame)', 
            'CKIN(4)=15.    ! Outgoing parton max Pt (c.m. Frame)', 
            'CKIN(9)=5.0    ! Outgoing parton (with max rapidity) min rapidity (c.m. Frame)', 
            'CKIN(10)=6.5   ! Outgoing parton (with max rapidity) max rapidity (c.m. Frame)', 
            'CKIN(11)=5.0   ! Outgoing parton (with min rapidity) min rapidity (c.m. Frame)', 
            'CKIN(12)=6.5   ! Outgoing parton (with min rapidity) max rapidity (c.m. Frame)')
    )
)


