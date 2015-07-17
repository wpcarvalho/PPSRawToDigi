import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("Phojet",
    verbosity = cms.untracked.uint32(12),
    bufferSize = cms.untracked.uint32(10), # number of particles to generate in one run of phojet
    cmsEnergy = cms.double( 8000.0 ), # center mass system energy
    process = cms.string('DPE'), # process, from the list (DPE, DD, SD, MB) or specified manually (ex. "0 0 0 0 1 1 0 1")
    phojetExecutable = cms.string('IOMC/Phojet/data/main1') # compiled phojet executable
)

# phojet process specification "1 2 3 4 5 6 7 8", where 
#   1  non-diffractive inelastic process
#   2  elastic scattering
#   3  quasi-elastic rho/omega/phi and pi+/pi-production
#   4  central diffraction (DPE)
#   5  single diffraction of PARTICLE1  
#   6  single diffraction  of PARTICLE2
#   7  double diffraction
#   8  single-resolved / direct processes / photon-hadron, photon-photon  
