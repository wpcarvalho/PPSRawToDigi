import FWCore.ParameterSet.Config as cms

def removeFromPaths(process,module):
    for pathName in process.paths:
        getattr(process,pathName).remove(getattr(process,module))

from UABaseTree_forward_cfg import process
process.load("UATree.UABaseTree.UABaseTree_forward_MC_cfi")
process.load('UATree.UABaseTree.UABaseTree_MC_cfi')
process.myTTRHBuilderWithoutAngle4PixelTriplets.ComputeCoarseLocalPositionFromDisk = True

removeFromPaths(process,'hltPhysicsDeclared')
#removeFromPaths(process,'noscraping')
removeFromPaths(process,'UEAnalysisTracks')
removeFromPaths(process,'ueSisCone5TracksJet500')
removeFromPaths(process,'UEAnalysisJetsAkOnlyReco')
