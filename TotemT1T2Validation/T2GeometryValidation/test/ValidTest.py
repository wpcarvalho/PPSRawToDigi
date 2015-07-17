import FWCore.ParameterSet.Config as cms

process = cms.Process("valRPT1T2pythiaSDbeta2")

########################### COMMON PART ##########################################

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/home/mirko/SL/WorkingArea/CMSSW311_TEST/CMSSW_3_1_1/src/TotemAlignment/T2TrkBasedInternalAlignment/test/VisualizationSimu/10MuPlNear.root')
                            #('file:/media/usbdisk/RecoBeforeHit2_Merged_OptorRX1_4And6clk.root')
)

# logging to txt files 
process.load("Configuration.TotemCommon.LoggerMin_cfi")


process.DAQInformationSourceXML = cms.ESSource("DAQInformationSourceXML",
    xmlFileName = cms.string('TotemT1T2Validation/T2GeometryValidation/test/T2GeoMapIP5_4quarter_vmea_cmssw.xml')
                                                #T2MapIP5_D3_D2_NewFormatBIS.xml
)

########################### T2 VALIDATION ##########################################


############# T2 RECO ###########

#process.load("TotemT1T2Validation.T2RecoValidation.T2RecoValidation_cfi")
#process.T2RecoAn.OutputFile = 'file:10MuPlNear.root'
#process.T2RecoAn.singleparticle = False     # true only for particle gun
#process.T2RecoAn.FitdzEnergy = 50.0         # use 30.0 or 50.0

process.load("TotemT1T2Validation.T2GeometryValidation.T2GeoValidation2_cfi")
process.T2GeoAn2.OutputFile = 'file:10MuPlNear.root'
process.T2GeoAn2.singleparticle = False    # true only for particle gun
process.T2GeoAn2.FitdzEnergy = 50.0        # use 30.0 or 50.0
process.T2GeoAn2.DoClusterFromVFatChannel = False
#process.p1 = cms.Path(process.RPInelProtRecVal*process.RPRecTracksVal*process.RPHitDists*process.TransportValidation*process.T1Val*process.T2G4An*process.T2DigiAn*process.T2RecoAn*process.T2GeoAn2)
process.p1 = cms.Path(process.T2GeoAn2)
