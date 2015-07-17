import FWCore.ParameterSet.Config as cms

process = cms.Process("valT1T2pi")

########################### COMMON PART ##########################################

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/tmp/berretti/run_prodT2_Pythia8Inelastic_E3.5TeV_ONLYGEANT_Nov2011_91mmSamp4_RELEASE424-95.root__startEvt0.root')
# )
#QGSET_8TeV_Test_REL445SimAndRec-I Py8_8TeV_Test_REL445SimAndRec-L
#pool_output_dir ='/castor/cern.ch/user/b/berretti/QGSJET2-01_8TeV_REL445Pure_SimAndRec-S-NoRECO_Vtx11250_445PureRecoT1T2GhostSuppressionAlignV1-GeomEffiSTD7'
fname='QGSJET2-01_8TeV_REL445Pure_SimAndRec-T-NoRECO_Vtx0_Sample0_445PureRecoT1T2GhostSuppressionAlignV1-Perf_Effi'
fname='Py8_TuneCMS_8TeV_REL445Pure_SimAndRec-T-NoRECO_Vtx0_Sample0_445PureRecoT1T2GhostSuppressionAlignV1-Perf_EffiNumPadTrackerSTD'
fname='EPOS_LHCRetune_8TeV-Vtx0-Pure445-CorrectCalo-LargeSamp0_445PureRecoT1T2GhostSuppressionAlignV2-Tuned8372Effi'

pool_output_dir ='/castor/cern.ch/user/b/berretti/'+fname
finalCastordir = pool_output_dir 

commandls = 'nsls ' + finalCastordir + ' | grep .root | grep -v .py' 
print '@@@@@@@@@@@@@@@ \n'
#echo $finalCastordir
import commands
stat, out = commands.getstatusoutput(commandls)
if not stat:
    print out


#Here I split the list of files

import string
lines = string.split(out, '\n')

new_list = []

for oneline in lines:
 #print 'rfio:', oneline
 OnecastorFilePath = 'rfio:'  + finalCastordir + '/' + oneline
 new_list += [OnecastorFilePath]

process.source = cms.Source("PoolSource",
 
    fileNames = cms.untracked.vstring(new_list),                        
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')                        
)


process.load("TotemT1T2Validation.T2RecoValidation.T2RecoValidation2_cfi")
process.T2RecoAn2.OutputFile = fname+'.root'
process.p1 = cms.Path(process.T2RecoAn2)
