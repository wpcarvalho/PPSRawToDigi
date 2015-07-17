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
pool_output_dir = '/castor/cern.ch/totem/offline/mcdata/QGSJET2-01_8TeV_REL445SimAndRec-R-WindowFlange/'
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
process.T2RecoAn2.OutputFile = 'QGSJET2-01_8TeV_REL445SimAndRec-R-WindowFlange.root'
process.p1 = cms.Path(process.T2RecoAn2)
