import FWCore.ParameterSet.Config as cms

process = cms.Process("valT1T2pi")

########################### COMMON PART ##########################################

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#pool_output_dir = '/castor/cern.ch/user/b/berretti/prodT2_Pythia8Inelastic_E3.5TeV_ONLYGEANT_Nov2011_91mmSamp3_NEWRELEASE_Samp0_GhostSuppression_BB_4pad_12DecGENRECO_T1T2Reco/'  

pool_output_dir = '/castor/cern.ch/user/b/berretti/prodT2_PhojetInelastic_E3.5TeV_NEWRELEASE_Samp0XSect_NoCHPartInT2_OnlyGEN/'

pool_output_dir = '/castor/cern.ch/totem/offline/mcdata/Py8_GenOnlyND/'


process.load("TotemT1T2Validation.GeneratorAnalysis.T2Generator_cfi")
process.T2Gen.OutputFile = 'H2ON.root'
process.T2Gen.FullInfos = cms.bool(False)





finalCastordir = pool_output_dir 

commandls = 'nsls ' + finalCastordir + ' | grep .root | grep -v .py' 

print (commandls)
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

#print '@@@@@@@@@@@@@@@ \n'

#for oneline2 in new_list:
# print oneline2


#Here I feed the cmssw process
#/castor/cern.ch/user/j/jwelti/sherpamirko.root 
process.source = cms.Source("PoolSource",
   # fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/j/jwelti/sherpajan.root'),#new_list
    fileNames = cms.untracked.vstring(new_list),                        
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')                        
)





# logging to txt files
process.load("Configuration.TotemCommon.LoggerMax_cfi")

# random number generators seeds
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

########################### T2 VALIDATION ##########################################
#/castor/cern.ch/totem/offline/mcdata/



process.p1 = cms.Path(process.T2Gen)
