import FWCore.ParameterSet.Config as cms

# raw data source
source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/j/jkaspar/public/run268608_ls0001_streamA_StorageManager.root')
)

# raw-to-digi conversion
from CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi import *
TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/ctpps_210_mapping.xml")

# in the emulated data the trigger block contains non-sense
#from EventFilter.TotemRawToDigi.totemTriggerRawToDigi_cfi import *
#totemTriggerRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")
#totemTriggerRawToDigi.fedId = 577

from EventFilter.TotemRawToDigi.totemRPRawToDigi_cfi import *
totemRPRawToDigi.rawDataTag = cms.InputTag("rawDataCollector")
totemRPRawToDigi.fedIds = cms.vuint32(578, 579, 580) # in the emulated data one OptoRx was not functional
totemRPRawToDigi.RawToDigi.testID = 0
totemRPRawToDigi.RawToDigi.printErrorSummary = 1
totemRPRawToDigi.RawToDigi.printUnknownFrameSummary = 1

totemRawToDigi = cms.Sequence(
  #totemTriggerRawToDigi *
  totemRPRawToDigi
)
