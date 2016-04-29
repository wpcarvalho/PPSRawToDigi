import FWCore.ParameterSet.Config as cms

# raw data source
from TotemRawData.Readers.TotemStandaloneRawDataSource_cfi import *
source.verbosity = 0
source.printProgressFrequency = 0
source.fileNames.append('/afs/cern.ch/user/j/jkaspar/public/run_9987_EVB11_1.003.srs')

# raw-to-digi conversion
from CondFormats.TotemReadoutObjects.TotemDAQMappingESSourceXML_cfi import *
TotemDAQMappingESSourceXML.mappingFileNames.append("CondFormats/TotemReadoutObjects/xml/totem_rp_210far_220_mapping.xml")

from EventFilter.TotemRawToDigi.totemTriggerRawToDigi_cfi import *
totemTriggerRawToDigi.rawDataTag = cms.InputTag("source")
totemTriggerRawToDigi.fedId = 0x29c

from EventFilter.TotemRawToDigi.totemRPRawToDigi_cfi import *
totemRPRawToDigi.rawDataTag = cms.InputTag("source")
totemRPRawToDigi.fedIds = cms.vuint32(0x1a1, 0x1a2, 0x1a9, 0x1aa, 0x1b5, 0x1bd)

totemRawToDigi = cms.Sequence(
  totemTriggerRawToDigi *
  totemRPRawToDigi
)
