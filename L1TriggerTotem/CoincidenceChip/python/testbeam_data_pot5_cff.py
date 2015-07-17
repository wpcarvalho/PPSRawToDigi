import FWCore.ParameterSet.Config as cms

# raw data source
from TotemRawData.Readers.RawDataSource_cfi import *
source.verbosity = 3
source.fileNames.append('/data2/totem/RPtest/data/pot5/rp5_cosmic12.vme')
source.fileNames.append('/data2/totem/RPtest/data/pot5/rp5_cosmic13.vme')
source.fileNames.append('/data2/totem/RPtest/data/pot5/rp5_cosmic14.vme')

# raw to digi conversion
DAQInformationSourceXML = cms.ESSource("DAQInformationSourceXML",
    xmlFileName = cms.string('/data2/totem/RPtest/xml/POT5_planes.xml')
)

RawToDigi = cms.EDProducer("RPDataDigiProducer",
	verbosity = cms.untracked.uint32(0)
)

#RawToCC = cms.EDProducer("RPDataCCProducer",
#	verbosity = cms.untracked.uint32(5),
#    productLabelRaw = cms.string('')
#    #productLabelRaw = cms.string('RPCCRawBits')
#)

from TotemRawData.RawToDigi.RPDataCCProducer_cfi import *
RawToCC           = RPDataCCProducer
RawToCC.verbosity = 5
#RawToCC.productLabelRaw = RPDataCCProducer.productLabelRaw


