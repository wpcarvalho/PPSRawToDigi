import FWCore.ParameterSet.Config as cms
import os
#import sys

#process.load("L1TriggerTotem.CoincidenceChip.testbeam_data_cff")
from L1TriggerTotem.CoincidenceChip.runs_data_cff import sourceFile

import Configuration.Info.job_info as job_info
job_info.getInfo()

#################################################################
from Configuration.TotemCommon.totemCLParser import parser

parser.add_argument(
    "+run","++runNumber",
    dest    = 'runNumber',
    default = "run_9_A",
    action  = "store",
    choices =  [str(key) for key in sorted(sourceFile.keys())],
    help    = "Roman Pot under test"
)

parser.add_argument(
    "+nplanes",
    dest    = 'nplanes',
    default = 3,
    type    = int,
    action  = "store",
    help    = "minimum number of planes needed for road finding and track fitting"
)
parser.add_argument(
    "+chosenPotsId",
    dest    = 'chosenPotsId',
    nargs = '*',
    default = [],
    type    = int,
    action  = "store",
    help    = "chose RPs for further analysis by specifying their ID if no RP is specified all RPs are chosen"
)

parser.add_argument(
    "+outputdir",
    dest    = 'outputDir',
#    nargs = '*',
    default = "trigger_studies_def",
#    type    = int,
    action  = "store",
    help    = "directory where the output will be stored"
)

parser.add_argument(
    "+tracks",
    dest    = 'tracks',
#    nargs = '*',
#    default = [],
#    type    = bool,
    action  = "store_true",
    help    = "do not skip track finding and fitting"
)

parser.add_argument(
    "+reco",
    dest    = 'reco',
#    nargs = '*',
#    default = [],
#    type    = bool,
    action  = "store_true",
    help    = "make reconstruction"
)

parser.add_argument(
    "+trigger",
    dest    = 'trigger',
#    nargs = '*',
    default = "",
#    type    = bool,
    choices =  ["","UorV", "UandV", "UorV1sectorOnly","UandV1sectorOnly","False"],
    action  = "store",
    help    = "Definition of sending trigger signal. If trigger=\"\" then the value has to by specified somewhere in a configuration file. If trigger!=\"\" then this trigger is used."
)

parser.set_defaults(
    nevents   = -1  ,
#    energy    = 3500,
#    beta      = "",
    verbosity = 1
)
##################################################################
args = parser.parse_args()
try:
  args.runNumber = int(args.runNumber)
finally:
  pass
##################################################################
triggerFromFile = sourceFile[args.runNumber].params.trigger
#print triggerFromFile.configValue
#print triggerFromFile.pythonValue
if args.trigger == "":
  if triggerFromFile == "":
    print("error: no trigger defined!")
    assert(0)
  else:
#    print("overwriting trigger def")
    args.trigger = triggerFromFile.value()
print("args.trigger="+args.trigger)
##################################################################
_chosenPotsFromFile = sourceFile[args.runNumber].params.chosenPotsId
#print triggerFromFile.configValue
#print triggerFromFile.pythonValue
if len(args.chosenPotsId) == 0:
  if len(_chosenPotsFromFile):
    args.chosenPotsId = _chosenPotsFromFile.value()
print("args.chosenPotsId = ",args.chosenPotsId)
##################################################################

preprocessName="triggerAnalyzerProcess"
if args.reco:
  process = cms.Process(preprocessName+"Reco")
else:
  process = cms.Process(preprocessName)

#if args.printRuns == True:
  #for sf in sourceFile:
#  print(sourceFile.keys())
  #  print sf.keys()
#  sys.exit(0)

#process.load("Configuration.TotemCommon.LoggerMax_cfi")
process.load("Configuration.TotemCommon.LoggerMin_cfi")

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(args.nevents)
)

# Geometry - we assume beta*=90m for analysis of good tracks
process.load("Configuration.TotemCommon.geometryRP_cfi")
#process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_'+args.beta+'/RP_Dist_Beam_Cent.xml')
process.XMLIdealGeometryESSource.geomXMLFiles = sourceFile[args.runNumber].params.geomXMLFiles


process.DAQInformationSourceXML = cms.ESSource("DAQInformationSourceXML",
  xmlFileName = cms.string('')
)
process.DAQInformationSourceXML.xmlFileName = sourceFile[args.runNumber].params.xmlFileName


process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")
# process.RPClustProd.Verbosity = 1

process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")
# process.RPHecoHitProd.Verbosity = 1

# takes raw data and calculates trigger bits which may serve as input to CC
process.TriggerBits = cms.EDProducer("RPTriggerBitsProducer",
  verbose = cms.bool(False)
)

# raw to digi conversion
process.load('TotemRawData.RawToDigi.RPDataDigiProducer_cfi')
#process.RPDataDigiProducer.verbosity = 0

#process.RawToDigi = cms.EDProducer("RPDataDigiProducer",
#	verbosity = cms.untracked.uint32(0)
#)

process.load('TotemRawData.RawToDigi.RPDataCCProducer_cfi')
process.RawToCC = process.RPDataCCProducer
process.RawToCC.verbosity = 5 

process.load("L1TriggerTotem.CoincidenceChip.RPCoincidenceProducer_cfi")
process.RPCC.coincidenceChipConfig = sourceFile[args.runNumber].params.coincidenceChipConfig 

#process.load("L1TriggerTotem.CoincidenceChip.RPNonParallelTrackCandidateFinder2_cfi")
process.load("RecoTotemRP.RPNonParallelTrackCandidateFinder.RPNonParallelTrackCandidateFinder_cfi")
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = args.nplanes
process.NonParallelTrackFinder.minPlanesPerProjectionToFit    = args.nplanes

#process.load("RecoTotemRP.RPSingleCandidateTrackFinder.RPSingleTrackCandFindConf_cfi")
# process.RPSinglTrackCandFind.Verbosity = 1

process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
# process.RPSingleTrackCandCollFit.Verbosity = 1

process.RPTriggerAnalyzer = cms.EDAnalyzer("RPTriggerAnalyzer",
  verbosity             = cms.uint32(args.verbosity),
  planes                = cms.uint32(args.nplanes),
  outputFile            = cms.string('run_'+str(args.runNumber)+'_triggerStats.root'),
  outputDir             = cms.string(args.outputDir+'/run_'+str(args.runNumber)),
  chosenPotsId          = cms.vuint32(args.chosenPotsId),
  # tagRaw = cms.InputTag("RawToCC"),
  # tagSimu = cms.InputTag("RPCC"),
  modulLabelRaw         = cms.string('RawToCC'),
  productLabelRaw       = cms.string(''),
  modulLabelSimu        = cms.string('RPCC'),
  productLabelSimu      = cms.string(''),
  tracks                = cms.bool(args.tracks), # should we produce some statistics for tracks & trigger?
  trigger               = cms.string(args.trigger),
  params                = sourceFile[args.runNumber].params
)

#process.RPTriggerAnalyzer.tagSimu.setProductInstanceLabel(str(process.RPCC.productLabelSimu))
#process.RPTriggerAnalyzer.tagRaw.setProductInstanceLabel(str(process.RawToCC.productLabelRaw))

process.RPTriggerAnalyzer.productLabelSimu = process.RPCC.productLabelSimu
process.RPTriggerAnalyzer.productLabelRaw  = process.RawToCC.productLabelRaw

recoTrackSequence = cms.Sequence(process.NonParallelTrackFinder*process.RPSingleTrackCandCollFit)
digiSequence = cms.Sequence(process.RPDataDigiProducer*process.RawToCC*process.TriggerBits*process.RPClustProd*process.RPHecoHitProd*process.RPCC)

if args.tracks:
  finalSequence = cms.Sequence(digiSequence*recoTrackSequence)
else:
  finalSequence = cms.Sequence(digiSequence)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

def processRawSource(process):
  #############################
  # raw data source
  process.load('TotemRawData.Readers.RawDataSource_cfi')
  process.source.verbosity = 3
  process.source.skipCorruptedEvents = True
  process.source.printProgressFrequency = 0
  process.source.skipCorruptedEvents = False
  process.source.fileNames   = sourceFile[args.runNumber].params.sourceFileNames


recoFileName = './reco/run_'+str(args.runNumber)+'_reco_triggerAnalyzer_.root'
os.system("mkdir -p `dirname "+recoFileName+"`")

def recoOutput(process):
  process.output = cms.OutputModule("PoolOutputModule",
      fileName = cms.untracked.string('file:'+recoFileName),
      outputCommands = cms.untracked.vstring(
           'keep *',
           'drop TotemRawEvent_*_*_*'
  #        'drop *',
  #        'keep RPDigClusteredmDetSetVector_*_*_*',
  #        'keep RPFittedTrackCollection_*_*_*',
  #        'keep RPRecoHitedmDetSetVector_*_*_*',
  #        'keep RPRecognizedPatternsCollection_*_*_*',
  #        'keep RPStripDigiedmDetSetVector_*_*_*',
  #        'keep RPTrackCandidateCollection_*_*_*',
     )
  )
  process.outpath = cms.EndPath(process.output)

if args.reco:
  processRawSource(process)
  recoOutput(process)
  # process.p = cms.Path(finalSequence*process.dump)
  process.p = cms.Path(finalSequence)
else:
  if os.path.exists(recoFileName):
    print("using reco file: "+recoFileName)
    #############################
    process.source = cms.Source("PoolSource",
      verbosity = cms.untracked.uint32(0),
      fileNames = cms.untracked.vstring(),
      noEventSort = cms.untracked.bool(True),
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    )
    process.source.fileNames.append('file:'+recoFileName)
    process.p1 = cms.Path(process.RPTriggerAnalyzer)
  else:
    processRawSource(process)
    #recoOutput(process)
    process.p1 = cms.Path(finalSequence*process.RPTriggerAnalyzer)

#process.p1 = cms.Path(finalSequence*process.dump*process.RPTriggerAnalyzer)
#process.p1 = cms.Path(finalSequence*process.RPTriggerAnalyzer)

#process.p1 = cms.Path(process.RawToDigi*process.RawToCC*process.TriggerBits*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit*process.RPCC)

#process.o1 = cms.OutputModule("PoolOutputModule",
#  fileName = cms.untracked.string('file:output2.root')
#)
#process.outpath = cms.EndPath(process.o1)
####################################################################################
