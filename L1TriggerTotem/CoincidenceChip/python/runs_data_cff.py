import FWCore.ParameterSet.Config as cms

"""
import sys
sys.path.insert(0,'/home/jprochaz/workspace/cmssw/311/src/L1TriggerTotem/CoincidenceChip/pymodules/common')
sys.path.append('./pymodules/common')
sys.path.append('$(CMSSW_BASE)/src/L1TriggerTotem/CoincidenceChip/pymodules/common')
sys.path.append('./pymodules')
print sys.path
import collections
import bidict
"""

rpMappingCCUtoCMSSW = {"7a":122,"7c":120,"7f":123,"7b":121,"7d":124,"7e":125,"6e":25,"6f":23,"6a":22,"6b":21,"6c":20,"6d":24}
rpMappingH8toCCU = {1:"6f",2:"7f",3:"7e",4:"7d", 5:"7c",6:"7b", 7:"7a",8:"6a",9:"6d",10:"6e",11:"6c",12:"6b"}
rpMappingCMSSWNumberCMSSWName = {
    20:"rp_45_220_nr_tp",
    21:"rp_45_220_nr_bt",
    22:"rp_45_220_nr_hr",
    23:"rp_45_220_fr_hr",
    24:"rp_45_220_fr_tp",
    25:"rp_45_220_fr_bt",
    120:"rp_56_220_nr_tp",
    121:"rp_56_220_nr_bt",
    122:"rp_56_220_nr_hr",
    123:"rp_56_220_fr_hr",
    124:"rp_56_220_fr_tp",
    125:"rp_56_220_fr_bt"
}

horizontalPots = ['rp_45_220_fr_hr','rp_45_220_nr_hr','rp_56_220_fr_hr','rp_56_220_nr_hr']

class MyPSet(cms.PSet):
  def __init__(self, *args, **kwargs):
    #print args
    #print kwargs
    # add these parameters 
    mydict = dict(
          xmlFileName = cms.string(''), 
          geomXMLFiles = cms.vstring(),
          sourceFileNames   = cms.untracked.vstring(),
          chosenPotsId   = cms.untracked.vuint32(),
          runInfo   = cms.string(''),
          runName   = cms.string(''),
          trigger = cms.string(''),
          maskCC = cms.vstring()
    )

    newkwargs = dict(mydict,**kwargs)
    #print newkwargs
    cms.PSet.__init__(self, *args, **newkwargs)



class SourceFile:
    def __init__(self,params):
       # cloning of params is esential. Each SourceFile has its own params. 
       self.params = params.clone()

#    def setCCcontrolRegisters(self, reg1, reg2, reg3):
#       self.ccControlRegister1 = reg1
#       self.ccControlRegister2 = reg2
#       self.ccControlRegister3 = reg3

sourceFile = {}

testBeamParams = MyPSet(
  #useValidTrackCondition = cms.bool(True),
  coincidenceChipConfig = cms.PSet(
    useControlRegisters = cms.bool(True),
    controlRegister1    = cms.uint32(18),
    controlRegister2    = cms.uint32(240),
    controlRegister3    = cms.uint32(136),
    useLogicWithWrongNP = cms.uint32(0)
  ),
  trigger = cms.string("UorV")
)

runs_testBeamH8 = []
runLabel = "run_9_A"
runs_testBeamH8.append(runLabel)
sourceFile[runLabel] = SourceFile(testBeamParams)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.xmlFileName = '/data2/totem/RPtest2009/xml/POT9_planes.xml' #'/data2/totem/RPtest/xml/POT9_planes.xml'
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot9/pot9_2ndtest_beam_2.vme') #/data2/totem/RPtest/data/pot9/pot9_2ndtest_beam_2.vme')
sourceFile[runLabel].params.runInfo = '  CC settings: 3 18 240 136 128\n  TM latency: 33\n'

runLabel = "run_9_B"
runs_testBeamH8.append(runLabel)
sourceFile[runLabel] = SourceFile(testBeamParams)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot9/pot9_2ndtest_beam_3.vme')
sourceFile[runLabel].params.xmlFileName = '/data2/totem/RPtest2009/xml/POT9_planes.xml'
sourceFile[runLabel].params.runInfo = '  CC settings: 3 18 240 136 0\n  TM latency: 34\n' 

runLabel = "run_10_A"
runs_testBeamH8.append(runLabel)
sourceFile[runLabel] = SourceFile(testBeamParams)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot10/pot10_beam_10.vme')
sourceFile[runLabel].params.xmlFileName = '/data2/totem/RPtest2009/xml/POT10_planes.xml'
sourceFile[runLabel].params.runInfo = '  CC settings: 3 18 240 136 128\n  TM latency: 33 (2 clk)\n  VFAT latency: 38\n  VFAT TH= 10\n'

runLabel = "run_10_Bvth15"
runs_testBeamH8.append(runLabel)
sourceFile[runLabel] = SourceFile(testBeamParams)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot10/pot10_beam_1.vme')
sourceFile[runLabel].params.xmlFileName = '/data2/totem/RPtest2009/xml/POT10_planes.xml'
sourceFile[runLabel].params.runInfo = '  CC settings: 3 18 240 136 0\n  TM latency: 34 (2 clk)\n  VFAT latency: 38\n  VFAT TH= 15\n'

runLabel = "run_10_Bvth10"
runs_testBeamH8.append(runLabel)
sourceFile[runLabel] = SourceFile(testBeamParams)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot10/pot10_beam_9.vme')
sourceFile[runLabel].params.xmlFileName = '/data2/totem/RPtest2009/xml/POT10_planes.xml'
sourceFile[runLabel].params.runInfo = '  CC settings: 3 18 240 136 0\n  TM latency: 34 (2 clk)\n  VFAT latency: 38\n  VFAT TH= 10\n'

runLabel = "run_11_A"
runs_testBeamH8.append(runLabel)
sourceFile[runLabel] = SourceFile(testBeamParams)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot11/pot11_beam_1.vme')
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot11/pot11_beam_2.vme')
sourceFile[runLabel].params.xmlFileName = '/data2/totem/RPtest2009/xml/POT11_planes.xml'
sourceFile[runLabel].params.runInfo = '  CC settings: 3 18 240 136 128\n  TM latency=33'

runLabel = "run_12_A"
runs_testBeamH8.append(runLabel)
sourceFile[runLabel] = SourceFile(testBeamParams)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot12/pot12_beam1_005.vme')
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot12/pot12_beam1_006.vme')
sourceFile[runLabel].params.xmlFileName = '/data2/totem/RPtest2009/xml/POT12_planes.xml'
sourceFile[runLabel].params.runInfo = '  CC settings: 3 18 240 136 128\n' 

runLabel = "run_12_B"
runs_testBeamH8.append(runLabel)
sourceFile[runLabel] = SourceFile(testBeamParams)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot12/pot12_beam1_001.vme')
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RPtest2009/data/pot12/pot12_beam1_002.vme')
sourceFile[runLabel].params.xmlFileName = '/data2/totem/RPtest2009/xml/POT12_planes.xml'
sourceFile[runLabel].params.runInfo = '  CC settings: 3 18 240 136 0\n' 

testBeamRunPot = {}
for run in runs_testBeamH8:
  pot = ""
  for char in run.replace("run_",""):
    if char == "_": break
    pot+=char
  testBeamRunPot[run] = pot

testDataIP5Params = MyPSet(
  #useValidTrackCondition = cms.bool(False),
  coincidenceChipConfig = cms.PSet(
    useControlRegisters = cms.bool(True), 
    controlRegister1 = cms.uint32(20),
    controlRegister2 = cms.uint32(240),
    controlRegister3 = cms.uint32(128),
    useLogicWithWrongNP = cms.uint32(0)
  )
)

runLabel = "run_X246"
sourceFile[runLabel] = SourceFile(testDataIP5Params)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RP_IP5/data/run_246.vmea')
sourceFile[runLabel].params.xmlFileName = '/afs/cern.ch/user/j/jprochaz/scratch/jprochaz/cmssw_users/user/jprochaz/mapping/RP_all_new.xml'
#sourceFile[runLabel].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
sourceFile[runLabel].params.runInfo = '  Vfat Latency= 157 4 clk\n  TM lat = 155 \n  run 246 221 ev in 26 min'

runLabel = "run_X247"
sourceFile[runLabel] = SourceFile(testDataIP5Params)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RP_IP5/data/run_247.vmea')
sourceFile[runLabel].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
sourceFile[runLabel].params.runInfo = ' MASK all trigger sectors except  pot 7a (56_nr_hr)\n  Vfat Latency= 157 8 clk\n  TM lat = 155\n  CC: V=3\n  run 247  27 ev in 25 min'

runLabel = "run_X248"
sourceFile[runLabel] = SourceFile(testDataIP5Params)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/RP_IP5/data/run_248.vmea')
sourceFile[runLabel].params.xmlFileName = '/afs/cern.ch/user/j/jprochaz/scratch/jprochaz/cmssw_users/user/jprochaz/mapping/RP_all_new.xml'
#sourceFile[runLabel].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
sourceFile[runLabel].params.runInfo = '  Vfat Latency= 157 4 clk\n  TM lat = 155\n  CC: V=4\n  run 248' 

latencyScanParams = MyPSet(
  #useValidTrackCondition = cms.bool(False),
  coincidenceChipConfig = cms.PSet(
    useControlRegisters = cms.bool(True), 
    controlRegister1 = cms.uint32(19),
    controlRegister2 = cms.uint32(240),
    controlRegister3 = cms.uint32(128),
    useLogicWithWrongNP = cms.uint32(0)
  )
)


"""
vth=10
V/NP=3
n clocks = 4  (vfat, TMC, CC)

trigger  |7a(122)|7c(120)|7f(123)|7b(121)|7d(124)|7e(125)|6e(25)|6f(23)|6a(22)|6b(21)|6c(20)|6d(24)|
----------------------------------------------------------------------------------------------------
lat  153 | 695   |  699  |   -   |  708  |   -  |   -    |  -   |  -   |  744 |  770 |  773 |  -   |
     155 | 693   |  698  |  705  |  706  |  713 |  717   |  -   |  -   |  743 |  768 |  774 | 779  |
     157 |  -    |  697  |  704  |   -   |  710 |  715   | 730  | 736  |   -  |   -  |   -  | 778  |
     159 |  -    |   -   |  702  |   -   |   -  |  719   | 731  | 738  |   -  |   -  |   -  |  -   |
"""
#lat153runs = [695,699,708,744,770,773]
#lat155runs = [693,698,705,706,713,717,743,768,774,779] 
#lat157runs = [697,704,710,715,730,736,778]
#lat159runs = [702,719,731,738]


latRuns = dict([
    (695,(153,"7a")),(693,(155,"7a")),
    (699,(153,"7c")),(698,(155,"7c")),(697,(157,"7c")),
    (705,(155,"7f")),(704,(157,"7f")),(702,(159,"7f")),
    (708,(153,"7b")),(706,(155,"7b")),
    (713,(155,"7d")),(710,(157,"7d")),
    (717,(155,"7e")),(715,(157,"7e")),(719,(159,"7e")),
    (730,(157,"6e")),(731,(159,"6e")),
    (736,(157,"6f")),(738,(159,"6f")),
    (744,(153,"6a")),(743,(155,"6a")),
    (770,(153,"6b")),(768,(155,"6b")),
    (773,(153,"6c")),(774,(155,"6c")),
    (779,(155,"6d")),(778,(157,"6d"))
])

#runs = lat153runs + lat155runs + lat157runs + lat159runs
#latRuns = map(lambda x: "run_"+str(x),runs)

#for run in runs:
for run, latpot in latRuns.iteritems():
  lat,pot = latpot
  #print run, lat, pot
  runLabel = "run_"+str(run)
  sourceFile[runLabel] = SourceFile(latencyScanParams)
  sourceFile[runLabel].params.runName = runLabel
  sourceFile[runLabel].params.sourceFileNames.append('/data2/totem/usr/jprochaz/latency_scan/data/run_'+str(run)+'.vmea')
  #sourceFile[runLabel].params.sourceFileNames.append('vmeastream:///data2/totem/usr/jprochaz/latency_scan/data/run_'+str(run)+'.vmea')
  sourceFile[runLabel].params.xmlFileName = '/afs/cern.ch/user/j/jprochaz/scratch/jprochaz/cmssw_users/user/jprochaz/mapping/RP_all_new.xml'
  #sourceFile[runLabel].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
  #if run in  lat153runs: info = "153"
  #if run in  lat155runs: info = "155"
  #if run in  lat157runs: info = "157"
  #if run in  lat159runs: info = "159"
  sourceFile[runLabel].params.runInfo = ' vth=10 \n V/NP=3 \n n clocks = 4  (vfat, TMC, CC) \n lat '+str(lat)+'\n  triggering on pot='+pot+'('+str(rpMappingCCUtoCMSSW[pot])+')'

"""  #                   |++             |++                                  |++           |++
trigger  |7a(122)|7c(120)|7f(123)|7b(121)|7d(124)|7e(125)|6e(25)|6f(23)|6a(22)|6b(21)|6c(20)|6d(24)|
----------------------------------------------------------------------------------------------------
lat  153 | +695  | -699  |   -   | +708  |   -  |   -    |  -   |  -   | +744 | +770 |+-773 |  -   |
     155 | -693  | +698  | -705  | +706  | +713 | -717   |  -   |  -   | +743 | +768 | +774 |+779  |
     157 |  -    | +697  | +704  |   -   | +710 | +715   | -730 | +-736|   -  |   -  |   -  |+778  |
     159 |  -    |   -   | +702  |   -   |   -  | -719   |-+731 | +-738|   -  |   -  |   -  |  -   |
"""                                        

########### runs june 25 2010 #################### 
# exercise with collimators and RPs
for run in range(1338,1347):
  #print run, lat, pot
  runLabel = "run_"+str(run)
  sourceFile[runLabel] = SourceFile(latencyScanParams)
  sourceFile[runLabel].params.runName = runLabel
  sourceFile[runLabel].params.sourceFileNames.append('file:/afs/cern.ch/exp/totem/scratch/data/RP/2010_06_25/reco/reco_'+str(run)+'.root')
#  sourceFile[runLabel].params.sourceFileNames.append('/afs/cern.ch/exp/totem/scratch/data/RP/2010_06_25/reco/reco_'+str(run)+'.root')
  #sourceFile[runLabel].params.sourceFileNames.append('vmeastream:///data2/totem/usr/jprochaz/latency_scan/data/run_'+str(run)+'.vmea')
  sourceFile[runLabel].params.xmlFileName = '/afs/cern.ch/exp/totem/scratch/data/RP/2010_06_25/xml/RP_all_new.xml'
  sourceFile[runLabel].params.runInfo = "no info"
#  sourceFile[runLabel].params.geomXMLFiles.append("Geometry/TotemRPData/data/2010_06_25/RP_Dist_Beam_Cent.xml")

########### runs july 13 2010 #################### 
# trigger on 1 RP only 
runs_13_06_2010_Params = MyPSet(
  coincidenceChipConfig = cms.PSet(
    useControlRegisters = cms.bool(True), 
    controlRegister1 = cms.uint32(19),
    controlRegister2 = cms.uint32(240),
    controlRegister3 = cms.uint32(128),
    useLogicWithWrongNP = cms.uint32(0)
  ),
  trigger = cms.string("UandV1sectorOnly")
)

runs_13_06_2010_Params_triggerOnPot = {
  1827:24, 1828:25, 1829:23, 1830:21, 1831:22, 1832:20, 1834:20,
  1833:125, 1835:123, 1836:124, 1837:121, 1838:122, 1839:120
}

#runs_13_06_2010 = range(1827,1839+1)
runs_13_06_2010 = runs_13_06_2010_Params_triggerOnPot.keys()
for runNumber in runs_13_06_2010:
  sourceFile[runNumber] = SourceFile(runs_13_06_2010_Params)
  sourceFile[runNumber].params.runName = str(runNumber)
  sourceFile[runNumber].params.sourceFileNames.append('/castor/cern.ch/totem/LHCRawData/2010/RP/run_'+str(runNumber)+'.vmea')
  sourceFile[runNumber].params.xmlFileName = '/afs/cern.ch/user/j/jprochaz/scratch/jprochaz/cmssw_users/user/jprochaz/mapping/RP_all_new.xml'
  sourceFile[runNumber].params.runInfo = "no info"
#  sourceFile[runLabel].params.trigger = "UandV1sectorOnly"

#######################################################
# Tue Jul 13 14:55:45 2010 https://www.ba.infn.it/elog/TOTEM-IP5/173 (174)
"""
mask CC in horizontals pots , V=3
trigger all

run 1941

increased hit multiplicity

rate 160 hz

------------------------------------------------------------------
mask CC in horizontals pots , V=5
trigger all

run 1942

rate 120hz 

90k events
"""

run1941Params = MyPSet(
  coincidenceChipConfig = cms.PSet(
    useControlRegisters = cms.bool(True), 
    controlRegister1 = cms.uint32(19),
    controlRegister2 = cms.uint32(240),
    controlRegister3 = cms.uint32(128),
    useLogicWithWrongNP = cms.uint32(0)
  ),
  trigger = cms.string("False"),
  maskCC = cms.vstring(horizontalPots)
)

runLabel = "run_1941"
sourceFile[runLabel] = SourceFile(run1941Params)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/castor/cern.ch/totem/LHCRawData/2010/RP/'+runLabel+".000"+'.vmea')
sourceFile[runLabel].params.xmlFileName = '/afs/cern.ch/user/j/jprochaz/scratch/jprochaz/cmssw_users/user/jprochaz/mapping/RP_all_new.xml'
#sourceFile[runLabel].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
sourceFile[runLabel].params.runInfo = ''

run1942Params = MyPSet(
  coincidenceChipConfig = cms.PSet(
    useControlRegisters = cms.bool(True), 
    controlRegister1 = cms.uint32(21),
    controlRegister2 = cms.uint32(240),
    controlRegister3 = cms.uint32(128),
    useLogicWithWrongNP = cms.uint32(0)
  ),
  trigger = cms.string("False"),
  maskCC = cms.vstring(horizontalPots)
)

runLabel = "run_1942"
sourceFile[runLabel] = SourceFile(run1942Params)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/castor/cern.ch/totem/LHCRawData/2010/RP/'+runLabel+".000"+'.vmea')
sourceFile[runLabel].params.sourceFileNames.append('/castor/cern.ch/totem/LHCRawData/2010/RP/'+runLabel+".001"+'.vmea')
sourceFile[runLabel].params.xmlFileName = '/afs/cern.ch/user/j/jprochaz/scratch/jprochaz/cmssw_users/user/jprochaz/mapping/RP_all_new.xml'
#sourceFile[runLabel].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
sourceFile[runLabel].params.runInfo = ''
#################################################################
#trigger : 45 or 56 / U or V
run1845Params = MyPSet(
  coincidenceChipConfig = cms.PSet(
    useControlRegisters = cms.bool(True), 
    controlRegister1 = cms.uint32(19),
    controlRegister2 = cms.uint32(240),
    controlRegister3 = cms.uint32(128),
    useLogicWithWrongNP = cms.uint32(0)
  ),
  trigger = cms.string("False")
)

runLabel = "run_1845"
sourceFile[runLabel] = SourceFile(run1845Params)
sourceFile[runLabel].params.runName = runLabel
sourceFile[runLabel].params.sourceFileNames.append('/castor/cern.ch/totem/LHCRawData/2010/RP/'+runLabel+".000"+'.vmea')
sourceFile[runLabel].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
sourceFile[runLabel].params.runInfo = ''

#trigger 45 or 56 / U & V
run1846Params = MyPSet(
  coincidenceChipConfig = cms.PSet(
    useControlRegisters = cms.bool(True), 
    controlRegister1 = cms.uint32(19),
    controlRegister2 = cms.uint32(240),
    controlRegister3 = cms.uint32(128),
    useLogicWithWrongNP = cms.uint32(0)
  ),
  trigger = cms.string("False")
)

runLabel = "run_1846"
sourceFile[runLabel] = SourceFile(run1846Params)
sourceFile[runLabel].params.runName = runLabel
for i in range(3):
  sourceFile[runLabel].params.sourceFileNames.append('/castor/cern.ch/totem/LHCRawData/2010/RP/'+runLabel+".00"+str(i)+'.vmea')
sourceFile[runLabel].params.xmlFileName = '/afs/cern.ch/user/j/jprochaz/scratch/jprochaz/cmssw_users/user/jprochaz/mapping/RP_all_new.xml'
#sourceFile[runLabel].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
sourceFile[runLabel].params.runInfo = ''

runLabel = "run_1847"
sourceFile[runLabel] = SourceFile(run1846Params)
sourceFile[runLabel].params.runName = runLabel
for i in range(6):
  sourceFile[runLabel].params.sourceFileNames.append('/castor/cern.ch/totem/LHCRawData/2010/RP/'+runLabel+".00"+str(i)+'.vmea')
sourceFile[runLabel].params.xmlFileName = '/afs/cern.ch/user/j/jprochaz/scratch/jprochaz/cmssw_users/user/jprochaz/mapping/RP_all_new.xml'
#sourceFile[runLabel].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
sourceFile[runLabel].params.runInfo = ''

###########################################################################
# https://www.ba.infn.it/elog/TOTEM-IP5/189
# Fri Jul 30 05:37:29 2010; trigger test
# trigger on 1 RP only  Trigger only 1 sector

runs_30_06_2010_Params = MyPSet(
  coincidenceChipConfig = cms.PSet(
    useControlRegisters = cms.bool(True), 
    controlRegister1 = cms.uint32(19),
    controlRegister2 = cms.uint32(240),
    controlRegister3 = cms.uint32(128),
    useLogicWithWrongNP = cms.uint32(0)
  ),
  trigger = cms.string("UandV1sectorOnly")
)

runs_30_06_2010_Params_triggerOnPot = {
  2200:25, 2201:23, 2202:24, 2203:21, 2204:22, 2205:20, 2206:125,
  2207:123, 2208:124, 2209:121, 2210:122, 2211:120
}

"""
    20:"rp_45_220_nr_tp",
    21:"rp_45_220_nr_bt",
    22:"rp_45_220_nr_hr",
    23:"rp_45_220_fr_hr",
    24:"rp_45_220_fr_tp",
    25:"rp_45_220_fr_bt",
    120:"rp_56_220_nr_tp",
    121:"rp_56_220_nr_bt",
    122:"rp_56_220_nr_hr",
    123:"rp_56_220_fr_hr",
    124:"rp_56_220_fr_tp",
    125:"rp_56_220_fr_bt"
45_far_bt run 2200
45_far_hr run 2201
45_far_tp run 2202
45_nr_bt run 2203
45_nr_hr run 2204
45_nr_tp run 2205
56_fr_bt run 2206
56_fr_hr run 2207
56_fr_tp run 2208
56_nr_bt run 2209
56_nr_hr run 2210
56_nr_tp run 2211
"""

runs_30_06_2010 = range(2200,2211+1)
for runNumber in runs_30_06_2010:
  #runLabel = "run_"+str(runNumber)
  sourceFile[runNumber] = SourceFile(runs_30_06_2010_Params)
  sourceFile[runNumber].params.runName = str(runNumber)
  sourceFile[runNumber].params.sourceFileNames.append('/castor/cern.ch/totem/LHCRawData/2010/RP/run_'+str(runNumber)+'.000.vmea')
  sourceFile[runNumber].params.xmlFileName = '/afs/cern.ch/user/j/jprochaz/scratch/jprochaz/cmssw_users/user/jprochaz/mapping/RP_all_new.xml'
  #sourceFile[runNumber].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
#  sourceFile[runNumber].params.xmlFileName = '/afs/cern.ch/exp/totem/scratch/data/RP/test_06_25/xml/RP_all_new.xml'
  sourceFile[runNumber].params.runInfo = 'no info'
#  sourceFile[runNumber].params.trigger = "UandV1sectorOnly"

######################################################################
runs_14_10_2010_triggerParams = MyPSet(
  coincidenceChipConfig = cms.PSet(
    useControlRegisters = cms.bool(True),
    controlRegister1    = cms.uint32(19),
    controlRegister2    = cms.uint32(240),
    controlRegister3    = cms.uint32(128),
    useLogicWithWrongNP = cms.uint32(0)
  ),
  trigger = cms.string("UandV")
)

runToFileMap_14_10_2010 = {3462:[0,1],3463:[0],3464:[0,1],3465:[0,1], 3466:[0,1], 3467:[0,1], 3468:[0,1], 3469:[0,1], 3470:[0,1], 3578:[0,1],3580:[0,1],3581:[0,1]}
runToChosenTriggerPotIdMap =  {3462:[20,120],3463:[20,120],3464:[20,120],3465:[21,121], 3466:[24,124], 3467:[25,125], 3468:[20,121], 3469:[21,120], 3470:[24,125],3578:[21,120],3580:[25,124],3581:[24,125]}
for runNumber,fileNumbers in runToFileMap_14_10_2010.items(): 
  #runLabel = "run_"+str(run)
  #runs_september_2010_test_trigger.append(runLabel)
  sourceFile[runNumber] = SourceFile(runs_14_10_2010_triggerParams)
  sourceFile[runNumber].params.runName = str(runNumber)
  #for fileNumber in runToFileMap_14_10_2010[runNumber]:
  for fileNumber in fileNumbers:
    #print extraNum
    sourceFile[runNumber].params.sourceFileNames.append('/castor/cern.ch/totem/LHCRawData/2010/RP/run_'+str(runNumber)+'.00'+str(fileNumber)+'.vmea')
  sourceFile[runNumber].params.xmlFileName = '/afs/cern.ch/user/j/jprochaz/scratch/jprochaz/cmssw_users/user/jprochaz/mapping/RP_all_new.xml'
  #sourceFile[runNumber].params.xmlFileName = '/data2/totem/RP_IP5/xml/RP_all_new.xml'
  #sourceFile[runNumber].params.xmlFileName = '/afs/cern.ch/exp/totem/scratch/data/RP/test_06_25/xml/RP_all_new.xml'
  sourceFile[runNumber].params.runInfo = 'no info'
  sourceFile[runNumber].params.chosenPotsId = runToChosenTriggerPotIdMap[runNumber]

  #sourceFile[runLabel].params.trigger = "UandV1sectorOnly"
 #####################################################################

