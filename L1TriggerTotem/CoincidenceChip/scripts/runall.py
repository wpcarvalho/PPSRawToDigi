import subprocess
import sys,os
#sys.path.append("/data/usr/jiri/mnt/pctotem31/workspace/cmssw/311/src/L1TriggerTotem/CoincidenceChip/python/")

import L1TriggerTotem.CoincidenceChip.runs_data_cff as runs_data_cff
#from testbeam_data_cff import latRuns,rpMapping
#from L1TriggerTotem.CoincidenceChip.TriggerAnalyzerParser_cfi import parser
#from L1TriggerTotem.CoincidenceChip.test.TestBeamDataAnalysis_potX_cfg import parser

cfgFile = os.environ['CMSSW_BASE']+'/src/L1TriggerTotem/CoincidenceChip/test/triggerAnalyzer_cfg.py'

def create_runs_13_06_2010():
  for runNumber in runs_data_cff.runs_13_06_2010:
    print 'runtime for run ', runNumber,'is (may take same time...):'
    outputDir = "trigger_studies_"+"runs_13_06_2010"
    runOutputDir = outputDir+"/run_"+str(runNumber)
    runLogFile = runOutputDir+"/run_"+str(runNumber)+".log"
    chosenPot = str(runs_data_cff.runs_13_06_2010_Params_triggerOnPot[runNumber])
    subprocess.call("mkdir -p "+runOutputDir, shell=True)
    test = subprocess.call("time cmsRun "+cfgFile+" +nevents=-1 +run="+str(runNumber)+" +verbosity=0 +chosenPotsId "+chosenPot+" +outputdir="+outputDir+" &> "+runLogFile+" &", shell=True)
    #if not test: break

def create_runs_30_06_2010():
  for runNumber in runs_data_cff.runs_30_06_2010:
    print 'runtime for run ', runNumber,'is (may take same time...):'
    # if int(run.replace('run_','')) != 2211: continue
    outputDir = "trigger_studies_"+"runs_30_06_2010"
    runOutputDir = outputDir+"/run_"+str(runNumber)
    runLogFile = runOutputDir+"/run_"+str(runNumber)+".log"
    chosenPot = str(runs_data_cff.runs_30_06_2010_Params_triggerOnPot[runNumber])
    subprocess.call("mkdir -p "+runOutputDir, shell=True)
    test = subprocess.call("time cmsRun "+cfgFile+" +nevents=-1 +run="+str(runNumber)+" +verbosity=0 +chosenPotsId "+chosenPot+" +outputdir="+outputDir+" &> "+runLogFile+" &", shell=True)

def create_runs_14_10_2010():   
  for runNumber in runs_data_cff.runToFileMap_14_10_2010.keys():
    print 'runtime for run ', runNumber,'is (may take same time...):'
    #if runNumber != 3462: continue
    #if int(run.replace('run_','')) <= 3470: continue
    outputDir = "trigger_studies_"+"runs_14_10_2010"
    #subprocess.Popen([r"rm","-r",outputDir]).wait()
    runOutputDir = outputDir+"/run_"+str(runNumber)
    runLogFile = runOutputDir+"/run_"+str(runNumber)+".log"
    #chosenPot = str(runs_data_cff.run2200Params_triggerOnPot[run])
    subprocess.call("mkdir -p "+runOutputDir, shell=True)
    test = subprocess.call("time cmsRun "+cfgFile+" +nevents=-1 +run="+str(runNumber)+" +verbosity=0 +outputdir="+outputDir+" &> "+runLogFile+" &", shell=True)

if 1 or subprocess.Popen([r"scram","b"]).wait()==0:
  create_runs_13_06_2010()
  create_runs_30_06_2010()
  #create_runs_14_10_2010()

##################################################
  """
  for run in ["run_1941", "run_1942"]:
    chosenPot = ""#str(runs_data_cff.run1827Params_triggerOnPot[run])
    subprocess.call("mkdir -p "+outputDir+"/"+run, shell=True)
    test = subprocess.call("time cmsRun test/triggerAnalyzer_cfg.py +nevents=-1 +run="+run+" +verbosity=2 +chosenPotsId "+chosenPot+" +outputdir="+outputDir+" | tee "+runLogFile, shell=True)
  #produceTable()  
  """

"""
  for run in runs_data_cff.runs_testBeamH8:
    if run== "run_11_A": continue
    if run== "run_12_A": continue
    if run== "run_12_B": continue
#    if run== "run_12_B": continue
    chosenPot = "" #str(run1827Params_triggerOnPot[run])
    subprocess.call("mkdir -p "+outputDir+"/"+run, shell=True)
    test = subprocess.call("time cmsRun test/triggerAnalyzer_cfg.py +nevents=-1 +run="+run+" +verbosity=2 +reco +chosenPotsId "+chosenPot+" +outputdir="+outputDir+" 2>&1 | tee "+outputDir+"/"+run+"/"+run+".log", shell=True)
  #  break
"""

###################################################################

"""
import Configuration.TotemCommon.argparse as argparse
#import argparse
parser = argparse.ArgumentParser(description=" a description of runall3")

parser.add_argument(
    "-table",
#    dest    = 'outputDir',
#    metavar    = 'outputDir',
#    nargs = '*',
    default = False,
#    type    = int,
    action  = "store_true",
    help    = "produce a table?"
)

parser.add_argument(
    "-plots",
#    dest    = 'outputDir',
#    metavar    = 'outputDir',
#    nargs = '*',
    default = False,
#    type    = int,
    action  = "store_true",
    help    = "produce plots?"
)

args = parser.parse_args()

# print outputDir
outputDir = "def3"
#print args.nevents
outputDir = "trigger_studies_"+outputDir

def producePlots():
  global outputDir
  subprocess.call("find "+outputDir+" -name \*.root  | xargs -n1 -t python scripts/makeplots.py -inputrootfile $@", shell=True)

def produceTable():
  global outputDir
  #subprocess.call("cat `find "+outputDir+" -name \*info -print` | grep '&*&*&' | grep -v condition | sort > "+outputDir+r"/table.tex; sed 's/_/\\_/g' "+outputDir+"/table.tex -i; pdflatex -shell-escape '\\newcommand{\mytable}{"+outputDir+"/table.tex} \input{scripts/trigger_table.tex}'; mv trigger_table.pdf "+outputDir, shell=True)
  subprocess.call("pdflatex -shell-escape '\\newcommand{\mytable}{"+outputDir+"} \input{scripts/trigger_table.tex}'; mv trigger_table.pdf "+outputDir, shell=True)

if args.table or args.plots:
  if args.table:
    produceTable()
  if args.plots:
    producePlots()
  exit(0)

#subprocess.Popen([r"rm","-r",outputDir]).wait()

#def _runOutputDir(runNumber):
#  return 
"""


##################################################
## very old...
"""
import subprocess
import sys
from L1TriggerTotem.CoincidenceChip.testbeam_data_cff import latRuns,rpMapping
#from L1TriggerTotem.CoincidenceChip.TriggerAnalyzerParser_cfi import parser
#from L1TriggerTotem.CoincidenceChip.test.TestBeamDataAnalysis_potX_cfg import parser

#args = parser.parse_args()

runs = latRuns
scenario = sys.argv[1]
print "scenario = ", sys.argv[1]
if scenario == "0":
  outputDir = "onepot"
else:
  outputDir = "allpots"

print outputDir
#print args.nevents
outputDir = "trigger_studies_"+outputDir
subprocess.Popen([r"rm","-r",outputDir]).wait()
if subprocess.Popen([r"sc","b"]).wait()==0:
  #print subprocess.Popen(["pwd"], stdout=subprocess.PIPE).communicate()[0]
  #subprocess.Popen(["cmsRun", "test/TestBeamDataAnalysis_potX_cfg.py", "+nevents=-1", "+run=run_693", "+verbosity=2", "+chosenPotsId", "&>", "test.test"])
  for run,latpot in runs.iteritems():
    lat,pot = latpot
    run = "run_"+str(run)
    if scenario=="0":
      xpots = str(rpMapping[pot])
    else: 
      xpots = " "
    subprocess.call("mkdir -p "+outputDir+"/"+run, shell=True)
    #subprocess.call("cmsRun test/TestBeamDataAnalysis_potX_cfg.py +nevents=-1 +run="+run+" +verbosity=2 +chosenPotsId "+str(rpMapping[pot])+" | tee trigger_studies/"+run+"/"+run+".info", shell=True)
    test = subprocess.call("cmsRun test/TestBeamDataAnalysis_potX_cfg.py +nevents=-1 +run="+run+" +verbosity=2 +chosenPotsId "+xpots+" +outputdir="+outputDir+" | tee "+outputDir+"/"+run+"/"+run+".info", shell=True)
    break
    #if not test: break
  subprocess.call("cat `find "+outputDir+" -name \*info -print` | grep '&*&*&' | grep -v condition > "+outputDir+r"/table.tex; sed 's/_/\\_/g' "+outputDir+"/table.tex -i; pdflatex '\\newcommand{\mytable}{"+outputDir+"/table.tex} \input{scripts/trigger_table.tex}'; mv trigger_table.pdf "+outputDir, shell=True)
"""
