import ROOT
from odict import odict

#print runs_data_cff.extraRunToFileMap.keys()


#ROOT.gSystem.Load( "libReflex.so" )
#ROOT.gSystem.Load( "libCintex.so" )
#ROOT.Cintex.Enable()
#ROOT.gSystem.Load("libFWCoreFWLite.so")
#ROOT.AutoLibraryLoader.enable()

ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.AutoLibraryLoader.enable()
ROOT.gSystem.Load("libDataFormatsFWLite")
##ROOT.gSystem.Load("libDataFormatsPatCandidates")

ROOT.gSystem.Load("libL1TriggerTotemCoincidenceChip")
#ROOT.gSystem.AddIncludePath(" -I/home/jprochaz/workspace/cmssw/311/src/L1TriggerTotem/CoincidenceChip/interface ")
###ROOT.gSystem.AddIncludePath(" -I/home/jprochaz/workspace/cmssw/311/lib/slc4_ia32_gcc345/ ")
##ROOT.gSystem.AddIncludePath(" -I/home/jprochaz/workspace/cmssw/311/src/ ")
#ROOT.gSystem.Load("/data/usr/jiri/mnt/pctotem31/workspace/cmssw/311/lib/slc4_ia32_gcc345/libL1TriggerTotemCoincidenceChip")
##ROOT.gROOT.ProcessLine(".L ./src/L1TriggerTotemCoincidenceChipLinkDef.h++g")

import L1TriggerTotem.CoincidenceChip.runs_data_cff as runs_data_cff

_input_dir = ''

def returnFile(runNumber):
  #print "!@#$!@#$!@#$!@#$"
  #print [ int(runN.replace('run_','')) for runN in runs_data_cff.runs_13_06_2010]
  #print "!@#$!@#$!@#$!@#$"
  #if runNumber >= 2200 and runNumber <=2211:
  #  typeOut = 'allruns'
  #if runNumber in [ int(runN.replace('run_','')) for runN in runs_data_cff.runs_13_06_2010]:
  if runNumber in runs_data_cff.runs_13_06_2010:
    typeOut = 'runs_13_06_2010'
  #elif runNumber in [ int(runN.replace('run_','')) for runN in runs_data_cff.runs_30_06_2010]:
  elif runNumber in runs_data_cff.runs_30_06_2010:
    #typeOut = 'def2'
    typeOut = 'runs_30_06_2010'
  elif runNumber in runs_data_cff.runToFileMap_14_10_2010.keys():
    typeOut = 'runs_14_10_2010'
    #typeOut = 'def3'
  else:
    typeOut = 'def'
  #typeOut +='old'

  #print runNumber, '  ', typeOut
  return ROOT.TFile("file:"+_input_dir+'trigger_studies_'+typeOut+'/run_'+str(runNumber)+'/run_'+str(runNumber)+'_triggerStats.root',"READ")
  #return ROOT.TFile("file:"+'trigger_studies_'+typeOut+'/run_'+str(runNumber)+'/run_'+str(runNumber)+'_triggerStats.root',"READ")
  #return ROOT.TFile("file:"+'trigger_studies_'+typeOut+'/run_'+str(runNumber)+'/run_'+str(runNumber)+'_trigger_studies_def3.root',"READ")

def getInfoCollector(runNumber):
  f = returnFile(runNumber)
  for key in f.GetListOfKeys():
    #print "key        = ", key
    #print "name       = ", key.ReadObj().GetName()
    #print "class name = ", key.ReadObj().ClassName()
    #print "key name   = ", key.GetName()
    if key.ReadObj().ClassName() == "RPTriggerAnalyzerInfoCollector":
      infoCollector = key.ReadObj()
    if key.ReadObj().ClassName() == "PotCollection":
      potCollection = key.ReadObj() 
  #print ""

  return infoCollector

def getPotCollection(runNumber):
  f = returnFile(runNumber)
  for key in f.GetListOfKeys():
    #print "key        = ", key
    #print "name       = ", key.ReadObj().GetName()
    #print "class name = ", key.ReadObj().ClassName()
    #print "key name   = ", key.GetName()
    if key.ReadObj().ClassName() == "RPTriggerAnalyzerInfoCollector":
      infoCollector = key.ReadObj()
    if key.ReadObj().ClassName() == "PotCollection":
      potCollection = key.ReadObj() 
  #print ""
  return potCollection


def runInfo(runNumber):
  infoCollector =  getInfoCollector(runNumber)
  print "#############################################################"
  print "Stats global:" 
  print "NEvents         = ", infoCollector.fNEvents
  print "NEmptyEvents    = ", infoCollector.fNEmptyEvents 
  print "NonEmptyEvents  = ", infoCollector.fNEvents - infoCollector.fNEmptyEvents 
  
  for pot in getPotCollection(runNumber):
    print "#############################################################"
    print "RP ", pot.GetName(), "   (", pot.GetDecRPIdFull(),")"
    print ""
    print "NEmptyEvents = ", pot.fNEmptyEvents
    for trigStat in pot.fConditionalTriggerStat:
      print "condition                                      = ", trigStat.GetName()
      print "NEvents satisfaing condition                   = ", trigStat.fConditions                                     
      print "condition and RawSimuCCoutputTheSame           = ", trigStat.fConditionANDRawSimuCCoutputTheSame                                 
      print "condition and RawSimuCCoutputNotTheSame        = ", trigStat.fConditionANDRawSimuCCoutputNotTheSame           
      print "condition and RawSimuCCoutputNotTheSameEven    = ", trigStat.fConditionANDRawSimuCCoutputNotTheSameEven         
      print "condition and RawSimuCCoutputNotTheSameOdd     = ", trigStat.fConditionANDRawSimuCCoutputNotTheSameOdd   
      print "condition and RawSimuCCoutputNotTheSameEvenOdd = ", trigStat.fConditionANDRawSimuCCoutputNotTheSameEvenOdd     
      print "condition and TriggerRaw0Simu0                 = ", trigStat.fConditionANDTriggerRaw0Simu0 
      print "condition and TriggerRaw1Simu0                 = ", trigStat.fConditionANDTriggerRaw1Simu0                 
      print "condition and TriggerRaw0Simu1                 = ", trigStat.fConditionANDTriggerRaw0Simu1                 
      print "condition and TriggerRaw1Simu1                 = ", trigStat.fConditionANDTriggerRaw1Simu1                 
      print ""
    for potType in (ROOT.ERawOrSimu.kRaw,ROOT.ERawOrSimu.kSimu):  
      print "Pot ", pot.GetName(), ", ", pot.GetPot(potType).GetType(), " trigger"
      for trigStat in pot.GetPot(potType).fConditionalTriggerStat:
       # if trigStat.GetName() != "NoTrue":
        print "condition                = ", trigStat.GetName() 
        print "NEvents satisfaing cond. = ", trigStat.GetConditionNEvents()
        print "condition and Trigger    = ", trigStat.fConditionANDTrigger 
        print "condition and NoTrigger  = ", trigStat.fConditionANDNoTrigger
        print ""
      print ""
    print ""


#import subprocess
##import sys
##sys.path.append('./pymodules/common')
##sys.path.append('./pymodules')
#from L1TriggerTotem.CoincidenceChip.pymodules.common.odict.py import odict

#print subprocess.Popen("ls", stdout=subprocess.PIPE).communicate()[0].replace("_","\_")
#outputDir = "trigger_studies_allruns"
#outputDir = "trigger_studies_1827-1839"
#outputDir = "trigger_studies_def"
#command = "cat `find "+outputDir+" -name \*info -print` | grep 'FORTEX' | grep -v condition | sort"
"""
command = "find "+outputDir+" -name *log -print"
infoFiles = subprocess.Popen(command.split(), stdout=subprocess.PIPE ).communicate()[0].split()
#print infoFiles.replace("_","\_")
#"{0:.{precision}f}".format(float(a),precision=6)
print infoFiles
"""

def returnQuantities(runNumber):
  infoCollector =  getInfoCollector(runNumber)
  quantities =[]
#  f = returnFile(runNumber)
#  for key in f.GetListOfKeys():
#    if key.ReadObj().ClassName() == "PotCollection":
#      potCollection = key.ReadObj() 

  for pot in getPotCollection(runNumber):
    for trigStat in pot.fConditionalTriggerStat:
      quantities.append(odict())
      #quantities[-1]["runName"] = "testRunName"
      quantities[-1]["runName"] = infoCollector.runName.replace('run_','')
      quantities[-1]["potName"] = pot.GetName()
      quantities[-1]["condition"] = trigStat.GetName()
      quantities[-1]["nevents"] = trigStat.fConditions
      quantities[-1][""] = trigStat.fConditionANDRawSimuCCoutputTheSame                                 
      quantities[-1]["fConditionANDRawSimuCCoutputNotTheSame"] = trigStat.fConditionANDRawSimuCCoutputNotTheSame           
      quantities[-1]["fConditionANDRawSimuCCoutputNotTheSameEven"] = trigStat.fConditionANDRawSimuCCoutputNotTheSameEven         
      quantities[-1]["fConditionANDRawSimuCCoutputNotTheSameOdd"] = trigStat.fConditionANDRawSimuCCoutputNotTheSameOdd   
      quantities[-1]["fConditionANDRawSimuCCoutputNotTheSameEvenOdd"] = trigStat.fConditionANDRawSimuCCoutputNotTheSameEvenOdd     
      quantities[-1]["fConditionANDTriggerRaw0Simu0"] = trigStat.fConditionANDTriggerRaw0Simu0 
      quantities[-1]["fConditionANDTriggerRaw1Simu0"] = trigStat.fConditionANDTriggerRaw1Simu0                 
      quantities[-1]["fConditionANDTriggerRaw0Simu1"] = trigStat.fConditionANDTriggerRaw0Simu1                 
      quantities[-1]["fConditionANDTriggerRaw1Simu1"] = trigStat.fConditionANDTriggerRaw1Simu1 
      for potType in (ROOT.ERawOrSimu.kRaw,ROOT.ERawOrSimu.kSimu):  
        for trigStat in pot.GetPot(potType).fConditionalTriggerStat:
          if  quantities[-1]["condition"] == trigStat.GetName(): 
            #print "fConditionANDTrigger"+str(ROOT.PreTriggerStat.RawOrSimuString(potType))
            quantities[-1]["fConditionANDTrigger"+str(ROOT.PreTriggerStat.RawOrSimuString(potType))] = trigStat.fConditionANDTrigger
  assert(len(quantities)!=0)
  return quantities

#print quantities
"""
trigger_statistics_table = []
for quant in quantities:
  trigger_statistics_table.append(odict([
    ("Run"           , "%s" % (quant["runName"].replace("_",r"\_"))),
    ("Pot"           , "%s" % (quant["potName"].replace("_220","").replace("_",r"\_"))),
    ("condition"     , "%s" % (quant["condition"])),
    ("NEvents"       , "%u" % (quant["nevents"])),
    (r"$\frac{TrigTM_{raw}}{NEvents}$", "%u" % (quant["fConditionANDTriggerRaw"])),
    (r"$\frac{TrigTM_{sim}}{NEvents}$", "%u" % (quant["fConditionANDTriggerSimu"])),
  
    (r"$\frac{RSTM_{!=}}{NEvents}$", "%u" % (quant["fConditionANDRawSimuCCoutputNotTheSame"])),
    (r"$\frac{RSTM_{!=Even}}{NEvents}$", "%u" % (quant["fConditionANDRawSimuCCoutputNotTheSameEven"])),
    (r"$\frac{RSTM_{!=Odd}}{NEvents}$", "%u" % (quant["fConditionANDRawSimuCCoutputNotTheSameOdd"])),
    (r"$\frac{RSTM_{!=EvenOdd}}{NEvents}$", "%u" % (quant["fConditionANDRawSimuCCoutputNotTheSameEvenOdd"])),
  
    (r"$\frac{TRaw0Simu0}{NEvents}$", "%u" % (quant["fConditionANDTriggerRaw0Simu0"])),
    (r"$\frac{TRaw1Simu0}{NEvents}$", "%u" % (quant["fConditionANDTriggerRaw1Simu0"])),
    (r"$\frac{TRaw0Simu1}{NEvents}$", "%u" % (quant["fConditionANDTriggerRaw0Simu1"])),
    (r"$\frac{TRaw1Simu1}{NEvents}$", "%u" % (quant["fConditionANDTriggerRaw1Simu1"]))
  ])) 

"""    

def tableRaw(runNumber):
  trigger_statistics_table2 = []
  for quant in returnQuantities(runNumber): 
    #print type(quant)
    nevents = float(quant["nevents"])
    if nevents:
      forFormat="%.3g"
    else: 
      forFormat="%.3f"
    if nevents == 0:
       nevents = float('nan')
  
    #print type(nevents)
    #print type(float(quant["fConditionANDTriggerRaw"]))
    #print (float(quant["fConditionANDTriggerRaw"])/nevents)
    trigger_statistics_table2.append(odict([
      ("Run"           , "%s" % (quant["runName"].replace("_",r"\_"))),
      ("Pot"           , "%s" % (quant["potName"].replace("_220","").replace("rp_","").replace("_",r"\_"))),
      ("condition"     , "%s" % (quant["condition"])),
      ("NEvents"       , "%u" % (quant["nevents"])),
  
      (r"$\frac{TrigTM_{raw}}{NEvents}$", forFormat % (float(quant["fConditionANDTriggerRaw"])/nevents)),
      (r"$\frac{TrigTM_{sim}}{NEvents}$", forFormat % (float(quant["fConditionANDTriggerSimu"])/nevents)),
    
      (r"$\frac{RSTM_{!=}}{NEvents}$", forFormat % (float(quant["fConditionANDRawSimuCCoutputNotTheSame"])/nevents)),
      (r"$\frac{RSTM_{!=Even}}{NEvents}$", forFormat % (float(quant["fConditionANDRawSimuCCoutputNotTheSameEven"])/nevents)),
      (r"$\frac{RSTM_{!=Odd}}{NEvents}$", forFormat % (float(quant["fConditionANDRawSimuCCoutputNotTheSameOdd"])/nevents)),
      (r"$\frac{RSTM_{!=EvenOdd}}{NEvents}$", forFormat % (float(quant["fConditionANDRawSimuCCoutputNotTheSameEvenOdd"])/nevents)),
    
      (r"$\frac{TRaw0Simu0}{NEvents}$", forFormat % (float(quant["fConditionANDTriggerRaw0Simu0"])/nevents)),
      (r"$\frac{TRaw1Simu0}{NEvents}$", forFormat % (float(quant["fConditionANDTriggerRaw1Simu0"])/nevents)),
      (r"$\frac{TRaw0Simu1}{NEvents}$", forFormat % (float(quant["fConditionANDTriggerRaw0Simu1"])/nevents)),
      (r"$\frac{TRaw1Simu1}{NEvents}$", forFormat % (float(quant["fConditionANDTriggerRaw1Simu1"])/nevents))
    ])) 

  return trigger_statistics_table2



############################################################
############################################################
############################################################
"""
trigger_statistics_table = []
for quant in quantities:
  trigger_statistics_table.append(odict([
    ("Run"           , "{0:^}".format(quant["runName"].replace("_",r"\_"))),
    ("Pot"           , "{0:^}".format(quant["potName"].replace("_220","").replace("_",r"\_"))),
    ("condition"     , "{0:^}".format(quant["condition"])),
    ("NEvents"       , "{0:^}".format(quant["nevents"])),
    (r"$\frac{TrigTM_{raw}}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw"])),
    (r"$\frac{TrigTM_{sim}}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerSimu"])),
  
    (r"$\frac{RSTM_{!=}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSame"])),
    (r"$\frac{RSTM_{!=Even}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSameEven"])),
    (r"$\frac{RSTM_{!=Odd}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSameOdd"])),
    (r"$\frac{RSTM_{!=EvenOdd}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSameEvenOdd"])),
  
    (r"$\frac{TRaw0Simu0}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw0Simu0"])),
    (r"$\frac{TRaw1Simu0}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw1Simu0"])),
    (r"$\frac{TRaw0Simu1}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw0Simu1"])),
    (r"$\frac{TRaw1Simu1}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw1Simu1"]))
  ])) 

trigger_statistics_table2 = []
for quant in quantities:
  nevents = float(quant["nevents"])
  if nevents:
    forFormat="{0:^.3g}"
  else: 
    forFormat="{0:^.3}"
  trigger_statistics_table2.append(odict([
    ("Run"           , "{0:^}".format(quant["runName"].replace("_",r"\_"))),
    ("Pot"           , "{0:^}".format(quant["potName"].replace("_220","").replace("rp_","").replace("_",r"\_"))),
    ("condition"     , "{0:^}".format(quant["condition"])),
    ("NEvents"       , "{0:^}".format(quant["nevents"])),

    (r"$\frac{TrigTM_{raw}}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerRaw"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{TrigTM_{sim}}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerSimu"])/nevents if nevents!=0 else 'nan')),
  
    (r"$\frac{RSTM_{!=}}{NEvents}$", forFormat.format(float(quant["fConditionANDRawSimuCCoutputNotTheSame"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{RSTM_{!=Even}}{NEvents}$", forFormat.format(float(quant["fConditionANDRawSimuCCoutputNotTheSameEven"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{RSTM_{!=Odd}}{NEvents}$", forFormat.format(float(quant["fConditionANDRawSimuCCoutputNotTheSameOdd"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{RSTM_{!=EvenOdd}}{NEvents}$", forFormat.format(float(quant["fConditionANDRawSimuCCoutputNotTheSameEvenOdd"])/nevents if nevents!=0 else 'nan')),
  
    (r"$\frac{TRaw0Simu0}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerRaw0Simu0"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{TRaw1Simu0}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerRaw1Simu0"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{TRaw0Simu1}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerRaw0Simu1"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{TRaw1Simu1}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerRaw1Simu1"])/nevents if nevents!=0 else 'nan'))
  ])) 

#table_print(trigger_statistics_table[0], mode="latex", onlyTitle=True)
table_print(trigger_statistics_table, mode="latex", noTitle=True)
table_print(trigger_statistics_table2, mode="latex", noTitle=True)
"""

"""
trigger_statistics_table = []
for quant in quantities:
  trigger_statistics_table.append(odict([
    ("Run"           , "{0:^}".format(quant["runName"].replace("_",r"\_"))),
    ("Pot"           , "{0:^}".format(quant["potName"].replace("_220","").replace("_",r"\_"))),
    ("condition"     , "{0:^}".format(quant["condition"])),
    ("NEvents"       , "{0:^}".format(quant["nevents"])),
    (r"$\frac{TrigTM_{raw}}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw"])),
    (r"$\frac{TrigTM_{sim}}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerSimu"])),
  
    (r"$\frac{RSTM_{!=}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSame"])),
    (r"$\frac{RSTM_{!=Even}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSameEven"])),
    (r"$\frac{RSTM_{!=Odd}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSameOdd"])),
    (r"$\frac{RSTM_{!=EvenOdd}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSameEvenOdd"])),
  
    (r"$\frac{TRaw0Simu0}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw0Simu0"])),
    (r"$\frac{TRaw1Simu0}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw1Simu0"])),
    (r"$\frac{TRaw0Simu1}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw0Simu1"])),
    (r"$\frac{TRaw1Simu1}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw1Simu1"]))
  ])) 
"""
"""
import subprocess
import sys
sys.path.append('./pymodules/common')
sys.path.append('./pymodules')
#from L1TriggerTotem.CoincidenceChip.pymodules.common.odict.py import odict
from odict import odict
from tables import *

#print subprocess.Popen("ls", stdout=subprocess.PIPE).communicate()[0].replace("_","\_")
outputDir = "trigger_studies_allruns"
#outputDir = "trigger_studies_def"
#command = "cat `find "+outputDir+" -name \*info -print` | grep 'FORTEX' | grep -v condition | sort"
command = "find "+outputDir+" -name *log -print"
infoFiles = subprocess.Popen(command.split(), stdout=subprocess.PIPE ).communicate()[0].split()
#print infoFiles.replace("_","\_")
#"{0:.{precision}f}".format(float(a),precision=6)
print infoFiles
quantities =[]
for infoFile in infoFiles:
  f = open(infoFile)
  for line in f:
#    print line
    forTex = "FORTEX"
    if forTex in line:
       if "BEGIN_LABEL" in line: 
         quantities.append(odict())
         continue
       #print line
       quantity, value = line.replace(forTex,"").split()
       quantities[-1][quantity] = value

trigger_statistics_table = []
for quant in quantities:
  trigger_statistics_table.append(odict([
    ("Run"           , "{0:^}".format(quant["runName"].replace("_",r"\_"))),
    ("Pot"           , "{0:^}".format(quant["potName"].replace("_220","").replace("_",r"\_"))),
    ("condition"     , "{0:^}".format(quant["condition"])),
    ("NEvents"       , "{0:^}".format(quant["nevents"])),
    (r"$\frac{TrigTM_{raw}}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw"])),
    (r"$\frac{TrigTM_{sim}}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerSimu"])),
  
    (r"$\frac{RSTM_{!=}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSame"])),
    (r"$\frac{RSTM_{!=Even}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSameEven"])),
    (r"$\frac{RSTM_{!=Odd}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSameOdd"])),
    (r"$\frac{RSTM_{!=EvenOdd}}{NEvents}$", "{0:^}".format(quant["fConditionANDRawSimuCCoutputNotTheSameEvenOdd"])),
  
    (r"$\frac{TRaw0Simu0}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw0Simu0"])),
    (r"$\frac{TRaw1Simu0}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw1Simu0"])),
    (r"$\frac{TRaw0Simu1}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw0Simu1"])),
    (r"$\frac{TRaw1Simu1}{NEvents}$", "{0:^}".format(quant["fConditionANDTriggerRaw1Simu1"]))
  ])) 

trigger_statistics_table2 = []
for quant in quantities:
  nevents = float(quant["nevents"])
  if nevents:
    forFormat="{0:^.3g}"
  else: 
    forFormat="{0:^.3}"
  trigger_statistics_table2.append(odict([
    ("Run"           , "{0:^}".format(quant["runName"].replace("_",r"\_"))),
    ("Pot"           , "{0:^}".format(quant["potName"].replace("_220","").replace("_",r"\_"))),
    ("condition"     , "{0:^}".format(quant["condition"])),
    ("NEvents"       , "{0:^}".format(quant["nevents"])),

    (r"$\frac{TrigTM_{raw}}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerRaw"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{TrigTM_{sim}}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerSimu"])/nevents if nevents!=0 else 'nan')),
  
    (r"$\frac{RSTM_{!=}}{NEvents}$", forFormat.format(float(quant["fConditionANDRawSimuCCoutputNotTheSame"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{RSTM_{!=Even}}{NEvents}$", forFormat.format(float(quant["fConditionANDRawSimuCCoutputNotTheSameEven"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{RSTM_{!=Odd}}{NEvents}$", forFormat.format(float(quant["fConditionANDRawSimuCCoutputNotTheSameOdd"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{RSTM_{!=EvenOdd}}{NEvents}$", forFormat.format(float(quant["fConditionANDRawSimuCCoutputNotTheSameEvenOdd"])/nevents if nevents!=0 else 'nan')),
  
    (r"$\frac{TRaw0Simu0}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerRaw0Simu0"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{TRaw1Simu0}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerRaw1Simu0"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{TRaw0Simu1}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerRaw0Simu1"])/nevents if nevents!=0 else 'nan')),
    (r"$\frac{TRaw1Simu1}{NEvents}$", forFormat.format(float(quant["fConditionANDTriggerRaw1Simu1"])/nevents if nevents!=0 else 'nan'))
  ])) 

#table_print(trigger_statistics_table[0], mode="latex", onlyTitle=True)
table_print(trigger_statistics_table, mode="latex", noTitle=True)
table_print(trigger_statistics_table2, mode="latex", noTitle=True)
"""
