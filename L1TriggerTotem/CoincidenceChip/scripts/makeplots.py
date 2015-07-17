import os, time, subprocess
import root2matplotlib as r2m
import ROOT
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import make_axes_locatable

import argparse 

import pylab
######################################
parser  = argparse.ArgumentParser(description=" a description")

parser.add_argument(
    "-outputdir",
#    dest    = 'outputDir',
#    metavar    = 'outputDir',
#    nargs = '*',
    default = "",
#    type    = int,
    action  = "store",
    help    = "directory where the output will be stored"
)

parser.add_argument(
    "-formats",
#    dest    = 'outputDir',
#    metavar    = 'outputDir',
#    nargs = '*',
     default = ["pdf"],
#    type    = int,
#    choices = ["pdf","png","eps"],
    action  = "store",
    help    = "formats of produced figures"
)

parser.add_argument(
    "-inputrootfile",
#    dest    = 'outputDir',
#    metavar    = 'outputDir',
#    nargs = '*',
    default = "./trigger_studies_def/run_1827/run_1827_trigger_studies_def.root",
#    type    = int,
    action  = "store",
    help    = "directory where the output will be stored"
)
#./trigger_studies_def/run_1827
args = parser.parse_args()

print args
###################################
#inputDir = "./outputDir2"
#inputDir = "./trigger_studies_def/run_1827"
outputDir = os.path.dirname(args.inputrootfile)


#print os.path.abspath()
#f =  ROOT.TFile("file:/afs/cern.ch/exp/totem/scratch/data/RP/test_06_25/conditions/lvdt.root","READ")
#f = ROOT.TFile("file:"+os.path.abspath(inputDir)+"/run_1827_trigger_studies_def.root","READ")
f = ROOT.TFile("file:"+args.inputrootfile,"READ")
#TIter next(GetListOfPrimitives());
#while ((TObject *obj = next()))
#obj->Draw(next.GetOption());
#potLabelList = ["rp_45_220_fr_bt", "rp_45_220_fr_tp", "rp_45_220_nr_tp"]
potLabelList = ["rp_45_220_fr_tp"]
            
#print "WERWE : ", time.strptime("25 Jun 10", "%d %b %y").tm_sec
#print "WERWE2: ", time.mktime(time.strptime("30 Nov 00", "%d %b %y"))
#  minTime = 1277491200 #// date -d "Fri Jun 25 20:40:00 2010" +%s
#  maxTime = 1277497237 #// date -d "Fri Jun 25 22:20:37 2010" +%s
ROOT.gROOT.Macro("$(HOME)/rootlogon.C")
ROOT.gStyle.SetOptTitle(1)

#assert(0)
can = {}
hist2D = {}
projYhist = {}
projXhist = {}
hist2D = {}
RPNPlanesDividedByTwo = 5
NSectors = 16
pyfig = []

counter = 0

def newHist(rootTH2D):
      global pyfig
      pyfig += [plt.figure(figsize=(20,10))]
      ax = plt.subplot(111)

      #ax.set_xticks(range(NSectors),minor=True)
      #plt.xticks(range(NSectors))
      #ax.set_aspect(1.)
      #ax = plt.axes()
      #hist = r2m.Hist2D(hist2D[prehist2D.GetName()])
      hist = r2m.Hist2D(rootTH2D)
      plt.suptitle(hist.title,fontsize="large")
      pot = ""
#      for s in hist.title:
#        if s != "_":
#          pot+=s
#        else:
#          break
#      print "POT=", pot
#      print " th2d name = ", rootTH2D.GetName()
      for rpName in potLabelList:
        if rootTH2D.GetName().startswith(rpName): pot = rpName
      print "POT2=", pot
    



    #  tx = pylab.text(binz, scont, cont,  **alignment1)
      #pylab.text(0.5, 0.5, 'matplotlib'+str(int(sum(hist.content))),
      #pylab.text(1.3, 1.3, 'Entries: 2776',
      print 
      pylab.text(1.15, 1.3, 'Entries: '+str(int(np.sum(hist.content))) ,
        horizontalalignment='center',
        verticalalignment='center',
        transform = ax.transAxes
      )
      hist.colz()
      print hist

#      ty = pylab.text(cont, binz+1, cont,rotation=-90 , **alignment2)
      
      #plt.tick_params(axis="x",which="both")
      
      projectionY = [0]*hist.nbinsy
      projectionX = [0]*hist.nbinsx

      for binx in range(hist.nbinsx):
        for biny in range(hist.nbinsy):
         # print binx, biny
          xycontent = hist.content[biny][binx]
          if xycontent == 0. or xycontent % (xycontent) == 0.: 
            xycontent = int(xycontent)
          t = pylab.text(binx,  biny+1,xycontent , ha = 'center', va = 'center')
          projectionY[biny] += xycontent 
          projectionX[binx] += xycontent 
      #projectionY += [ sum(hist.content[biny])]
      divider = make_axes_locatable(ax)
      axHistx = divider.append_axes("top", 2., pad=0.1, sharex=ax)
      plt.xticks(range(NSectors))
      axHisty = divider.append_axes("right", 2.2, pad=0.1, sharey=ax)
      axHistx.set_title("Sector projection")
      alignment0 = {'horizontalalignment':'center', 'verticalalignment':'top'}
      #axHisty.set_title("Plane projection title",rotation=-90,**alignment0)
      axHisty.set_title("Plane projection")

      #print plt.gca()
      #plt.figure()
      #ax2 = plt.gca()
      #ax2.barh(hist.yedges,projectionY)
      #print plt.gca()
      plt.sca(axHistx)
      axHistx.bar(hist.xedges[:-1], projectionX, width=1)
      
      #print number of entries for each plane/sector projection
      for binz in range(hist.nbinsx):
        cont = projectionX[binz] 
        if cont == 0. or cont % int(cont) == 0.:
          cont = int(cont)
        tx = pylab.text(binz, cont, cont,  ha='center', va='bottom')

      plt.sca(axHisty)
      axHisty.barh(hist.yedges[:-1], projectionY, height=1)

      for binz in range(hist.nbinsy):
        cont = projectionY[binz] 
        if cont == 0. or cont % int(cont) == 0.:
          cont = int(cont)
        ty = pylab.text(cont, binz+1, cont,rotation=-90 , va = 'center', ha='left')

      # rotate labels
      for label in axHisty.get_xticklabels():
        label.set_rotation(-90)

      # set limits
      ax.set_xlim(hist.xedges[0],hist.xedges[-1:][0])
      ax.set_ylim(hist.yedges[0],hist.yedges[-1:][0])
      axHistx.set_ylim(0,max(projectionX)*1.2)
      axHisty.set_xlim(0,max(projectionY)*1.2)

      # set labels
      ax.set_xlabel(hist.xlabel)
      ax.set_ylabel(hist.ylabel)
      
      ax.set_xticks(range(NSectors),minor=True)
 
      #axHisty.set_ylabel("xPlane projection title 3", x=2*axHisty.get_xlim()[1])
      #axHisty.set_ylabel("xPlane projection title 3")

      # Do not display ticks
      axHisty.yaxis.set_ticks_position('none')
      axHistx.xaxis.set_ticks_position('none') #this doesnt work, why??
      for tl in axHistx.get_xticklabels():
          tl.set_visible(False)

      plt.show()
     # """
      #for figformat in ["pdf","png","eps"]:
      #for figformat in ["pdf"]:
      for figformat in args.formats:
        myfig = outputDir+"/"+pot+'/'+hist.title+"."+figformat
        subprocess.Popen("mkdir --parents `dirname "+myfig+'`', shell=True)
        plt.savefig(myfig)
        print myfig
     # """
         # can[key.GetName()].SaveAs(outputDir+"/"+can[key.GetName()].GetName()+"."+figformat)
      
      #plt.savefig('colz')

      #plt.clf() # clear figure
      #ax = plt.axes(aspect='equal')
      #hist.contour()
      #plt.savefig('contour')


    #  build = False # False means do not draw default canvas...
    #  can = ROOT.TCanvas(build) 
    #  f.GetObject(key.GetName(),can)



#for potLabel in potLabelList:
for key in f.GetListOfKeys():
  #if key.GetName() == potLabel+"_plane0_occupancy_position_beamProfile":
  # print("class name :"+key.ReadObj().ClassName())
  #print("      name :"+key.GetName())
  #print("     title :"+key.GetTitle())
  # key.Print()
  #if key.ReadObj().ClassName() == "TH2D":
  if key.ReadObj().InheritsFrom("TH2D"):
#    for rpName in potLabelList:
#      if !rootTH2D.GetName().startswith(rpName): break
    #print key.ReadObj().InheritsFrom("TH2D")
    counter +=1
    #if counter > 1: break
    #build = False # False means do not draw default canvas...
    prehist2D = ROOT.TH2D() 
    f.GetObject(key.GetName(),prehist2D)
    #prehist2D.SetName("new"+prehist2D.GetName())
    newHist(prehist2D)
    #break
"""
      can[key.GetName()] = ROOT.TCanvas("can"+key.GetName(),key.GetTitle())
      can[key.GetName()].Divide(2,2)
      can[key.GetName()].cd(1)
      hist2D[prehist2D.GetName()] = prehist2D
      hist2D[prehist2D.GetName()].Draw("")
      hist2D[prehist2D.GetName()].SetDrawOption("COLZTEXT")
      #hist2D[prehist2D.GetName()].GetXaxis()->SetNdivisions(0*10000+0*100+fPotCollection[j].GetNSectors())
      #hist2D[prehist2D.GetName()].GetYaxis()->SetNdivisions(0*10000+0*100+fPotCollection[j].GetNPlanes())
      hist2D[prehist2D.GetName()].GetXaxis().SetNdivisions(NSectors)
      hist2D[prehist2D.GetName()].GetYaxis().SetNdivisions(RPNPlanesDividedByTwo)

      can[key.GetName()].cd(2)
      projYhist[prehist2D.GetName()] = hist2D[prehist2D.GetName()].ProjectionY("ProjectionY")
      projYhist[prehist2D.GetName()].SetTitle(hist2D[prehist2D.GetName()].GetTitle()) # this line should not be needed but because of a bug probably in ROOT we need it...
      projYhist[prehist2D.GetName()].Draw("")
      projYhist[prehist2D.GetName()].GetYaxis().SetRangeUser(0,projYhist[prehist2D.GetName()].GetMaximum()*1.1)
      #fHist2D[i].Write()

      can[key.GetName()].cd(3)
      projXhist[prehist2D.GetName()] = hist2D[prehist2D.GetName()].ProjectionX("ProjectionX")
      projXhist[prehist2D.GetName()].SetTitle(hist2D[prehist2D.GetName()].GetTitle())  # this line should not be needed but because of a bug probably in ROOT we need it...
      projXhist[prehist2D.GetName()].Draw("")

      print("hist2="+hist2D[prehist2D.GetName()].GetTitle())
      print("projY="+projYhist[prehist2D.GetName()].GetTitle())
      print("projX="+projXhist[prehist2D.GetName()].GetTitle())
      #hist2D.Write()

      #can[key.GetName()].Draw()
      for figformat in ["pdf","png"]:
        pass
         # can[key.GetName()].SaveAs(outputDir+"/"+can[key.GetName()].GetName()+"."+figformat)
      #plt.figure()
      #fig = plt.figure(1, figsize=(20,10))
"""
      
