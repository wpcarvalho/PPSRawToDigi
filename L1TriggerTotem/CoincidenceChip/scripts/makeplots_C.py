import os, time, subprocess
import ROOT

#inputDir = "./outputDir2"
inputDir = "./trigger_studies_def/run_1827"
outputDir = "./plots"

#print os.path.abspath()
#f =  ROOT.TFile("file:/afs/cern.ch/exp/totem/scratch/data/RP/test_06_25/conditions/lvdt.root","READ")
f = ROOT.TFile("file:"+os.path.abspath(inputDir)+"/run_1827_trigger_studies_def.root","READ")
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
print "after load"
ROOT.gStyle.SetOptTitle(1)
print "after set opt"

subprocess.Popen("mkdir --parents "+outputDir, shell=True)
#assert(0)
can = {}
hist2D = {}
projYhist = {}
projXhist = {}
hist2D = {}
RPNPlanesDividedByTwo = 5
NSectors = 16

counter = 0
for potLabel in potLabelList:
  for key in f.GetListOfKeys():
    #if key.GetName() == potLabel+"_plane0_occupancy_position_beamProfile":
    # print("class name :"+key.ReadObj().ClassName())
    #print("      name :"+key.GetName())
    #print("     title :"+key.GetTitle())
    # key.Print()
    if key.ReadObj().ClassName() == "TH2D":
      counter +=1
      if counter > 3: break
      #build = False # False means do not draw default canvas...
      can[key.GetName()] = ROOT.TCanvas("can"+key.GetName(),key.GetTitle())
      can[key.GetName()].Divide(2,2)
      can[key.GetName()].cd(1)
      prehist2D = ROOT.TH2D() 
      f.GetObject(key.GetName(),prehist2D)
      #key.ReadObj()
      prehist2D.SetName("new"+prehist2D.GetName())
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
      can[key.GetName()].Draw()
      for figformat in ["pdf","png"]:
         pass
         # can[key.GetName()].SaveAs(outputDir+"/"+can[key.GetName()].GetName()+"."+figformat)

    #  build = False # False means do not draw default canvas...
    #  can = ROOT.TCanvas(build) 
    #  f.GetObject(key.GetName(),can)


