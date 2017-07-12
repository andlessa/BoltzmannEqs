#!/usr/bin/env python

import sys
sys.path.append('../pyCode')
from ROOT import *
import AuxPlot
import logging as logger
from AuxFuncs import getDataFrom

DoPrint = True

dataFile = 'test.out'

def getPlot():
    
    parDict,summaryDict,dataDict = getDataFrom(dataFile)
    
    ylabel = "Yield"
    
    xprint = ('1/[T (GeV)]','T^{-1} (1/GeV)')
    yprint = ('[n_{DM} (GeV^{3})]/[s (GeV^{3})]', 'Y_{DM}')
    grDM = AuxPlot.getTGraph('test.out', xprint, yprint)
    grDM.SetLineColor(kRed)
    grDM.SetLineWidth(3)
    yprint = ('[n_{Mediator} (GeV^{3})]/[s (GeV^{3})]', 'Y_{M}')
    grMed = AuxPlot.getTGraph('test.out', xprint, yprint)
    grMed.SetLineColor(kOrange)
    grMed.SetLineWidth(3)
    
    yprint = ('[s (GeV^{3})]','Y_{M}')
    grS = AuxPlot.getTGraph('test.out', xprint, yprint)
    grS.SetLineColor(kBlue)
    grS.SetLineWidth(3)    

    plane = TCanvas("c1", "c1",0,0,900,600)    
    AuxPlot.Default(plane,"TCanvas")
    
    
    plane.SetLeftMargin(0.13)
    plane.SetLogx()
    plane.SetLogy()
    plane.cd()
    base = TMultiGraph()
    base.Add(grDM,'L')
    base.Add(grMed,'L')
#     base.Add(grS,'L')
    base.Draw('AL')
    
    #Get all graphs:
    allGraphs = []
    nnext = TIter(base.GetListOfGraphs())
    gr = nnext.Next()
    while gr:
        allGraphs.append(gr)
        gr = nnext.Next()
    
    xmin = min([TMath.MinElement(gr.GetN(),gr.GetX()) for gr in allGraphs if TMath.MinElement(gr.GetN(),gr.GetX())])
    xmax = max([TMath.MaxElement(gr.GetN(),gr.GetX()) for gr in allGraphs if TMath.MaxElement(gr.GetN(),gr.GetX())])
    ymin = min([TMath.MinElement(gr.GetN(),gr.GetY()) for gr in allGraphs if TMath.MinElement(gr.GetN(),gr.GetY())])
    ymax = max([TMath.MaxElement(gr.GetN(),gr.GetY()) for gr in allGraphs if TMath.MaxElement(gr.GetN(),gr.GetY())])
    
    xmin *= 0.9
    ymin *= 0.9
    xmax *= 1.5
    ymax *= 1.5
    
    base.GetYaxis().SetTitle(ylabel)
    base.GetXaxis().SetTitle(xprint[1])
    base.GetYaxis().SetRangeUser(ymin,ymax)
    base.GetXaxis().SetLimits(xmin,xmax)
    AuxPlot.Default(base,"TGraph")
    
    pars = TPaveText(0.2,0.2,0.4,0.5,'brNDC')
    pars.SetFillColor(0)
    pars.SetTextFont(12)
    pars.AddText("#bf{Parameters:}")
    for key,val in parDict.items():
        pars.AddText(key+' = %2.3g'%val)    
    pars.Draw()
    
    #Legend     
    nplots = len(allGraphs)
    maxTitle = max([len(gr.GetYaxis().GetTitle()) for gr in allGraphs])
    leg = TLegend(min(0.85,0.99-0.0135*maxTitle),0.6,0.99,min(1.,0.6+0.1*nplots))
    AuxPlot.Default(leg,"Legend")
    for gr in allGraphs:
        leg.AddEntry(gr,gr.GetYaxis().GetTitle(),'L')
    
    leg.SetTextSize(0.0427337)
    leg.SetMargin(0.4)
    leg.Draw()

    if DoPrint:
        plane.Print(dataFile+'.png')
    
    ans = raw_input("Hit any key to close\n")
    
    return
    
    
getPlot()
