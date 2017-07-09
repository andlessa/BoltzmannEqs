#!/usr/bin/env python

from ROOT import *
import AuxPlot
import logging as logger
from pyCode.AuxFuncs import getDataFrom

DoPrint = False

def getPlot():
    
    parDict,summaryDict,dataDict = getDataFrom(dataFile)
    
    xprint = 'R'
    yprint = '#rho_{DM} (GeV^{2})'
    grDM = AuxPlot.getTGraph('test.out', xprint, yprint)
    grDM.SetLineColor(kRed)
    grDM.SetLineWidth(3)
    yprint = '#rho_{Mediator} (GeV^{2})'
    grMed = AuxPlot.getTGraph('test.out', xprint, yprint)
    grMed.SetLineColor(kOrange)
    grMed.SetLineWidth(3)

    plane = TCanvas("c1", "c1",0,0,900,600)    
    AuxPlot.Default(plane,"TCanvas")
    
    
    plane.SetLeftMargin(0.13)
    plane.SetLogx()
    plane.SetLogy()
    plane.cd()
    base = TMultiGraph()
    base.Add(grDM,'L')
    base.Add(grMed,'L')
    base.Draw('AL')
    base.GetYaxis().SetTitle('#rho')
    base.GetXaxis().SetTitle(xprint)
    base.GetYaxis().SetRangeUser(1e-41,1e25)
    AuxPlot.Default(base,"TGraph")
    
    pars = TPaveText(0.2,0.2,0.4,0.5,'brNDC')
    pars.SetFillColor(0)
    pars.SetTextFont(12)
    pars.AddText("#bf{Parameters:}")
    for key,val in parDict.items():
        pars.AddText(key+' = %2.3g'%val)    
    pars.Draw()
    
    #Legend
    nplots = int(base.GetListOfGraphs().GetSize())
    leg = TLegend(0.77,0.6,0.98,min(1.,0.6+0.1*nplots))
    AuxPlot.Default(leg,"Legend")
    nnext = TIter(base.GetListOfGraphs())
    gr = nnext.Next()
    while gr:
        leg.AddEntry(gr,gr.GetYaxis().GetTitle(),'L')
        gr = nnext.Next()
    
    leg.SetTextSize(0.0427337)
    leg.SetMargin(0.2)
    leg.Draw()

    if DoPrint:
        plane.Print(dataFile+'.png')
    
    ans = raw_input("Hit any key to close\n")
    
    return
    
    
getPlot()