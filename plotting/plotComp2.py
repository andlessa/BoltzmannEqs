#!/usr/bin/env python

import sys
sys.path.append('../pyCode')
from ROOT import *
import AuxPlot
import logging as logger
from AuxFuncs import getDataFrom

DoPrint = False

def getPlot():
    
    parDict,summaryDict,dataDict = getDataFrom('test.out')
    
    xprint = 'T^{-1} (1/GeV)'
    yprint = '#rho_{DM} (GeV^{2})'
    grDM = AuxPlot.getTGraph('test.out', xprint, yprint)
    grDM.SetLineColor(kAzure+2)
    grDM.SetLineWidth(3)
    yprint = '#rho_{Mediator} (GeV^{2})'
    grMed = AuxPlot.getTGraph('test.out', xprint, yprint)
    grMed.SetLineColor(kOrange)
    grMed.SetLineWidth(3)
    
    yprint = '#rho_{DM} (GeV^{2})'
    grDMlow = AuxPlot.getTGraph('test2.out', xprint, yprint)
    grDMlow.SetLineColor(kAzure+7)
    grDMlow.SetLineWidth(3)
    grDMlow.SetLineStyle(9)
    yprint = '#rho_{Mediator} (GeV^{2})'
    grMedlow = AuxPlot.getTGraph('test2.out', xprint, yprint)
    grMedlow.SetLineColor(kOrange+7)
    grMedlow.SetLineWidth(3)
    grMedlow.SetLineStyle(9)

    plane = TCanvas("c1", "c1",0,0,900,600)    
    AuxPlot.Default(plane,"TCanvas")
    
    
    plane.SetLeftMargin(0.13)
    plane.SetLogx()
    plane.SetLogy()
    plane.cd()
    base = TMultiGraph()
    base.Add(grDM,'L')
    base.Add(grMed,'L')
    base.Add(grDMlow,'L')
    base.Add(grMedlow,'L')    
    base.Draw('AL')
    base.GetYaxis().SetTitle('#rho')
    base.GetXaxis().SetTitle(xprint)
    base.GetYaxis().SetRangeUser(1e-41,1e25)
    AuxPlot.Default(base,"TGraph")
    
    pars = TPaveText(0.2,0.2,0.4,0.5,'brNDC')
    pars.SetFillColor(0)
    pars.SetTextFont(12)
    pars.AddText("#bf{Parameters} (GeV)")
    for key,val in parDict.items():
        pars.AddText(key+' = %3.4g'%val)    
    pars.Draw()
    
    #Legend
    nplots = int(base.GetListOfGraphs().GetSize())
    leg = TLegend(0.7,0.6,0.98,min(1.,0.6+0.1*nplots))
    AuxPlot.Default(leg,"Legend")
    nnext = TIter(base.GetListOfGraphs())
    gr = nnext.Next()
    while gr:
        leg.AddEntry(gr,gr.GetYaxis().GetTitle(),'L')
        gr = nnext.Next()
    
    leg.SetTextSize(0.0427337)
    leg.SetMargin(0.4)
    leg.Draw()
    
    legB = TLegend(0.68,0.2,0.95,0.35)
    legB.AddEntry(grDM,'T_{RH} = 10^{5} GeV','L')
    legB.AddEntry(grDMlow,'T_{RH} = 100 GeV','L')
    legB.SetTextSize(0.0427337)
    legB.SetMargin(0.3)
    AuxPlot.Default(legB,"Legend")
    legB.Draw()


    if DoPrint:
        plane.Print(dataFile+'.png')
    
    ans = raw_input("Hit any key to close\n")
    
    return
    
    
getPlot()
