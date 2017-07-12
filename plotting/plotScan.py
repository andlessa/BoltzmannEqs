#!/usr/bin/env python

import sys,os
sys.path.append('../pyCode')
from ROOT import *
import AuxPlot
import logging as logger
from AuxFuncs import getDataFrom
from array import array

DoPrint = True


def getPlot():
    
    allFiles = []
    for root, dirs, files in os.walk("../note/yscan"):
        for f in files:
            allFiles.append(os.path.abspath(os.path.join(root,f)))
    allFiles = sorted(allFiles)

    Tdecay = []
    Tdecouple = []
    y = []
    omegaNum = []
    omegaAna = []
    for f in allFiles:    
        parDict,summaryDict,dataDict = getDataFrom(f)
        y.append(parDict['yCoupling'])
        Tdecay.append(summaryDict['Mediator']['T(decay)~'])
        Tdecouple.append(summaryDict['Mediator']['T(decouple)~'])
        omegaNum.append(summaryDict['DM']['Omega h^2 (@TF)'])
        #Get analytical result:
        ff = open(f,'r')
        lines = [l.replace('\n','').strip() for l in ff.readlines()]
        lines = [l for l in lines if l] #Remove empty lines
        omegaAna.append(eval([l for l in lines if '(analytic)' in l][0].split('=')[-1]))

    #Sort according to y:
    omegaAna = sorted(omegaAna, key=lambda x_el: y[omegaAna.index(x_el)])
    omegaNum = sorted(omegaNum, key=lambda x_el: y[omegaNum.index(x_el)])
    Tdecay = sorted(Tdecay, key=lambda x_el: y[Tdecay.index(x_el)])
    Tdecouple = sorted(Tdecouple, key=lambda x_el: y[Tdecouple.index(x_el)])
    y = sorted(y)
    
    
#     grAna = TGraph(len(y),array('d',y),array('d',omegaAna))
#     grAna.SetTitle("Analytical Solution")
#     grNum = TGraph(len(y),array('d',y),array('d',omegaNum))
#     grNum.SetTitle("Numerical Solution")
#     AuxPlot.Default(grAna,'TGraph')
#     grAna.SetLineColor(kOrange)
#     grAna.SetLineWidth(3)
#     AuxPlot.Default(grNum,'TGraph')
#     grNum.SetLineColor(kRed)
#     grNum.SetLineWidth(3)
#     grNum.SetLineStyle(9)

    grAna = TGraph(len(y),array('d',y),array('d',Tdecay))
    grAna.SetTitle("T_{decay}")
    grNum = TGraph(len(y),array('d',y),array('d',Tdecouple))
    grNum.SetTitle("T_{freeze-out}")
    AuxPlot.Default(grAna,'TGraph')
    grAna.SetLineColor(kOrange)
    grAna.SetLineWidth(3)
    AuxPlot.Default(grNum,'TGraph')
    grNum.SetLineColor(kRed)
    grNum.SetLineWidth(3)




    plane = TCanvas("c1", "c1",0,0,900,600)    
    AuxPlot.Default(plane,"TCanvas")
    
    
    plane.SetLeftMargin(0.13)
    plane.SetRightMargin(0.13)
    plane.SetLogx()
    plane.SetLogy()
    plane.cd()
    base = TMultiGraph()
    base.Add(grAna,'L')
    base.Add(grNum,'L')
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
#     base.GetYaxis().SetTitle('#Omega_{DM} h^{2}')
    base.GetYaxis().SetTitle('T (GeV)')
    base.GetXaxis().SetTitle('y')
    base.GetYaxis().SetRangeUser(ymin,ymax)
    base.GetXaxis().SetLimits(xmin,xmax)
    AuxPlot.Default(base,"TGraph")
    
    nplots = len(allGraphs)
    maxTitle = max([len(gr.GetTitle()) for gr in allGraphs])
    leg = TLegend(min(0.85,0.99-0.0135*maxTitle),0.6,0.99,min(1.,0.6+0.1*nplots))
    AuxPlot.Default(leg,"Legend")
    for gr in allGraphs:
        leg.AddEntry(gr,gr.GetTitle(),'L')
    
    leg.SetTextSize(0.0427337)
    leg.SetMargin(0.4)
    leg.Draw()

    if DoPrint:
        plane.Print('scan.png')
    
    
    ans = raw_input("Hit any key to close\n")    
    
    return
    
    
getPlot()
