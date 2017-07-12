#!/usr/bin/env python
import os,sys,copy
sys.path.append('../')
from ROOT import TTree,TColor,TCanvas,TF1,TGraph,Double,TFile,gDirectory,TNamed
import logging as logger
from pyCode import AuxFuncs
from numpy import pi
  

def getTGraph(dataFile,xprint,yprint):
    
    dataDict = AuxFuncs.getDataFrom(dataFile)[-1]

    xlabel = xprint
    ylabel = yprint    
    if xprint == 'T^{-1} (1/GeV)':
        xlabel = 'T (GeV)'

    if not xlabel in dataDict:
        logger.error('Variable %s not found in dataDict. \n Variables available:\n %s' %(xlabel,dataDict.keys()))
        return None
    if not ylabel in dataDict:
        logger.error('Variable %s not found in dataDict. \n Variables available:\n %s' %(ylabel,dataDict.keys()))
        return None
    
    xpts, ypts = dataDict[xlabel],dataDict[ylabel]
    gr = TGraph()
    for i,x in enumerate(xpts):
#         y = ypts[i]/(2.*pi*AuxFuncs.gSTARS(x)*x**3/45.)  ##Yield
        y = ypts[i]
        if xprint == 'T^{-1} (1/GeV)':
            x = 1./x
        gr.SetPoint(gr.GetN(),x,y)

    gr.GetXaxis().SetTitle(xlabel)
    gr.GetYaxis().SetTitle(ylabel)
    Default(gr,"TGraph")
    
    return gr
    
def Default(obj,Type):
  
    if Type == "TCanvas":
        obj.SetLeftMargin(0.1097891)
        obj.SetRightMargin(0.02700422)
        obj.SetTopMargin(0.02796053)
        obj.SetBottomMargin(0.14796053)
        obj.SetFillColor(0)
        obj.SetBorderSize(0)
        obj.SetFrameBorderMode(0)
    elif "TGraph" in Type or "TH" in Type:
        obj.GetYaxis().SetTitleFont(132)
        obj.GetYaxis().SetTitleSize(0.065)
        obj.GetYaxis().CenterTitle(True)
        obj.GetYaxis().SetTitleOffset(0.9)
        obj.GetXaxis().SetTitleFont(52)
        obj.GetXaxis().SetTitleSize(0.065)
        obj.GetXaxis().CenterTitle(True)
        obj.GetXaxis().SetTitleOffset(1.0)
        obj.GetYaxis().SetLabelFont(132)
        obj.GetXaxis().SetLabelFont(132)
        obj.GetYaxis().SetLabelSize(0.05)
        obj.GetXaxis().SetLabelSize(0.05)
    if "TGraph2D" in Type or "TH2" in Type:
        obj.GetZaxis().SetTitleFont(132)
        obj.GetZaxis().SetTitleSize(0.06)
        obj.GetZaxis().CenterTitle(True)
        obj.GetZaxis().SetTitleOffset(0.7)
        obj.GetZaxis().SetLabelFont(132)
        obj.GetZaxis().SetLabelSize(0.05)
    elif "Leg" in Type:
        obj.SetBorderSize(1)
        obj.SetTextFont(132)
        obj.SetTextSize(0.05)
        obj.SetLineColor(1)
        obj.SetLineStyle(1)
        obj.SetLineWidth(1)
        obj.SetFillColor(0)
        obj.SetFillStyle(1001)


    
    
def set_palette(gStyle,name="none", ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""
    
    from array import array

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)
