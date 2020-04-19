import numpy as np
import yaml
import argparse
import pandas as pd
from ROOT import TFile, TGraphAsymmErrors, TCanvas, TLatex, TLegend, TH1F, TLine
from StyleFormatter import SetGlobalStyle, SetObjectStyle

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('infileFONLL', metavar='text', default='FONLL.txt', help='input FONLL txt file name')
args = parser.parse_args()

SetGlobalStyle(padleftmargin=0.16, padbottommargin=0.14)

dfModel = pd.read_csv(args.infileFONLL, sep=' ')

gFONLL = TGraphAsymmErrors(len(dfModel))
for iPt, (ptMin, ptMax, crossSec, crossSecMin, crossSecMax) in enumerate(zip(\
    dfModel['ptmin'].values, dfModel['ptmax'].values, dfModel['central'].values, dfModel['min'].values, dfModel['max'].values)):

    print((ptMin+ptMax)/2, crossSec)
    gFONLL.SetPoint(iPt, (ptMin+ptMax)/2, crossSec / (ptMax-ptMin))
    gFONLL.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, (crossSec-crossSecMin) / (ptMax-ptMin), (crossSecMax-crossSec) / (ptMax-ptMin))

SetObjectStyle(gFONLL, color=864, fillalpha=0.2, markerstyle=20)

cFONLL = TCanvas('cFONLL', '', 500, 500)
cFONLL.DrawFrame(0., 1.e4, 24., 1.e8, ';d#it{p}_{T} (GeV/#it{c});(d^{2}#sigma/d#it{p}_{T}d#it{y}) (pb^{-1} GeV^{-1} #it{c})')
cFONLL.SetLogy()
gFONLL.Draw('2')
gFONLL.Draw('PX')

input()