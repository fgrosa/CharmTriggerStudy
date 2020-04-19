'''
Script for the computation of predictions for D+, D*+, Ds, Lc
'''

import numpy as np
import pandas as pd
from ROOT import TFile, TCanvas, TGraphAsymmErrors, TF1
from ROOT import gStyle, kBlack, kAzure


def SetGraphStyle(graph, color, marker, linewidth=2):
    '''
    Helper method to set graph style
    '''
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(1)
    graph.SetLineWidth(linewidth)


gStyle.SetPadTopMargin(0.035)
gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadRightMargin(0.035)
gStyle.SetTitleSize(0.050, 'xy')
gStyle.SetLabelSize(0.045, 'xy')
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

dfD0FONLL = pd.read_csv('models/FONLL-D0-pp14TeV-y05-measbins.txt', sep=' ')
ratiosDfile = TFile.Open('HEPData-ins1716440-v1-root.root')

# D-meson ratios
gDplusOverD0 = ratiosDfile.Get('Table 6/Graph1D_y1')
gDstarOverD0 = ratiosDfile.Get('Table 7/Graph1D_y1')
gDsOverD0 = ratiosDfile.Get('Table 8/Graph1D_y1')
SetGraphStyle(gDplusOverD0, kBlack, 20)
SetGraphStyle(gDstarOverD0, kBlack, 20)
SetGraphStyle(gDsOverD0, kBlack, 20)
fDplusOverD0 = TF1('fDplusOverD0', 'pol0', 0, 50)
fDstarOverD0 = TF1('fDstarOverD0', 'pol0', 0, 50)
fDsOverD0 = TF1('fDsOverD0', 'pol0', 0, 50)
SetGraphStyle(fDplusOverD0, kAzure+4, -1, 3)
SetGraphStyle(fDstarOverD0, kAzure+4, -1, 3)
SetGraphStyle(fDsOverD0, kAzure+4, -1, 3)

# Lc/D0 ratio
dfLcOverD0 = pd.read_csv('Lc_to_D0_ratio_pp5TeV_y05.txt', sep=' ')
dfLcOverD0['totunclow'] = dfLcOverD0.apply(
    lambda x: np.sqrt((x['central']-x['minstat'])**2+(x['central']-x['minsyst'])**2), axis=1)
dfLcOverD0['totunchigh'] = dfLcOverD0.apply(
    lambda x: np.sqrt((x['maxstat']-x['central'])**2+(x['maxsyst']-x['central'])**2), axis=1)
gLcOverD0 = TGraphAsymmErrors(len(dfLcOverD0))
for iPt, (ptMin, ptMax, cent, uncLow, uncHigh) in enumerate(zip(dfLcOverD0['ptmin'], dfLcOverD0['ptmax'], \
    dfLcOverD0['central'], dfLcOverD0['totunclow'], dfLcOverD0['totunchigh'])):
    gLcOverD0.SetPoint(iPt, (ptMax+ptMin)/2, cent)
    gLcOverD0.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, uncLow, uncHigh)
SetGraphStyle(gLcOverD0, 1, 20)
fLcOverD0 = TF1('fLcOverD0', '[0]+[1]*TMath::Erfc(x*[2])', 0, 50)

fLcOverD0.SetParameters(0.11, 0.5, 1., 0.5)
fLcOverD0.FixParameter(0, 0.11)
# fLcOverD0.SetParLimits(2, 0.001, 1.e-1)
# fLcOverD0.SetParLimits(3, 0.9, 1.)
SetGraphStyle(fLcOverD0, kAzure+4, -1, 3)

cRatios = TCanvas('cRatios', '', 1920, 1080)
cRatios.Divide(2, 2)
hFrameDplus = cRatios.cd(1).DrawFrame(
    0., 0., 50., 1.5, ';#it{p}_{T} (GeV/#it{c}); D^{+} / D^{0}')
hFrameDplus.GetYaxis().SetDecimals()
gDplusOverD0.Draw('p')
gDplusOverD0.Fit('fDplusOverD0', 'rq')
hFrameDstar = cRatios.cd(2).DrawFrame(
    0., 0., 50., 1.5, ';#it{p}_{T} (GeV/#it{c}); D*^{+} / D^{0}')
hFrameDstar.GetYaxis().SetDecimals()
gDstarOverD0.Draw('p')
gDstarOverD0.Fit('fDstarOverD0', 'rq')
hFrameDs = cRatios.cd(3).DrawFrame(
    0., 0., 50., 1.5, ';#it{p}_{T} (GeV/#it{c}); D_{s}^{+} / D^{0}')
hFrameDs.GetYaxis().SetDecimals()
gDsOverD0.Draw('p')
gDsOverD0.Fit('fDsOverD0', 'rq')
hFrameLc = cRatios.cd(4).DrawFrame(
    0., 0., 50., 1.5, ';#it{p}_{T} (GeV/#it{c}); #Lambda_{c}^{+} / D^{0}')
hFrameLc.GetYaxis().SetDecimals()
gLcOverD0.Draw('p')
gLcOverD0.Fit('fLcOverD0', 'r')

dfDplus = dfD0FONLL.copy()
for col in dfDplus.columns:
    if col not in ['ptmin', 'ptmax']:
        dfDplus[col] = dfDplus.apply(lambda x: x[col]*gDplusOverD0.Eval((x['ptmax']+x['ptmin'])/2), axis=1)

dfDs = dfD0FONLL.copy()
for col in dfDs.columns:
    if col not in ['ptmin', 'ptmax']:
        dfDs[col] = dfDs.apply(lambda x: x[col]*gDsOverD0.Eval((x['ptmax']+x['ptmin'])/2), axis=1)

dfLc = dfD0FONLL.copy()
for col in dfLc.columns:
    if col not in ['ptmin', 'ptmax']:
        dfLc[col] = dfLc.apply(lambda x: x[col]*gLcOverD0.Eval((x['ptmax']+x['ptmin'])/2), axis=1)

dfDplus.to_csv('models/FONLLDataDriven-Dplus-pp14TeV-y05-measbins.txt', sep=' ')
dfDs.to_csv('models/FONLLDataDriven-Ds-pp14TeV-y05-measbins.txt', sep=' ')
dfLc.to_csv('models/FONLLDataDriven-Lc-pp14TeV-y05-measbins.txt', sep=' ')

input('Press enter to exit')
