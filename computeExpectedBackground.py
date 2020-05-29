'''
Script for the computation of the expected background
'''
import argparse
import os
import numpy as np
import pandas as pd
import uproot
import yaml
from ROOT import TH1F, TF1, TFile
from dfUtils import applyLinearCuts, writeTH1

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='config.yml',
                    help='input root file name')
args = parser.parse_args()

def getMeanSigmaValues(dfSignal, ptMins, ptMaxs):
    '''
    Method for the computation of the peak position and width
    '''
    means, sigmas, funcMassList, hMassSignalList = [], [], [], []
    for ptMin, ptMax in zip(ptMins, ptMaxs):
        dfPtBin = dfSignal.query('{0} < fPt < {1}'.format(ptMin, ptMax))
        massHist, massBins = np.histogram(dfPtBin['InvMass'], bins=200)
        hMass = TH1F(f'hMass_pT_{ptMin:.0f}_{ptMax:.0f}', '', len(massBins[:-1]), massBins)
        for iBin, counts in enumerate(massHist):
            hMass.SetBinContent(iBin+1, counts)
            hMass.SetBinError(iBin+1, np.sqrt(counts))
        funcMass = TF1(f'funcMass_pT_{ptMin:.0f}_{ptMax:.0f}', 'gaus', 0., 10.)
        hMass.Fit(f'funcMass_pT_{ptMin:.0f}_{ptMax:.0f}', 'Q0')
        means.append(funcMass.GetParameter(1))
        sigmas.append(funcMass.GetParameter(2))
        funcMassList.append(funcMass)
        hMassSignalList.append(hMass)

    return means, sigmas, funcMassList, hMassSignalList

with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)
indir = cfg['output']['directory']
recoBranch = cfg['recobranch']
ptmins = cfg['ptmin']
ptmaxs = cfg['ptmax']
ptlims = [ptmin for ptmin in ptmins]
ptlims.append(ptmaxs[-1])
if cfg['applysel']['doapplysel']:
    dicofcuts = cfg['applysel']['selections']
nExpEv = cfg['sigmaMB'] * cfg['luminosity'] / 3

treeRecoBkg = uproot.open(cfg['background']['inputfile'])['{0}/fRecoTree'.format(cfg['background']['inputdir'])]
dfEv = treeRecoBkg.pandas.df(['zVtxReco'])
dfEv = dfEv.query('abs(zVtxReco) < 10')
nEv = len(dfEv)

dfSignal = pd.read_parquet('{0}/RecoPrompt.parquet.gzip'.format(indir))
dfSignal = dfSignal.append(pd.read_parquet('{0}/RecoFD.parquet.gzip'.format(indir)))

meanList, sigmaList, funcMassList, hMassSignalList = getMeanSigmaValues(dfSignal, ptmins, ptmaxs)

outfile = TFile('{0}/Background.root'.format(indir), 'recreate')
hMassBkg = []

for iBkgFile, bkgFile in enumerate(os.listdir('{0}/RecoBkg.parquet.gzip'.format(indir))):
    dfBkg = pd.read_parquet('{0}/RecoBkg.parquet.gzip/{1}'.format(indir, bkgFile))

    if cfg['applysel']['doapplysel']:
        if iBkgFile == 0:
            dfBkgSel = applyLinearCuts(dfBkg, dicofcuts)[['fPt', 'InvMass']]
        else:
            dfBkgSel = dfBkgSel.append(applyLinearCuts(dfBkg, dicofcuts)[['fPt', 'InvMass']])
    else:
        if iBkgFile == 0:
            dfBkgSel = dfBkg[['fPt', 'InvMass']]
        else:
            dfBkgSel = dfBkgSel.append(dfBkg[['fPt', 'InvMass']])

    del dfBkg

for iPt, (ptMin, ptMax) in enumerate(zip(ptmins, ptmaxs)):
    dfBkgSelPt = dfBkgSel.query('{0} < fPt < {1}'.format(ptMin, ptMax))
    massBkgHist, massBkgBins = np.histogram(dfBkgSelPt['InvMass'], bins=100)

    del dfBkgSelPt

    hMassBkg.append(TH1F('hMassBkg_pT_%0.f_%0.f' % (ptMin, ptMax), '', len(massBkgBins[:-1]), massBkgBins))

    for iBin, counts in enumerate(massBkgHist):
        hMassBkg[iPt].SetBinContent(iBin+1, hMassBkg[iPt].GetBinContent(iBin+1)+counts)
        hMassBkg[iPt].SetBinError(iBin+1, np.sqrt(hMassBkg[iPt].GetBinContent(iBin+1)))


expBkg, expBkgUnc = [], []
for iPt, (hMass, mean, sigma, ptMin, ptMax) in enumerate(zip(hMassBkg, meanList, sigmaList, ptmins, ptmaxs)):
    funcBkgMass = TF1('funcBkgMass', 'pol1', 0., 10.)
    hMass.Fit('funcBkgMass', 'QL')
    bkg = funcBkgMass.Integral(mean-3*sigma, mean+3*sigma)/hMass.GetBinWidth(1)
    bkgrelunc = funcBkgMass.IntegralError(mean-3*sigma, mean+3*sigma)/hMass.GetBinWidth(1)/bkg
    expBkg.append(bkg / nEv * nExpEv)
    expBkgUnc.append(expBkg[-1] * bkgrelunc)
    hMass.Write()
    hMassBkgScaled = hMass.Clone(f'hMassBkgScaled_pT_{ptMin:.0f}_{ptMax:.0f}')
    hMassBkgScaledMin = hMass.Clone(f'hMassBkgScaledMin_pT_{ptMin:.0f}_{ptMax:.0f}')
    hMassBkgScaledMax = hMass.Clone(f'hMassBkgScaledMax_pT_{ptMin:.0f}_{ptMax:.0f}')
    for iBin in range(1, hMassBkgScaled.GetNbinsX()+1):
        minMass = hMassBkgScaled.GetBinLowEdge(iBin)
        maxMass = minMass+hMassBkgScaled.GetBinWidth(iBin)
        expCounts = funcBkgMass.Integral(minMass, maxMass)/hMass.GetBinWidth(1)
        expCountUnc = funcBkgMass.IntegralError(minMass, maxMass)/hMass.GetBinWidth(1)
        hMassBkgScaled.SetBinContent(iBin, expCounts / nEv * nExpEv)
        hMassBkgScaledMin.SetBinContent(iBin, (expCounts-expCountUnc) / nEv * nExpEv)
        hMassBkgScaledMax.SetBinContent(iBin, (expCounts+expCountUnc) / nEv * nExpEv)
        hMassBkgScaled.SetBinError(iBin, np.sqrt(hMassBkgScaled.GetBinContent(iBin)))
        hMassBkgScaledMin.SetBinError(iBin, np.sqrt(hMassBkgScaledMin.GetBinContent(iBin)))
        hMassBkgScaledMax.SetBinError(iBin, np.sqrt(hMassBkgScaledMax.GetBinContent(iBin)))
    hMassBkgScaled.Write()
    hMassBkgScaledMin.Write()
    hMassBkgScaledMax.Write()
    funcMassList[iPt].Write()
    hMassSignalList[iPt].Write()
outfile.Close()

writeTH1('{0}/Background.root'.format(indir), 'hBkg3Sigma', ptlims, expBkg, expBkgUnc, False)
writeTH1('{0}/Background.root'.format(indir), 'hBkgUnc3Sigma', ptlims, expBkgUnc, None, False)
