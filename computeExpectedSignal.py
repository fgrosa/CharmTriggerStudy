'''
Script for the computation of the expected signal
'''

import argparse
import numpy as np
import pandas as pd
import uproot
import yaml
from ROOT import TFile
from dfUtils import writeTH1


def computeExpSignal(ptMin, ptMax, yMin, yMax, prodCrossSec, effAcc, intLumi, BR, fF, fPrompt):
    '''
    Method for the computation of the expected signal
    '''
    return 2 * effAcc * (ptMax - ptMin) * (yMax - yMin) * intLumi * prodCrossSec * BR * fF/ fPrompt


parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='config.yml',
                    help='input root file name')
args = parser.parse_args()

with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)
indir = cfg['output']['directory']
ptmins = cfg['ptmin']
ptmaxs = cfg['ptmax']
ptlims = ptmins
ptlims.append(ptmaxs[-1])
BR = cfg['BR']
fF = cfg['fragfrac']
intLumi = cfg['luminosity']
fPrompt = cfg['fprompt']
model = cfg['model']
dfModel = pd.read_csv(model, sep=' ')
effAccHisto = uproot.open('{0}/Efficiency.root'.format(indir))['hEffAccPrompt']
effAccUncHisto = uproot.open('{0}/Efficiency.root'.format(indir))['hEffAccUncPrompt']
inFileBkg = TFile.Open(f'{indir}/Background.root')

outFile = TFile('{0}/Signal.root'.format(indir), 'recreate')
expSgn, expSgnMin, expSgnMax = [], [], []
expSgnUnc, expSgnUncMin, expSgnUncMax = [], [], []
for (ptMin, ptMax, effAcc, effAccUnc, crossSec, crossSecMin, crossSecMax) in zip(\
    ptmins, ptmaxs, effAccHisto.values, effAccUncHisto.values, \
        dfModel['central'].values, dfModel['min'].values, dfModel['max'].values):
    #cross section already multiplied by DeltaPt
    expSgn.append(computeExpSignal(0, 1, -0.5, 0.5, crossSec, effAcc, intLumi, BR, fF, fPrompt))
    expSgnMin.append(computeExpSignal(0, 1, -0.5, 0.5, crossSecMin, effAcc, intLumi, BR, fF, fPrompt))
    expSgnMax.append(computeExpSignal(0, 1, -0.5, 0.5, crossSecMax, effAcc, intLumi, BR, fF, fPrompt))
    expSgnUnc.append(expSgn[-1]*effAccUnc/effAcc)
    expSgnUncMin.append(expSgnMin[-1]*effAccUnc/effAcc)
    expSgnUncMax.append(expSgnMax[-1]*effAccUnc/effAcc)
    fMass = inFileBkg.Get(f'funcMass_pT_{ptMin:.0f}_{ptMax:.0f}')
    hMass = inFileBkg.Get(f'hMass_pT_{ptMin:.0f}_{ptMax:.0f}')
    hMassScaled = hMass.Clone(f'hMassSignalScaled_pT_{ptMin:.0f}_{ptMax:.0f}')
    hMassScaledMin = hMass.Clone(f'hMassSignalScaledMin_pT_{ptMin:.0f}_{ptMax:.0f}')
    hMassScaledMax = hMass.Clone(f'hMassSignalScaledMax_pT_{ptMin:.0f}_{ptMax:.0f}')
    hMass.SetDirectory(0)
    hMassScaled.SetDirectory(0)
    hMassScaledMin.SetDirectory(0)
    hMassScaledMax.SetDirectory(0)
    for iBin in range(1, hMassScaled.GetNbinsX()+1):
        minMass = hMassScaled.GetBinLowEdge(iBin)
        maxMass = minMass+hMassScaled.GetBinWidth(iBin)
        expFracCounts = fMass.Integral(minMass, maxMass)/fMass.Integral(0., 10)
        hMassScaled.SetBinContent(iBin, expFracCounts * expSgn[-1])
        hMassScaledMin.SetBinContent(iBin, expFracCounts * expSgnMin[-1])
        hMassScaledMax.SetBinContent(iBin, expFracCounts * expSgnMax[-1])
        hMassScaled.SetBinError(iBin, np.sqrt(hMassScaled.GetBinContent(iBin)))
        hMassScaledMin.SetBinError(iBin, np.sqrt(hMassScaledMin.GetBinContent(iBin)))
        hMassScaledMax.SetBinError(iBin, np.sqrt(hMassScaledMax.GetBinContent(iBin)))
    outFile.cd()
    hMassScaled.Write()
    hMassScaledMin.Write()
    hMassScaledMax.Write()

outFile.Close()

writeTH1('{0}/Signal.root'.format(indir), 'hSignal', ptlims, expSgn, expSgnUnc, False)
writeTH1('{0}/Signal.root'.format(indir), 'hSignalMin', ptlims, expSgnMin, expSgnUncMin, False)
writeTH1('{0}/Signal.root'.format(indir), 'hSignalMax', ptlims, expSgnMax, expSgnUncMax, False)
writeTH1('{0}/Signal.root'.format(indir), 'hSignalUnc', ptlims, expSgnUnc, None, False)
writeTH1('{0}/Signal.root'.format(indir), 'hSignalUncMin', ptlims, expSgnUncMin, None, False)
writeTH1('{0}/Signal.root'.format(indir), 'hSignalUncMax', ptlims, expSgnUncMax, None, False)
