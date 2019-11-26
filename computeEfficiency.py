'''
Script for the computation of efficiencies
'''

import argparse
import numpy as np
import pandas as pd
import yaml
from dfUtils import applyLinearCuts, writeTH1


def computeEffAndUnc(num, den):
    '''
    Method to compute efficiencies and uncertainty
    '''
    eff, unc = [], []
    for iBin, (nNum, nDen) in enumerate(zip(num, den)):
        if nDen > 0:
            eff.append(nNum / nDen)
            unc.append(np.sqrt((eff[iBin] * abs(1 - eff[iBin])) / nDen))
        else:
            eff.append(0.)
            unc.append(0.)

    return eff, unc


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
if cfg['applysel']['doapplysel']:
    dicofcuts = cfg['applysel']['selections']

#gen part
dfGenAccPrompt = pd.read_parquet('{0}/GenAccPrompt.parquet.gzip'.format(indir))
dfGenAccFD = pd.read_parquet('{0}/GenAccFD.parquet.gzip'.format(indir))
dfGenLimAccPrompt = pd.read_parquet('{0}/GenLimAccPrompt.parquet.gzip'.format(indir))
dfGenLimAccFD = pd.read_parquet('{0}/GenLimAccFD.parquet.gzip'.format(indir))

genAccPrompt, _ = np.histogram(dfGenAccPrompt['fPt'], bins=ptlims)
genAccFD, _ = np.histogram(dfGenAccFD['fPt'], bins=ptlims)

genLimAccPrompt, _ = np.histogram(dfGenLimAccPrompt['fPt'], bins=ptlims)
genLimAccFD, _ = np.histogram(dfGenLimAccFD['fPt'], bins=ptlims)

#reco part
dfSignalPrompt = pd.read_parquet('{0}/RecoPrompt.parquet.gzip'.format(indir))
dfSignalFD = pd.read_parquet('{0}/RecoFD.parquet.gzip'.format(indir))

if cfg['applysel']['doapplysel']:
    dfSignalPrompt = applyLinearCuts(dfSignalPrompt, dicofcuts)
    dfSignalFD = applyLinearCuts(dfSignalFD, dicofcuts)

recoPrompt, _ = np.histogram(dfSignalPrompt['fPt'], bins=ptlims)
recoFD, _ = np.histogram(dfSignalFD['fPt'], bins=ptlims)

effAccPrompt, effAccPromptUnc = computeEffAndUnc(recoPrompt, genAccPrompt)
effAccFD, effAccFDUnc = computeEffAndUnc(recoFD, genAccFD)

effPrompt, effPromptUnc = computeEffAndUnc(recoPrompt, genLimAccPrompt)
effFD, effFDUnc = computeEffAndUnc(recoFD, genLimAccFD)

accPrompt, effAccPromptUnc = computeEffAndUnc(genLimAccPrompt, genAccPrompt)
accFD, effAccFDUnc = computeEffAndUnc(genLimAccFD, genAccFD)

writeTH1('{0}/Efficiency.root'.format(indir), 'hEffAccPrompt', ptlims, effAccPrompt, effAccPromptUnc, True)
writeTH1('{0}/Efficiency.root'.format(indir), 'hEffAccFD', ptlims, effAccFD, effAccFDUnc, False)
writeTH1('{0}/Efficiency.root'.format(indir), 'hEffPrompt', ptlims, effPrompt, effPromptUnc, False)
writeTH1('{0}/Efficiency.root'.format(indir), 'hEffFD', ptlims, effFD, effFDUnc, False)
writeTH1('{0}/Efficiency.root'.format(indir), 'hAccPrompt', ptlims, accPrompt, effAccPromptUnc, False)
writeTH1('{0}/Efficiency.root'.format(indir), 'hAccFD', ptlims, accFD, effAccFDUnc, False)
