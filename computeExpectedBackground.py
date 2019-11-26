'''
Script for the computation of the expected background
'''
import argparse
import numpy as np
import pandas as pd
import uproot
import yaml
from ROOT import TH1F, TF1
from dfUtils import applyLinearCuts, writeTH1

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='config.yml',
                    help='input root file name')
args = parser.parse_args()

def getMeanSigmaValues(dfSignal, ptMins, ptMaxs):
    '''
    Method for the computation of the peak position and width
    '''
    means, sigmas = [], []
    for ptMin, ptMax in zip(ptMins, ptMaxs):
        dfPtBin = dfSignal.query('{0} < fPt < {1}'.format(ptMin, ptMax))
        massHist, massBins = np.histogram(dfPtBin['InvMass'], bins=200)
        hMass = TH1F('hMass', '', len(massBins[:-1]), massBins)
        for iBin, counts in enumerate(massHist):
            hMass.SetBinContent(iBin+1, counts)
            hMass.SetBinError(iBin+1, np.sqrt(counts))
        funcMass = TF1('funcMass', 'gaus', 0., 10.)
        hMass.Fit('funcMass', 'Q0')
        means.append(funcMass.GetParameter(1))
        sigmas.append(funcMass.GetParameter(2)/3)
        del hMass

    return means, sigmas

with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)
indir = cfg['output']['directory']
recoBranch = cfg['recobranch']
ptmins = cfg['ptmin']
ptmaxs = cfg['ptmax']
ptlims = ptmins
ptlims.append(ptmaxs[-1])
if cfg['applysel']['doapplysel']:
    dicofcuts = cfg['applysel']['selections']
nExpEv = cfg['sigmaMB'] * cfg['luminosity']

treeRecoBkg = uproot.open(cfg['background']['inputfile'])['{0}/fRecoTree'.format(cfg['background']['inputdir'])]
dfEv = treeRecoBkg.pandas.df(['zVtxReco'])
dfEv = dfEv.query('abs(zVtxReco) < 10')
nEv = len(dfEv)

dfSignalPrompt = pd.read_parquet('{0}/RecoPrompt.parquet.gzip'.format(indir))
meanList, sigmaList = getMeanSigmaValues(dfSignalPrompt, ptmins, ptmaxs)

dfBkg = pd.read_parquet('{0}/RecoBkg.parquet.gzip'.format(indir))

if cfg['applysel']['doapplysel']:
    dfBkg = applyLinearCuts(dfBkg, dicofcuts)

expBkg, expBkgUnc = [], []
for mean, sigma, ptMin, ptMax in zip(meanList, sigmaList, ptmins, ptmaxs):
    dfBkgSel = dfBkg.query('{0} < fPt < {1} & {2} < InvMass < {3}'.format(ptMin, ptMax, mean-3*sigma, mean+3*sigma))
    expBkg.append(len(dfBkgSel) / nEv * nExpEv)
    expBkgUnc.append(np.sqrt(expBkg[-1]))

writeTH1('{0}/Background.root'.format(indir), 'hBkg3Sigma', ptlims, expBkg, expBkgUnc, True)
