'''
Script for the visualization of distributions
'''

import argparse
import os
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import yaml
from dfUtils import applyLinearCuts
from StyleFormatter import formataxes

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='config.yml',
                    help='input root file name')
args = parser.parse_args()

with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)
varsToDraw = cfg['varstodraw']
varMins = cfg['varmins']
varMaxs = cfg['varmaxs']
ptMins = cfg['ptmin']
ptMaxs = cfg['ptmax']
indir = cfg['output']['directory']
if cfg['applysel']['doapplysel']:
    dicofcuts = cfg['applysel']['selections']

# signal
print('Reading signal dataframes')
for iPromptFile, promptFile in enumerate(os.listdir('{0}/RecoPrompt.parquet.gzip'.format(indir))):
    if iPromptFile == 0:
        dfSignalPrompt = pd.read_parquet('{0}/RecoPrompt.parquet.gzip/{1}'.format(indir, promptFile))
        if cfg['applysel']['doapplysel']:
            dfSignalPrompt = applyLinearCuts(dfSignalPrompt, dicofcuts)
    else:
        if len(dfSignalPrompt) >= 1000000:
            break
        dfSignalPrompt = dfSignalPrompt.append(pd.read_parquet('{0}/RecoPrompt.parquet.gzip/{1}'.format(\
            indir, promptFile)))
        if cfg['applysel']['doapplysel']:
            dfSignalPrompt = applyLinearCuts(dfSignalPrompt, dicofcuts)

for iFDFile, FDFile in enumerate(os.listdir('{0}/RecoFD.parquet.gzip'.format(indir))):
    if iFDFile == 0:
        dfSignalFD = pd.read_parquet('{0}/RecoFD.parquet.gzip/{1}'.format(indir, FDFile))
        if cfg['applysel']['doapplysel']:
            dfSignalFD = applyLinearCuts(dfSignalFD, dicofcuts)
    else:
        if len(dfSignalFD) >= 1000000:
            break
        dfSignalFD = dfSignalFD.append(pd.read_parquet('{0}/RecoFD.parquet.gzip/{1}'.format(indir, FDFile)))
        if cfg['applysel']['doapplysel']:
            dfSignalFD = applyLinearCuts(dfSignalFD, dicofcuts)

# background
print('Reading background dataframes')
for iBkgFile, bkgFile in enumerate(os.listdir('{0}/RecoBkg.parquet.gzip'.format(indir))):
    if iBkgFile == 0:
        dfBkg = pd.read_parquet('{0}/RecoBkg.parquet.gzip/{1}'.format(indir, bkgFile))
        if cfg['applysel']['doapplysel']:
            dfBkg = applyLinearCuts(dfBkg, dicofcuts)
    else:
        if len(dfBkg) >= 1000000:
            break
        dfBkg = dfBkg.append(pd.read_parquet('{0}/RecoBkg.parquet.gzip/{1}'.format(indir, bkgFile)))
        if cfg['applysel']['doapplysel']:
            dfBkg = applyLinearCuts(dfBkg, dicofcuts)

for (ptMin, ptMax) in zip(ptMins, ptMaxs):
    figVars, axesVars = plt.subplots(ncols=math.ceil(
        len(varsToDraw)/2), nrows=2, figsize=[12.5, 7])
    for (var, ax, varMin, varMax) in zip(varsToDraw, figVars.axes, varMins, varMaxs):
        ax.hist(dfBkg.query('{0} < fPt < {1}'.format(ptMin, ptMax))[var],
                density=True, color='steelblue', alpha=0.5, bins=100, range=(varMin, varMax), label='background')
        if len(dfSignalFD) > 0:
            ax.hist(dfSignalPrompt.query('{0} < fPt < {1}'.format(ptMin, ptMax))[var],
                    density=True, color='tomato', alpha=0.5, bins=100, range=(varMin, varMax), label='prompt')
            ax.hist(dfSignalFD.query('{0} < fPt < {1}'.format(ptMin, ptMax))[var],
                    density=True, color='seagreen', alpha=0.5, bins=100, range=(varMin, varMax), label='FD')
        else:
            ax.hist(dfSignalPrompt.query('{0} < fPt < {1}'.format(ptMin, ptMax))[var],
                    density=True, color='tomato', alpha=0.5, bins=100, range=(varMin, varMax), label='signal')
        ax.set_yscale('log')
        ax.set_xlabel(var)
        ax.set_ylabel('normalised entries')
        ax.grid()
        ax.legend(loc='best')
    plt.tight_layout()
    figVars.savefig(
        '{0}/VarsPt{1:.0f}_{2:.0f}.pdf'.format(indir, ptMin, ptMax))

if indir == 'B0':
    figPtScatter = plt.figure(figsize=(6, 5))
    hPtBvsPtD = plt.hist2d(dfSignalPrompt['fPt'], dfSignalPrompt['fPtD'], cmap='OrRd',
                           range=np.array([(0., 24.), (0., 24.)]), bins=(96, 96), norm=LogNorm(vmin=1.e-1))
    plt.tick_params(axis='both', which='both', direction="in", labelsize=12)
    axesPtBvsPtD = figPtScatter.get_axes()[0]
    formataxes(axesPtBvsPtD, 0., 24., 0., 24.,
               r'$p_{\rm T} ({\rm B^{0}})$ $({\rm GeV}/c)$',
               r'$p_{\rm T} ({\rm D^{-}})$ $({\rm GeV}/c)$',
               5, 1, 5., 1)
    barPtBvsPtB = figPtScatter.colorbar(hPtBvsPtD[3])
    barPtBvsPtB.ax.tick_params(labelsize=14)
    plt.tight_layout()
    figPtScatter.savefig('{0}/PtBvsPtD.pdf'.format(indir))
elif indir == 'Bplus':
    figPtScatter = plt.figure(figsize=(6, 5))
    hPtBvsPtD = plt.hist2d(dfSignalPrompt['fPt'], dfSignalPrompt['fPtD0'], cmap='OrRd',
                           range=np.array([(0., 24.), (0., 24.)]), bins=(96, 96), norm=LogNorm(vmin=1.e-1))
    plt.tick_params(axis='both', which='both', direction="in", labelsize=12)
    axesPtBvsPtD = figPtScatter.get_axes()[0]
    formataxes(axesPtBvsPtD, 0., 24., 0., 24.,
               r'$p_{\rm T} ({\rm B^{+}})$ $({\rm GeV}/c)$',
               r'$p_{\rm T} ({\rm \overline{D}^{0}})$ $({\rm GeV}/c)$',
               5, 1., 5., 1.)
    barPtBvsPtB = figPtScatter.colorbar(hPtBvsPtD[3])
    barPtBvsPtB.ax.tick_params(labelsize=14)
    plt.tight_layout()
    figPtScatter.savefig('{0}/PtBvsPtD.pdf'.format(indir))

plt.show()
