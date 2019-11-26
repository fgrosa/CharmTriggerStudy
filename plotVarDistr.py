'''
Script for the visualization of distributions
'''

import argparse
import math
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from dfUtils import applyLinearCuts

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

#signal
dfSignalPrompt = pd.read_parquet('{0}/RecoPrompt.parquet.gzip'.format(indir))
dfSignalFD = pd.read_parquet('{0}/RecoFD.parquet.gzip'.format(indir))
#background
dfBkg = pd.read_parquet('{0}/RecoBkg.parquet.gzip'.format(indir))

if cfg['applysel']['doapplysel']:
    dfBkg = applyLinearCuts(dfBkg, dicofcuts)
    dfSignalPrompt = applyLinearCuts(dfSignalPrompt, dicofcuts)
    dfSignalFD = applyLinearCuts(dfSignalFD, dicofcuts)

for (ptMin, ptMax) in zip(ptMins, ptMaxs):
    figVars, axesVars = plt.subplots(ncols=math.ceil(len(varsToDraw)/2), nrows=2, figsize=[12, 7])
    for (var, ax, varMin, varMax) in zip(varsToDraw, figVars.axes, varMins, varMaxs):
        ax.hist(dfBkg.query('{0} < fPt < {1}'.format(ptMin, ptMax))[var], \
            density=True, alpha=0.5, bins=100, range=(varMin, varMax), label='background')
        ax.hist(dfSignalPrompt.query('{0} < fPt < {1}'.format(ptMin, ptMax))[var], \
            density=True, alpha=0.5, bins=100, range=(varMin, varMax), label='prompt')
        ax.hist(dfSignalFD.query('{0} < fPt < {1}'.format(ptMin, ptMax))[var], \
            density=True, alpha=0.5, bins=100, range=(varMin, varMax), label='FD')
        ax.set_yscale('log')
        ax.set_xlabel(var)
        ax.set_ylabel('normalised entries')
        ax.grid()
    plt.tight_layout()
    figVars.savefig('{0}/VarsPt{1:.0f}_{2:.0f}.pdf'.format(indir, ptMin, ptMax))

plt.show()
