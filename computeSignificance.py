'''
Script for the computation of the expected signal
'''

import argparse
import numpy as np
import uproot
import yaml
from dfUtils import writeTH1


def computeExpSigificance(signal, background):
    '''
    Method for the computation of the expected significance
    '''
    signif = signal / np.sqrt(signal+background)
    signifunc = signif*np.sqrt((signal + background) / (4.*(signal+background)**2) + (
        background / (signal + background)) * signal / signal / signal)

    return signif, signifunc


parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='config.yml',
                    help='input root file name')
args = parser.parse_args()

with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)
indir = cfg['output']['directory']

signalCentHisto = uproot.open('{0}/Signal.root'.format(indir))['hSignal']
signalMinHisto = uproot.open('{0}/Signal.root'.format(indir))['hSignalMin']
signalMaxHisto = uproot.open('{0}/Signal.root'.format(indir))['hSignalMax']
backgroundHisto = uproot.open(
    '{0}/Background.root'.format(indir))['hBkg3Sigma']

signifCent, signifUncCent = computeExpSigificance(signalCentHisto.values, backgroundHisto.values)
signifMin, signifUncMin = computeExpSigificance(signalMinHisto.values, backgroundHisto.values)
signifMax, signifUncMax = computeExpSigificance(signalMaxHisto.values, backgroundHisto.values)

writeTH1('{0}/Significance.root'.format(indir), 'hSignificance', signalCentHisto.edges, signifCent, signifUncCent, True)
writeTH1('{0}/Significance.root'.format(indir), 'hSignificanceMin', signalCentHisto.edges, signifMin, signifUncMin, False)
writeTH1('{0}/Significance.root'.format(indir), 'hSignificanceMax', signalCentHisto.edges, signifMax, signifUncMax, False)
