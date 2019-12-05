'''
Script for the computation of the expected signal
'''

import argparse
import numpy as np
import uproot
import yaml
from dfUtils import writeTH1


def computeExpSigificance(signal, background, errsignal=0, errbkg=0):
    '''
    Method for the computation of the expected significance
    '''
    signif = signal / np.sqrt(signal+background)
    errsignalsq = errsignal**2
    errbkgsq = errbkg**2
    signifunc = signif*np.sqrt((errsignalsq + errbkgsq) / (4.*(signal+background)**2) + (
        background / (signal + background)) * errsignalsq / signal / signal)

    return signif, signifunc

def computeExpSoverB(signal, background, errsignal=0, errbkg=0):
    '''
    Method for the computation of the expected S/B ratio
    '''
    SoverB = signal / background
    SoverBunc = np.sqrt((errsignal/signal)**2 + (errbkg/background)**2) * SoverB
    return SoverB, SoverBunc

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
signalUncCentHisto = uproot.open('{0}/Signal.root'.format(indir))['hSignalUnc']
signalUncMinHisto = uproot.open('{0}/Signal.root'.format(indir))['hSignalUncMin']
signalUncMaxHisto = uproot.open('{0}/Signal.root'.format(indir))['hSignalUncMax']
backgroundHisto = uproot.open(
    '{0}/Background.root'.format(indir))['hBkg3Sigma']
backgroundUncHisto = uproot.open(
    '{0}/Background.root'.format(indir))['hBkgUnc3Sigma']

signifCent, signifUncCent = computeExpSigificance(signalCentHisto.values, backgroundHisto.values, \
    signalUncCentHisto.values, backgroundUncHisto.values)
signifMin, signifUncMin = computeExpSigificance(signalMinHisto.values, backgroundHisto.values, \
    signalUncMinHisto.values, backgroundUncHisto.values)
signifMax, signifUncMax = computeExpSigificance(signalMaxHisto.values, backgroundHisto.values, \
    signalUncMaxHisto.values, backgroundUncHisto.values)

SoverBCent, SoverBUncCent = computeExpSoverB(signalCentHisto.values, backgroundHisto.values, \
    signalUncCentHisto.values, backgroundUncHisto.values)
SoverBMin, SoverBUncMin = computeExpSoverB(signalMinHisto.values, backgroundHisto.values, \
    signalUncMinHisto.values, backgroundUncHisto.values)
SoverBMax, SoverBUncMax = computeExpSoverB(signalMaxHisto.values, backgroundHisto.values, \
    signalUncMaxHisto.values, backgroundUncHisto.values)

writeTH1('{0}/Significance.root'.format(indir), 'hSignificance', signalCentHisto.edges, signifCent, signifUncCent, True)
writeTH1('{0}/Significance.root'.format(indir), 'hSignificanceMin', signalCentHisto.edges, signifMin, signifUncMin, False)
writeTH1('{0}/Significance.root'.format(indir), 'hSignificanceMax', signalCentHisto.edges, signifMax, signifUncMax, False)

writeTH1('{0}/Significance.root'.format(indir), 'hSignificanceUnc', signalCentHisto.edges, signifUncCent, None, False)
writeTH1('{0}/Significance.root'.format(indir), 'hSignificanceUncMin', signalCentHisto.edges, signifUncMin, None, False)
writeTH1('{0}/Significance.root'.format(indir), 'hSignificanceUncMax', signalCentHisto.edges, signifUncMax, None, False)

writeTH1('{0}/SoverB.root'.format(indir), 'hSoverB', signalCentHisto.edges, SoverBCent, SoverBUncCent, True)
writeTH1('{0}/SoverB.root'.format(indir), 'hSoverBMin', signalCentHisto.edges, SoverBMin, SoverBUncMin, False)
writeTH1('{0}/SoverB.root'.format(indir), 'hSoverBMax', signalCentHisto.edges, SoverBMax, SoverBUncMax, False)

writeTH1('{0}/SoverB.root'.format(indir), 'hSoverBUnc', signalCentHisto.edges, SoverBUncCent, None, False)
writeTH1('{0}/SoverB.root'.format(indir), 'hSoverBUncMin', signalCentHisto.edges, SoverBUncMin, None, False)
writeTH1('{0}/SoverB.root'.format(indir), 'hSoverBUncMax', signalCentHisto.edges, SoverBUncMax, None, False)
