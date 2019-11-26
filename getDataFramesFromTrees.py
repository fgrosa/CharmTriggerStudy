'''
Script for the creation of the dataframes
'''

import argparse
import pandas as pd
import uproot
import yaml
from dfUtils import filterBitDf, flatColumnNames

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='config.yml',
                    help='input root file name')
parser.add_argument('--maxentries', type=int, help='max number of entries')
args = parser.parse_args()

with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)
genBranch = cfg['genbranch']
selDecays = cfg['seldecays']
massNames = cfg['massnames']
recoBranch = cfg['recobranch']
recoSelBits = cfg['recoselbits']
fidAccBit = cfg['fidaccbit']
outdir = cfg['output']['directory']

#gen signal
treeGen = uproot.open(cfg['signal']['inputfile'])['{0}/fGenTree'.format(cfg['signal']['inputdir'])]
if args.maxentries:
    dfGen = treeGen.pandas.df(entrystop=args.maxentries)
else:
    dfGen = treeGen.pandas.df()
dfGen = dfGen.query('abs(zVtxGen) < 10 & {0} > 0'.format(genBranch))
flatColumnNames(dfGen, genBranch)

dfGenPrompt = filterBitDf(dfGen, 'fCandType', [2])
dfGenFD = filterBitDf(dfGen, 'fCandType', [3])

dfGenPromptSelAcc, dfGenFDSelAcc = [], []
dfGenPromptSelLimAcc, dfGenFDSelLimAcc = [], []
for decay in selDecays:
    dfGenPromptSel = dfGenPrompt.query('fDecay == {0}'.format(decay))
    dfGenFDSel = dfGenFD.query('fDecay == {0}'.format(decay))

    dfGenPromptSelAcc.append(dfGenPromptSel.query('abs(fY) < 0.5'))
    dfGenFDSelAcc.append(dfGenFDSel.query('abs(fY) < 0.5'))

    dfGenPromptSelLimAcc.append(filterBitDf(dfGenPromptSel, 'fCandType', [4, 5], 'and'))
    dfGenFDSelLimAcc.append(filterBitDf(dfGenFDSel, 'fCandType', [4, 5], 'and'))

dfGenPromptSelAcc = pd.concat(dfGenPromptSelAcc, ignore_index=True)
dfGenFDSelAcc = pd.concat(dfGenFDSelAcc, ignore_index=True)
dfGenPromptSelLimAcc = pd.concat(dfGenPromptSelLimAcc, ignore_index=True)
dfGenFDSelLimAcc = pd.concat(dfGenFDSelLimAcc, ignore_index=True)

dfGenPromptSelAcc.to_parquet('{0}/GenAccPrompt.parquet.gzip'.format(outdir), compression='gzip')
dfGenFDSelAcc.to_parquet('{0}/GenAccFD.parquet.gzip'.format(outdir), compression='gzip')
dfGenPromptSelLimAcc.to_parquet('{0}/GenLimAccPrompt.parquet.gzip'.format(outdir), compression='gzip')
dfGenFDSelLimAcc.to_parquet('{0}/GenLimAccFD.parquet.gzip'.format(outdir), compression='gzip')

del dfGen, dfGenPrompt, dfGenFD, dfGenPromptSelAcc, dfGenFDSelAcc, dfGenPromptSelLimAcc, dfGenFDSelLimAcc
print('Gen signal saved.')

#reco signal
treeReco = uproot.open(cfg['signal']['inputfile'])['{0}/fRecoTree'.format(cfg['signal']['inputdir'])]
if args.maxentries:
    dfReco = treeReco.pandas.df(['zVtxReco', '{0}*'.format(recoBranch)], entrystop=args.maxentries)
else:
    dfReco = treeReco.pandas.df(['zVtxReco', '{0}*'.format(recoBranch)])

dfReco = dfReco.query('abs(zVtxReco) < 10 & {0} > 0'.format(recoBranch))
flatColumnNames(dfReco, recoBranch)

dfRecoPrompt = filterBitDf(dfReco, 'fCandType', [2])
dfRecoFD = filterBitDf(dfReco, 'fCandType', [3])

dfRecoPromptSelFidAcc, dfRecoFDSelFidAcc = [], []
for iMass, (recoSelBit, decay, massName) in enumerate(zip(recoSelBits, selDecays, massNames)):
    dfRecoPromptSel = dfRecoPrompt.query('fDecay == {0}'.format(decay))
    dfRecoFDSel = dfRecoFD.query('fDecay == {0}'.format(decay))

    dfRecoPromptSelFidAcc.append(filterBitDf(dfRecoPromptSel, 'fSelBit', [recoSelBit], 'and'))
    dfRecoFDSelFidAcc.append(filterBitDf(dfRecoFDSel, 'fSelBit', [recoSelBit], 'and'))
    dfRecoPromptSelFidAcc[iMass]['InvMass'] = dfRecoPromptSelFidAcc[iMass][massName]
    dfRecoFDSelFidAcc[iMass]['InvMass'] = dfRecoFDSelFidAcc[iMass][massName]

dfRecoPromptSelFidAcc = pd.concat(dfRecoPromptSelFidAcc, ignore_index=True)
dfRecoFDSelFidAcc = pd.concat(dfRecoFDSelFidAcc, ignore_index=True)

dfRecoPromptSelFidAcc.to_parquet('{0}/RecoPrompt.parquet.gzip'.format(outdir), compression='gzip')
dfRecoFDSelFidAcc.to_parquet('{0}/RecoFD.parquet.gzip'.format(outdir), compression='gzip')

del dfRecoPrompt, dfRecoFD, dfReco, dfRecoPromptSel, dfRecoFDSel, dfRecoPromptSelFidAcc, dfRecoFDSelFidAcc
print('Reco signal saved.')

#reco bkg
treeRecoBkg = uproot.open(cfg['background']['inputfile'])['{0}/fRecoTree'.format(cfg['background']['inputdir'])]
if args.maxentries:
    dfRecoBkg = treeRecoBkg.pandas.df(['zVtxReco', '{0}*'.format(recoBranch)], entrystop=args.maxentries)
else:
    dfRecoBkg = treeRecoBkg.pandas.df(['zVtxReco', '{0}*'.format(recoBranch)])

dfRecoBkg = dfRecoBkg.query('abs(zVtxReco) < 10 & {0} > 0'.format(recoBranch))
flatColumnNames(dfRecoBkg, recoBranch)
dfRecoBkgSel = filterBitDf(dfRecoBkg, 'fCandType', [1])

dfRecoBkgSelFidAcc = []
for iMass, (recoSelBit, massName) in enumerate(zip(recoSelBits, massNames)):
    dfRecoBkgSelFidAcc.append(filterBitDf(dfRecoBkgSel, 'fSelBit', [recoSelBit], 'and'))
    dfRecoBkgSelFidAcc[iMass]['InvMass'] = dfRecoBkgSelFidAcc[iMass][massName]

dfRecoBkgSelFidAcc = pd.concat(dfRecoBkgSelFidAcc, ignore_index=True)

dfRecoBkgSelFidAcc.to_parquet('{0}/RecoBkg.parquet.gzip'.format(outdir), compression='gzip')
print('Reco bkg saved.')
