'''
Script for the creation of the dataframes
'''

import os
import time
import argparse
import pandas as pd
import uproot
import yaml
from dfUtils import filterBitDf, flatColumnNames


def divideInChunks(nEv, nEvPerChunk):
    '''
    Method that divides events in chunks
    '''
    nChunks = int(nEv / nEvPerChunk)
    chunks = [nEvPerChunk for chunk in range(nChunks)]
    chunks.append(nEv-nChunks*nEvPerChunk)
    return chunks


parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='config.yml',
                    help='input root file name')
parser.add_argument('--nevperchunk', type=int, default=1000000, \
    help='number of events per chunk')
parser.add_argument('--convertgen', default=True, action='store_true', \
    help='flag to activate gen conversion')
parser.add_argument('--convertrecosignal', default=True, action='store_true', \
    help='flag to activate reco signal conversion')
parser.add_argument('--convertrecobkg', default=True, action='store_true', \
    help='flag to activate reco bkg conversion')
parser.add_argument('--no-convertgen', dest='convertgen', action='store_false')
parser.add_argument('--no-convertrecosignal', dest='convertrecosignal', action='store_false')
parser.add_argument('--no-convertrecobkg', dest='convertrecobkg', action='store_false')
parser.add_argument('--removepartial', default=False, action='store_true', \
    help='flag to activate deletion of partial parquet files')

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
if args.convertgen:
    startTimeGen = time.time()
    treeGen = uproot.open(cfg['signal']['inputfile'])['{0}/fGenTree'.format(cfg['signal']['inputdir'])]
    chunksGen = divideInChunks(treeGen.numentries, args.nevperchunk)
    print('Number of gen chunks to process: ', len(chunksGen))
    nProcessed = 0
    for iChunk, chunk in enumerate(chunksGen):
        if not os.path.exists('{0}/dfGen_chunk{1}.parquet'.format(outdir, iChunk)):
            dfGenPart = treeGen.pandas.df(entrystart=nProcessed, entrystop=nProcessed+chunk)
            dfGenPart.to_parquet('{0}/dfGen_chunk{1}.parquet'.format(outdir, iChunk))
            del dfGenPart
        nProcessed += chunk

    print('Time elapsed to convert gen tree: ', time.time()-startTimeGen)
    startTimeGen = time.time()

    dfGenChunkPromptSelAcc, dfGenChunkFDSelAcc = [], []
    dfGenChunkPromptSelLimAcc, dfGenChunkFDSelLimAcc = [], []
    for iChunk, _ in enumerate(chunksGen):
        dfGenChunk = pd.read_parquet('{0}/dfGen_chunk{1}.parquet'.format(outdir, iChunk))
        dfGenChunk = dfGenChunk.query('abs(zVtxGen) < 10 & {0} > 0'.format(genBranch))
        flatColumnNames(dfGenChunk, genBranch)

        dfGenChunkPrompt = filterBitDf(dfGenChunk, 'fCandType', [2])
        dfGenChunkFD = filterBitDf(dfGenChunk, 'fCandType', [3])
        for decay in selDecays:
            dfGenChunkPromptSel = dfGenChunkPrompt.query('fDecay == {0}'.format(decay))
            dfGenChunkFDSel = dfGenChunkFD.query('fDecay == {0}'.format(decay))

            dfGenChunkPromptSelAcc.append(dfGenChunkPromptSel.query('abs(fY) < 0.5'))
            dfGenChunkFDSelAcc.append(dfGenChunkFDSel.query('abs(fY) < 0.5'))

            dfGenChunkPromptSelLimAcc.append(filterBitDf(dfGenChunkPromptSel, 'fCandType', [4, 5], 'and'))
            dfGenChunkFDSelLimAcc.append(filterBitDf(dfGenChunkFDSel, 'fCandType', [4, 5], 'and'))

            del dfGenChunkPromptSel, dfGenChunkFDSel

        del dfGenChunk, dfGenChunkPrompt, dfGenChunkFD

    dfGenPromptSelAcc = pd.concat(dfGenChunkPromptSelAcc, ignore_index=True)
    dfGenFDSelAcc = pd.concat(dfGenChunkFDSelAcc, ignore_index=True)
    dfGenPromptSelLimAcc = pd.concat(dfGenChunkPromptSelLimAcc, ignore_index=True)
    dfGenFDSelLimAcc = pd.concat(dfGenChunkFDSelLimAcc, ignore_index=True)

    del dfGenChunkPromptSelAcc, dfGenChunkFDSelAcc, dfGenChunkPromptSelLimAcc, dfGenChunkFDSelLimAcc

    if args.removepartial:
        for iChunk, _ in enumerate(chunksGen):
            os.remove('{0}/dfGen_chunk{1}.parquet'.format(outdir, iChunk))

    dfGenPromptSelAcc.to_parquet('{0}/GenAccPrompt.parquet.gzip'.format(outdir), compression='gzip')
    dfGenFDSelAcc.to_parquet('{0}/GenAccFD.parquet.gzip'.format(outdir), compression='gzip')
    dfGenPromptSelLimAcc.to_parquet('{0}/GenLimAccPrompt.parquet.gzip'.format(outdir), compression='gzip')
    dfGenFDSelLimAcc.to_parquet('{0}/GenLimAccFD.parquet.gzip'.format(outdir), compression='gzip')

    del dfGenPromptSelAcc, dfGenFDSelAcc, dfGenPromptSelLimAcc, dfGenFDSelLimAcc

    print('Time elapsed to filter gen tree: ', time.time()-startTimeGen)

#reco signal
if args.convertrecosignal:
    startTimeReco = time.time()
    treeReco = uproot.open(cfg['signal']['inputfile'])['{0}/fRecoTree'.format(cfg['signal']['inputdir'])]
    chunksReco = divideInChunks(treeReco.numentries, args.nevperchunk)
    print('Number of reco signal chunks to process: ', len(chunksReco))
    nProcessed = 0
    for iChunk, chunk in enumerate(chunksReco):
        if not os.path.exists('{0}/dfReco_chunk{1}.parquet'.format(outdir, iChunk)):
            dfRecoPart = treeReco.pandas.df(['zVtxReco', '{0}*'.format(recoBranch)], \
                entrystart=nProcessed, entrystop=nProcessed+chunk)
            dfRecoPart.to_parquet('{0}/dfReco_chunk{1}.parquet'.format(outdir, iChunk))
            del dfRecoPart
        nProcessed += chunk

    print('Time elapsed to convert reco signal tree: ', time.time()-startTimeReco)
    startTimeReco = time.time()

    dfRecoChunkPromptSelFidAcc, dfRecoChunkFDSelFidAcc = [], []
    for iChunk, _ in enumerate(chunksReco):
        dfRecoChunk = pd.read_parquet('{0}/dfReco_chunk{1}.parquet'.format(outdir, iChunk))
        dfRecoChunk = dfRecoChunk.query('abs(zVtxReco) < 10 & {0} > 0'.format(recoBranch))
        flatColumnNames(dfRecoChunk, recoBranch)

        dfRecoChunkPrompt = filterBitDf(dfRecoChunk, 'fCandType', [2])
        dfRecoChunkFD = filterBitDf(dfRecoChunk, 'fCandType', [3])

        for recoSelBit, decay, massName in zip(recoSelBits, selDecays, massNames):
            dfRecoChunkPromptSel = dfRecoChunkPrompt.query('fDecay == {0}'.format(decay))
            dfRecoChunkFDSel = dfRecoChunkFD.query('fDecay == {0}'.format(decay))

            dfRecoChunkPromptSel['InvMass'] = dfRecoChunkPromptSel[massName]
            dfRecoChunkFDSel['InvMass'] = dfRecoChunkFDSel[massName]

            dfRecoChunkPromptSelFidAcc.append(filterBitDf(dfRecoChunkPromptSel, 'fSelBit', [recoSelBit], 'and'))
            dfRecoChunkFDSelFidAcc.append(filterBitDf(dfRecoChunkFDSel, 'fSelBit', [recoSelBit], 'and'))

            del dfRecoChunkPromptSel, dfRecoChunkFDSel

        del dfRecoChunk, dfRecoChunkPrompt, dfRecoChunkFD

    dfRecoPromptSelFidAcc = pd.concat(dfRecoChunkPromptSelFidAcc, ignore_index=True)
    dfRecoFDSelFidAcc = pd.concat(dfRecoChunkFDSelFidAcc, ignore_index=True)

    del dfRecoChunkPromptSelFidAcc, dfRecoChunkFDSelFidAcc

    if args.removepartial:
        for iChunk, _ in enumerate(chunksReco):
            os.remove('{0}/dfReco_chunk{1}.parquet'.format(outdir, iChunk))

    dfRecoPromptSelFidAcc.to_parquet('{0}/RecoPrompt.parquet.gzip'.format(outdir), compression='gzip')
    dfRecoFDSelFidAcc.to_parquet('{0}/RecoFD.parquet.gzip'.format(outdir), compression='gzip')

    del dfRecoPromptSelFidAcc, dfRecoFDSelFidAcc

    print('Time elapsed to filter reco signal tree: ', time.time()-startTimeReco)
    print('Reco signal saved.')

#reco bkg
if args.convertrecobkg:
    startTimeReco = time.time()
    treeRecoBkg = uproot.open(cfg['background']['inputfile'])['{0}/fRecoTree'.format(cfg['background']['inputdir'])]
    chunksReco = divideInChunks(treeRecoBkg.numentries, args.nevperchunk)
    print('Number of reco bkg chunks to process: ', len(chunksReco))
    nProcessed = 0
    for iChunk, chunk in enumerate(chunksReco):
        if not os.path.exists('{0}/dfRecoBkg_chunk{1}.parquet'.format(outdir, iChunk)):
            dfRecoPart = treeRecoBkg.pandas.df(['zVtxReco', '{0}*'.format(recoBranch)], \
                entrystart=nProcessed, entrystop=nProcessed+chunk)
            dfRecoPart.to_parquet('{0}/dfRecoBkg_chunk{1}.parquet'.format(outdir, iChunk))
            del dfRecoPart
        nProcessed += chunk

    print('Time elapsed to convert reco bkg tree: ', time.time()-startTimeReco)
    startTimeReco = time.time()

    dfRecoChunkBkgSelFidAcc = []
    for iChunk, _ in enumerate(chunksReco):
        dfRecoChunkBkg = pd.read_parquet('{0}/dfRecoBkg_chunk{1}.parquet'.format(outdir, iChunk))
        dfRecoChunkBkg = dfRecoChunkBkg.query('abs(zVtxReco) < 10 & {0} > 0'.format(recoBranch))
        flatColumnNames(dfRecoChunkBkg, recoBranch)
        dfRecoChunkBkgSel = filterBitDf(dfRecoChunkBkg, 'fCandType', [1])

        for recoSelBit, massName in zip(recoSelBits, massNames):

            dfRecoChunkBkgSel = filterBitDf(dfRecoChunkBkgSel, 'fSelBit', [recoSelBit], 'and')

            dfRecoChunkBkgSel['InvMass'] = dfRecoChunkBkgSel[massName]
            dfRecoChunkBkgSelFidAcc.append(dfRecoChunkBkgSel)

            del dfRecoChunkBkgSel

        del dfRecoChunkBkg

    dfRecoBkgSelFidAcc = pd.concat(dfRecoChunkBkgSelFidAcc, ignore_index=True)

    del dfRecoChunkBkgSelFidAcc

    if args.removepartial:
        for iChunk, _ in enumerate(chunksReco):
            os.remove('{0}/dfRecoBkg_chunk{1}.parquet'.format(outdir, iChunk))

    dfRecoBkgSelFidAcc.to_parquet('{0}/RecoBkg.parquet.gzip'.format(outdir), compression='gzip')

    del dfRecoBkgSelFidAcc

    print('Time elapsed to filter reco bkg tree: ', time.time()-startTimeReco)
    print('Reco bkg saved.')
