'''
Script for the creation of the dataframes
'''

import os
import time
import argparse
import pandas as pd
import uproot
import yaml
import pyarrow as pa
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
                    help='input yaml config name')
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

    for iChunk, _ in enumerate(chunksGen):
        dfGenChunk = pd.read_parquet('{0}/dfGen_chunk{1}.parquet'.format(outdir, iChunk))
        dfGenChunk = dfGenChunk.query('abs(zVtxGen) < 10 & {0} > 0'.format(genBranch))
        flatColumnNames(dfGenChunk, genBranch)

        dfGenChunkPrompt = filterBitDf(dfGenChunk, 'fCandType', [2])
        dfGenChunkFD = filterBitDf(dfGenChunk, 'fCandType', [3])

        print('process chunk number', iChunk)
        for decay in selDecays:
            dfGenChunkPromptSel = dfGenChunkPrompt.query('fDecay == {0}'.format(decay))
            dfGenChunkFDSel = dfGenChunkFD.query('fDecay == {0}'.format(decay))

            # Create a parquet table from your dataframe
            tbGenPromptSelAcc = pa.Table.from_pandas(dfGenChunkPromptSel.query('abs(fY) < 0.5'))
            tbGenFDSelAcc = pa.Table.from_pandas(dfGenChunkFDSel.query('abs(fY) < 0.5'))
            tbGenPromptSelLimAcc = pa.Table.from_pandas(filterBitDf(dfGenChunkPromptSel, 'fCandType', [4, 5], 'and'))
            tbGenFDSelLimAcc = pa.Table.from_pandas(filterBitDf(dfGenChunkFDSel, 'fCandType', [4, 5], 'and'))

            # Write direct to your parquet file
            pa.parquet.write_to_dataset(tbGenPromptSelAcc, root_path='{0}/GenAccPrompt.parquet.gzip'.format(outdir), \
                compression='gzip')
            pa.parquet.write_to_dataset(tbGenFDSelAcc, root_path='{0}/GenAccFD.parquet.gzip'.format(outdir), \
                compression='gzip')
            pa.parquet.write_to_dataset(tbGenPromptSelLimAcc, root_path='{0}/GenLimAccPrompt.parquet.gzip'.format(outdir), \
                compression='gzip')
            pa.parquet.write_to_dataset(tbGenFDSelLimAcc, root_path='{0}/GenLimAccFD.parquet.gzip'.format(outdir), \
                compression='gzip')

            del dfGenChunkPromptSel, dfGenChunkFDSel

        del dfGenChunk, dfGenChunkPrompt, dfGenChunkFD

    if args.removepartial:
        for iChunk, _ in enumerate(chunksGen):
            os.remove('{0}/dfGen_chunk{1}.parquet'.format(outdir, iChunk))

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

    for iChunk, _ in enumerate(chunksReco):
        dfRecoChunk = pd.read_parquet('{0}/dfReco_chunk{1}.parquet'.format(outdir, iChunk))
        dfRecoChunk = dfRecoChunk.query('abs(zVtxReco) < 10 & {0} > 0'.format(recoBranch))
        flatColumnNames(dfRecoChunk, recoBranch)

        dfRecoChunkPrompt = filterBitDf(dfRecoChunk, 'fCandType', [2])
        dfRecoChunkFD = filterBitDf(dfRecoChunk, 'fCandType', [3])

        print('process chunk number', iChunk)
        for recoSelBit, decay, massName in zip(recoSelBits, selDecays, massNames):
            dfRecoChunkPromptSel = dfRecoChunkPrompt.query('fDecay == {0}'.format(decay))
            dfRecoChunkFDSel = dfRecoChunkFD.query('fDecay == {0}'.format(decay))

            dfRecoChunkPromptSel['InvMass'] = dfRecoChunkPromptSel[massName]
            dfRecoChunkFDSel['InvMass'] = dfRecoChunkFDSel[massName]

            # Create a parquet table from your dataframe
            tbRecoPromptSelFidAcc = pa.Table.from_pandas(filterBitDf(dfRecoChunkPromptSel, 'fSelBit', [recoSelBit, fidAccBit], 'and'))
            tbRecoFDSelFidAcc = pa.Table.from_pandas(filterBitDf(dfRecoChunkFDSel, 'fSelBit', [recoSelBit, fidAccBit], 'and'))

            # Write direct to your parquet file
            pa.parquet.write_to_dataset(tbRecoPromptSelFidAcc, root_path='{0}/RecoPrompt.parquet.gzip'.format(outdir), \
                compression='gzip')
            pa.parquet.write_to_dataset(tbRecoFDSelFidAcc, root_path='{0}/RecoFD.parquet.gzip'.format(outdir), \
                compression='gzip')

            del dfRecoChunkPromptSel, dfRecoChunkFDSel

        del dfRecoChunk, dfRecoChunkPrompt, dfRecoChunkFD

    if args.removepartial:
        for iChunk, _ in enumerate(chunksReco):
            os.remove('{0}/dfReco_chunk{1}.parquet'.format(outdir, iChunk))

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

    dfRecoBkgSelFidAcc = []
    for iChunk, _ in enumerate(chunksReco):
        dfRecoChunkBkg = pd.read_parquet('{0}/dfRecoBkg_chunk{1}.parquet'.format(outdir, iChunk))
        dfRecoChunkBkg = dfRecoChunkBkg.query('abs(zVtxReco) < 10 & {0} > 0'.format(recoBranch))
        flatColumnNames(dfRecoChunkBkg, recoBranch)
        dfRecoChunkBkg = filterBitDf(dfRecoChunkBkg, 'fCandType', [1])

        print('process chunk number', iChunk)
        for recoSelBit, massName in zip(recoSelBits, massNames):

            dfRecoChunkBkgSel = filterBitDf(dfRecoChunkBkg, 'fSelBit', [recoSelBit, fidAccBit], 'and')

            dfRecoChunkBkgSel['InvMass'] = dfRecoChunkBkgSel[massName]

            # Create a parquet table from your dataframe
            tbRecoBkgSelFidAcc = pa.Table.from_pandas(dfRecoChunkBkgSel)
            # Write direct to your parquet file
            pa.parquet.write_to_dataset(tbRecoBkgSelFidAcc, root_path='{0}/RecoBkg.parquet.gzip'.format(outdir), \
                compression='gzip')

            del dfRecoChunkBkgSel

        del dfRecoChunkBkg

    if args.removepartial:
        for iChunk, _ in enumerate(chunksReco):
            os.remove('{0}/dfRecoBkg_chunk{1}.parquet'.format(outdir, iChunk))

    print('Time elapsed to filter reco bkg tree: ', time.time()-startTimeReco)
    print('Reco bkg saved.')
