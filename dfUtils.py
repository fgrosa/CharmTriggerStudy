import numpy as np
import array
import pandas as pd
from ROOT import TH1F, TFile

def flatten(dfToFlat, colToFlat):
    '''
    Method for df flattening
    '''
    idx = dfToFlat.index.repeat(dfToFlat[colToFlat[0]].str.len())
    dfFlat = pd.concat([
        pd.DataFrame({x: np.concatenate(dfToFlat[x].values)}) for x in colToFlat], axis=1)
    dfFlat.index = idx

    return dfFlat.join(dfToFlat.drop(colToFlat, 1), how='left')


def getMaskOfBits(bits):
    '''
    Method to get bit mask from bits
    '''
    mask = 0
    for bit in bits:
        mask += 2**bit

    return mask


def filterBitDf(dfToFilter, column, bitsToTest, logic='or'):
    '''
    Method to apply selection testing a bit
    '''
    maskOfBits = getMaskOfBits(bitsToTest)
    flags = dfToFilter[column] & maskOfBits
    if logic == 'or':
        flags = flags.astype('bool')
    elif logic == 'and':
        flags -= maskOfBits
        flags = ~flags.astype('bool')
    else:
        print('ERROR: only and and or logics are supported for bitwise operations')
        return None

    dfFilt = dfToFilter[flags.values]

    return dfFilt


def flatColumnNames(dfToRename, branchName):
    '''
    Method to remove struct name from columns
    '''
    for col in dfToRename.columns:
        if '{0}.'.format(branchName) in col:
            colnew = col.replace('{0}.'.format(branchName), '')
            dfToRename.rename(columns={col: colnew}, inplace=True)


def applyLinearCuts(dfToFilter, dicOfCuts):
    '''
    Method to apply linear selections to a pandas dataframe
    '''
    stringForOneBin = {}
    for var in dicOfCuts:
        stringForOneBin[var] = {}
        for iBin, (minVal, maxVal) in enumerate(zip(dicOfCuts[var]['min'], dicOfCuts[var]['max'])):
            stringForOneBin[var][iBin] = '{0} < {1} < {2}'.format(minVal, var, maxVal)


    stringTotCuts = ''
    for iBin in stringForOneBin[list(stringForOneBin.keys())[0]]:
        stringBin = ''
        for iVar, var in enumerate(stringForOneBin):
            if iVar > 0:
                stringBin += ' and '
            stringBin += stringForOneBin[var][iBin]

        if iBin > 0:
            stringTotCuts += ' or '
        stringTotCuts += '({0})'.format(stringBin)

    dfFilt = dfToFilter.query(stringTotCuts)

    return dfFilt

def writeTH1(filename, name, binLims, binContents, binErrors=None, reWriteFile=False):
    '''
    Method to save TH1F
    '''
    hToSave = TH1F(name, '', len(binLims[:-1]), array.array('f', binLims))
    if binErrors is None:
        binErrors = np.sqrt(binContents)
    for iBin, (content, error) in enumerate(zip(binContents, binErrors)):
        hToSave.SetBinContent(iBin+1, content)
        hToSave.SetBinError(iBin+1, error)

    if not reWriteFile and TFile.Open(filename):
        outfile = TFile.Open(filename, 'update')
    else:
        outfile = TFile(filename, 'recreate')
    hToSave.Write(name)
    outfile.Close()
