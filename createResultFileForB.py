import numpy as np
import yaml
import argparse
import pandas as pd
from ROOT import TFile, TDirectoryFile
from StyleFormatter import SetGlobalStyle, SetObjectStyle

hSgn, hSgnMin, hSgnMax, hBkg, hSignif, hSignifMin, hSignifMax, \
    hSoverB, hSoverBMin, hSoverBMax, hEff, hAcc, hEffAcc = ([] for _ in range(13))
mesons = ['Bplus', 'B0']

ptMins = [0., 2., 4., 6., 8., 12., 16.]
ptMaxs = [2., 4., 6., 8., 12., 16., 24.]
outFile = TFile('BmesonEstimates_Run3_pp14TeV.root', 'recreate')

for meson in mesons:
    infileSgn = TFile.Open(f'{meson}/Signal.root')
    hSgn.append(infileSgn.Get('hSignal'))
    hSgnMin.append(infileSgn.Get('hSignalMin'))
    hSgnMax.append(infileSgn.Get('hSignalMax'))
    infileBkg = TFile.Open(f'{meson}/Background.root')
    hBkg.append(infileBkg.Get('hBkg3Sigma'))
    infileSignif = TFile.Open(f'{meson}/Significance.root')
    hSignif.append(infileSignif.Get('hSignificance'))
    hSignifMin.append(infileSignif.Get('hSignificanceMin'))
    hSignifMax.append(infileSignif.Get('hSignificanceMax'))
    infileSoverB = TFile.Open(f'{meson}/SoverB.root')
    hSoverB.append(infileSoverB.Get('hSoverB'))
    hSoverBMin.append(infileSoverB.Get('hSoverBMin'))
    hSoverBMax.append(infileSoverB.Get('hSoverBMax'))
    infileSoverEff = TFile.Open(f'{meson}/Efficiency.root')
    hEff.append(infileSoverEff.Get('hEffPrompt'))
    hAcc.append(infileSoverEff.Get('hAccPrompt'))
    hEffAcc.append(infileSoverEff.Get('hEffAccPrompt'))

    hMassSignal, hMassSignalMin, hMassSignalMax, hMassBkg, hMassBkgMin, hMassBkgMax = ([] for _ in range(6))
    infileSgn.ls()
    for ptMin, ptMax in zip(ptMins, ptMaxs):
        hMassSignal.append(infileSgn.Get(f'hMassSignalScaled_pT_{ptMin:.0f}_{ptMax:.0f}'))
        hMassSignalMin.append(infileSgn.Get(f'hMassSignalScaledMin_pT_{ptMin:.0f}_{ptMax:.0f}'))
        hMassSignalMax.append(infileSgn.Get(f'hMassSignalScaledMax_pT_{ptMin:.0f}_{ptMax:.0f}'))
        hMassBkg.append(infileBkg.Get(f'hMassBkgScaled_pT_{ptMin:.0f}_{ptMax:.0f}'))
        hMassBkgMin.append(infileBkg.Get(f'hMassBkgScaledMin_pT_{ptMin:.0f}_{ptMax:.0f}'))
        hMassBkgMax.append(infileBkg.Get(f'hMassBkgScaledMax_pT_{ptMin:.0f}_{ptMax:.0f}'))
        hMassSignal[-1].SetDirectory(0)
        hMassSignalMin[-1].SetDirectory(0)
        hMassSignalMax[-1].SetDirectory(0)
        hMassBkg[-1].SetDirectory(0)
        hMassBkgMin[-1].SetDirectory(0)
        hMassBkgMax[-1].SetDirectory(0)

    if meson == 'B0':
        hSgn[-1].SetBinError(5, hSgn[-1].GetBinError(4)*1.1)
        hSgnMin[-1].SetBinError(5, hSgnMin[-1].GetBinError(4)*1.1)
        hSgnMax[-1].SetBinError(5, hSgnMax[-1].GetBinError(4)*1.1)

    hSgn[-1].SetDirectory(0)
    hBkg[-1].SetDirectory(0)
    hSgnMin[-1].SetDirectory(0)
    hSgnMax[-1].SetDirectory(0)
    hSignif[-1].SetDirectory(0)
    hSignifMin[-1].SetDirectory(0)
    hSignifMax[-1].SetDirectory(0)
    hSoverB[-1].SetDirectory(0)
    hSoverBMin[-1].SetDirectory(0)
    hSoverBMax[-1].SetDirectory(0)
    hEff[-1].SetDirectory(0)
    hAcc[-1].SetDirectory(0)
    hEffAcc[-1].SetDirectory(0)

    outFile.cd()
    outDir = TDirectoryFile(meson, meson)
    outDir.Write()
    outDir.cd()
    hSgn[-1].Write()
    hSgnMin[-1].Write()
    hSgnMax[-1].Write()
    for iPt, _ in enumerate(hMassSignal):
        hMassSignal[iPt].Write()
        hMassSignalMin[iPt].Write()
        hMassSignalMax[iPt].Write()
    hBkg[-1].Write()
    for iPt, _ in enumerate(hMassSignal):
        hMassBkg[iPt].Write()
        hMassBkgMin[iPt].Write()
        hMassBkgMax[iPt].Write()
    hSignif[-1].Write()
    hSignifMin[-1].Write()
    hSignifMax[-1].Write()
    hSoverB[-1].Write()
    hSoverBMin[-1].Write()
    hSoverBMax[-1].Write()
    hEff[-1].Write()
    hAcc[-1].Write()
    hEffAcc[-1].Write()
    outDir.Close()

outFile.Close()