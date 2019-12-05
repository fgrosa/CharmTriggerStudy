import argparse
import numpy as np
import yaml
from ROOT import TFile, TGraphAsymmErrors, TCanvas, TLatex, TLegend, TH1F, TLine
from StyleFormatter import SetGlobalStyle, SetObjectStyle

def getGraphFromHistos(histoCent, histoMin, histoMax):
    '''
    Method that takes in input 3 histos and returns a graph
    '''
    nPoints = histoCent.GetNbinsX()
    graph = TGraphAsymmErrors(nPoints)
    for iPt in range(nPoints):
        centC = histoCent.GetBinContent(iPt+1)
        minC = histoMin.GetBinContent(iPt+1)
        maxC = histoMax.GetBinContent(iPt+1)
        width = histoCent.GetBinWidth(iPt+1)/2
        graph.SetPoint(iPt, histoCent.GetBinCenter(iPt+1), centC)
        graph.SetPointError(iPt, width, width, centC-minC, maxC-centC)

    return graph

SetGlobalStyle(padbottommargin=0.14, padtopmargin=0.035)

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='config.yml', help='input root file name')
args = parser.parse_args()

with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

luminosity = cfg['luminosity']
ptbinvsLint = cfg['ptbinvsLint']

hSignificanceCent, hSignificanceMin, hSignificanceMax, gSignificance = ([] for iList in range(4))
hSoverBCent, hSoverBMin, hSoverBMax, gSoverB = ([] for iList in range(4))
hEffAcc = []

for indir, col, mark in zip(cfg['inputdirs'], cfg['colors'], cfg['markers']):

    infileSignif = TFile.Open('{0}/Significance.root'.format(indir))
    hSignificanceCent.append(infileSignif.Get('hSignificance'))
    hSignificanceMin.append(infileSignif.Get('hSignificanceMin'))
    hSignificanceMax.append(infileSignif.Get('hSignificanceMax'))
    hSignificanceCent[-1].SetDirectory(0)
    hSignificanceMin[-1].SetDirectory(0)
    hSignificanceMax[-1].SetDirectory(0)

    gSignificance.append(getGraphFromHistos(hSignificanceCent[-1], hSignificanceMin[-1], hSignificanceMax[-1]))
    SetObjectStyle(gSignificance[-1], color=col, fillalpha=0.2)
    SetObjectStyle(hSignificanceCent[-1], color=col, markerstyle=mark)

    infileSoverB = TFile.Open('{0}/SoverB.root'.format(indir))
    hSoverBCent.append(infileSoverB.Get('hSoverB'))
    hSoverBMin.append(infileSoverB.Get('hSoverBMin'))
    hSoverBMax.append(infileSoverB.Get('hSoverBMax'))
    hSoverBCent[-1].SetDirectory(0)
    hSoverBMin[-1].SetDirectory(0)
    hSoverBMax[-1].SetDirectory(0)

    gSoverB.append(getGraphFromHistos(hSoverBCent[-1], hSoverBMin[-1], hSoverBMax[-1]))
    SetObjectStyle(gSoverB[-1], color=col, fillalpha=0.2)
    SetObjectStyle(hSoverBCent[-1], color=col, markerstyle=mark)

    infileEff = TFile.Open('{0}/Efficiency.root'.format(indir))
    hEffAcc.append(infileEff.Get('hEffAccPrompt'))
    hEffAcc[-1].SetDirectory(0)
    SetObjectStyle(hEffAcc[-1], color=col, markerstyle=mark)

stepLumi = [0.1+lum*0.1 for lum in range(2000)]
hSignificanceVsLumi, gSignificanceVsLumi = [], []
for mes, col, mark, hSignif, hSignifMin, hSignifMax in zip(cfg['inputdirs'], cfg['colors'], cfg['markers'], \
    hSignificanceCent, hSignificanceMin, hSignificanceMax):
    gSignificanceVsLumi.append(TGraphAsymmErrors(0))
    hSignificanceVsLumi.append(TH1F("hSignificanceVsLumi{0}".format(mes), "", 2000, 0.05, 200.05))
    hSignificanceVsLumi[-1].SetDirectory(0)
    SetObjectStyle(hSignificanceVsLumi[-1], color=col, markerstyle=mark)
    SetObjectStyle(gSignificanceVsLumi[-1], color=col, fillalpha=0.2)
    for iPoint, Lumi in enumerate(stepLumi):
        signifCent = hSignif.GetBinContent(ptbinvsLint)*np.sqrt(Lumi/luminosity)
        signifCentUnc = hSignif.GetBinError(ptbinvsLint)*np.sqrt(Lumi/luminosity)
        signifMin = hSignifMin.GetBinContent(ptbinvsLint)*np.sqrt(Lumi/luminosity)
        signifMax = hSignifMax.GetBinContent(ptbinvsLint)*np.sqrt(Lumi/luminosity)
        gSignificanceVsLumi[-1].SetPoint(iPoint, Lumi, signifCent)
        gSignificanceVsLumi[-1].SetPointError(iPoint, 0.05, 0.05, signifCent-signifMin, signifMax-signifCent)
        hSignificanceVsLumi[-1].SetBinContent(iPoint+1, signifCent)
        hSignificanceVsLumi[-1].SetBinError(iPoint+1, signifCentUnc)

lat = TLatex()
lat.SetNDC()
lat.SetTextSize(0.05)
lat.SetTextColor(1)
lat.SetTextFont(42)

lineAtFive = TLine(0.1, 5., 200., 5.)
lineAtFive.SetLineColor(1)
lineAtFive.SetLineStyle(9)
lineAtFive.SetLineWidth(2)

gDummy = TGraphAsymmErrors(0)
SetObjectStyle(gDummy, fillcolor=1, fillalpha=0.2, linecolor=0)

legSignif = TLegend(0.4, 0.2, 0.8, 0.4)
legSignif.SetTextSize(0.05)
legSignif.SetFillStyle(0)
legSignifVsLumi = legSignif.Clone()
for (hist, legname) in zip(hEffAcc, cfg['legendnames']):
    legSignif.AddEntry(hist, legname, 'lp')
    legSignifVsLumi.AddEntry(hist, legname, 'l')
legSignif.AddEntry(gDummy, 'FONLL uncertainty', 'f')

legEff = TLegend(0.4, 0.25, 0.8, 0.4)
legEff.SetTextSize(0.05)
legEff.SetFillStyle(0)
for (hist, legname) in zip(hEffAcc, cfg['legendnames']):
    legEff.AddEntry(hist, legname, 'lp')

cSignificanceVsPt = TCanvas('cSignificanceVsPt', '', 500, 500)
hFrameVsPt = cSignificanceVsPt.DrawFrame(hSignificanceCent[0].GetBinLowEdge(1), 1.e-1, \
    hSignificanceCent[0].GetXaxis().GetBinUpEdge(hSignificanceCent[0].GetNbinsX()), 4.e3, \
        ';#it{p}_{T} (GeV/#it{c});Expected significance (3#sigma)')
hFrameVsPt.GetYaxis().SetTitleOffset(1.2)
cSignificanceVsPt.SetLogy()
for hist, gr in zip(hSignificanceCent, gSignificance):
    gr.Draw('2')
    hist.DrawCopy('same')
lat.DrawLatex(0.18, 0.87, 'pp, #sqrt{#it{s}} = 14 TeV, #it{L}_{int} = 200 pb^{-1}')
legSignif.Draw()

cSignificanceVsLumi = TCanvas('cSignificanceVsLumi', '', 500, 500)
hFrameVsLumi = cSignificanceVsLumi.DrawFrame(0.1, 1.e-2, 200, 4.e2, \
    ';#it{L}_{int} (pb^{-1});Expected significance (3#sigma)')
hFrameVsLumi.GetYaxis().SetTitleOffset(1.2)
cSignificanceVsLumi.SetLogy()
cSignificanceVsLumi.SetLogx()
for hist, gr in zip(hSignificanceVsLumi, gSignificanceVsLumi):
    gr.Draw('C4')
lineAtFive.Draw("same")
lat.DrawLatex(0.18, 0.87, 'pp, #sqrt{#it{s}} = 14 TeV, %0.f < #it{p}_{T} < %0.f GeV/#it{c}' \
    % (hSignificanceCent[0].GetBinLowEdge(ptbinvsLint), \
        hSignificanceCent[0].GetBinLowEdge(ptbinvsLint)+hSignificanceCent[0].GetBinWidth(ptbinvsLint)))
legSignifVsLumi.Draw()

cSoverBVsPt = TCanvas('cSoverBVsPt', '', 500, 500)
hFrameSoverBVsPt = cSoverBVsPt.DrawFrame(hSoverBCent[0].GetBinLowEdge(1), 1.e-3, \
    hSoverBCent[0].GetXaxis().GetBinUpEdge(hSoverBCent[0].GetNbinsX()), 1.e2, \
        ';#it{p}_{T} (GeV/#it{c});Expected S/B (3#sigma)')
hFrameSoverBVsPt.GetYaxis().SetTitleOffset(1.2)
cSoverBVsPt.SetLogy()
for hist, gr in zip(hSoverBCent, gSoverB):
    gr.Draw('2')
    hist.DrawCopy('Lsame')
lat.DrawLatex(0.18, 0.87, 'pp, #sqrt{#it{s}} = 14 TeV, #it{L}_{int} = 200 pb^{-1}')
legSignif.Draw()

cEff = TCanvas('cEff', '', 500, 500)
hFrameEffVsPt = cEff.DrawFrame(hEffAcc[0].GetBinLowEdge(1), 1.e-4, \
    hEffAcc[0].GetXaxis().GetBinUpEdge(hEffAcc[0].GetNbinsX()), 1., \
        ';#it{p}_{T} (GeV/#it{c});Acceptance times efficiency')
hFrameEffVsPt.GetYaxis().SetTitleOffset(1.2)
cEff.SetLogy()
for hist in hEffAcc:
    hist.DrawCopy('same')
legEff.Draw()

input('Press enter to exit')