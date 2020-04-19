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

SetGlobalStyle(padbottommargin=0.14, padleftmargin=0.14, padtopmargin=0.08, opttitle=1)

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='config.yml', help='input yaml file name')
args = parser.parse_args()

with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

luminosity = cfg['luminosity']
ptbinvsLint = cfg['ptbinvsLint']
outputDir = cfg['outputdir']

hSignificanceCent, hSignificanceMin, hSignificanceMax, gSignificance = ([] for iList in range(4))
hSoverBCent, hSoverBMin, hSoverBMax, gSoverB = ([] for iList in range(4))
hEffAccPrompt, hEffAccFD = [], []

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
    hEffAccPrompt.append(infileEff.Get('hEffAccPrompt'))
    hEffAccFD.append(infileEff.Get('hEffAccFD'))
    if 'Ds' in indir:
        hEffAccPrompt[-1].Scale(1./2)
        hEffAccFD[-1].Scale(1./2)
    hEffAccPrompt[-1].SetDirectory(0)
    hEffAccFD[-1].SetDirectory(0)
    SetObjectStyle(hEffAccPrompt[-1], linecolor=col, markercolor=col, markerstyle=mark)
    SetObjectStyle(hEffAccFD[-1], linecolor=col, markercolor=col, markerstyle=mark)

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
for (hist, legname) in zip(hEffAccPrompt, cfg['legendnames']):
    legSignif.AddEntry(hist, legname, 'lp')
    legSignifVsLumi.AddEntry(hist, legname, 'l')
legSignif.AddEntry(gDummy, 'FONLL uncertainty', 'f')
legSignifVsLumi.AddEntry(gDummy, 'FONLL uncertainty', 'f')

legEff = TLegend(0.4, 0.2, 0.8, 0.45)
legEff.SetTextSize(0.05)
legEff.SetFillStyle(0)
for (hist, legname) in zip(hEffAccPrompt, cfg['legendnames']):
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
lat.DrawLatex(0.18, 0.89, 'ALICE Upgrade projection')
lat.DrawLatex(0.18, 0.82, 'pp, #sqrt{#it{s}} = 14 TeV, #it{L}_{int} = 200 pb^{-1}')
legSignif.Draw()

cSignificanceVsLumi = TCanvas('cSignificanceVsLumi', '', 500, 500)
hFrameVsLumi = cSignificanceVsLumi.DrawFrame(0.1, 0.01, 200, 4.e2, \
    ';#it{L}_{int} (pb^{-1});Expected significance (3#sigma)')
hFrameVsLumi.GetYaxis().SetTitleOffset(1.2)
cSignificanceVsLumi.SetLogy()
cSignificanceVsLumi.SetLogx()
for hist, gr in zip(hSignificanceVsLumi, gSignificanceVsLumi):
    gr.Draw('C4')
lineAtFive.Draw("same")
lat.DrawLatex(0.18, 0.89, 'ALICE Upgrade projection')
lat.DrawLatex(0.18, 0.82, 'pp, #sqrt{#it{s}} = 14 TeV, %0.f < #it{p}_{T} < %0.f GeV/#it{c}' \
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
lat.DrawLatex(0.18, 0.89, 'ALICE Upgrade projection')
lat.DrawLatex(0.18, 0.82, 'pp, #sqrt{#it{s}} = 14 TeV, #it{L}_{int} = 200 pb^{-1}')
legSignif.Draw()

cEff = TCanvas('cEff', '', 1000, 500)
cEff.Divide(2, 1)
for iHist, (histP, histF) in enumerate(zip(hEffAccPrompt, hEffAccFD)):
    if iHist == 0:
        hFrameEffPromptVsPt = cEff.cd(1).DrawFrame(hEffAccPrompt[0].GetBinLowEdge(1), 2.e-4, \
            hEffAccPrompt[0].GetXaxis().GetBinUpEdge(hEffAccPrompt[0].GetNbinsX()), 1.5, \
            'Prompt;#it{p}_{T} (GeV/#it{c});Acceptance #times efficiency #times 2#it{y}_{fid}')
        hFrameEffPromptVsPt.GetYaxis().SetTitleOffset(1.2)
    cEff.cd(1).SetLogy()
    histP.DrawCopy('Esame')
    legEff.Draw()

    if iHist == 0:
        hFrameEffFDVsPt = cEff.cd(2).DrawFrame(hEffAccPrompt[0].GetBinLowEdge(1), 2.e-4, \
            hEffAccPrompt[0].GetXaxis().GetBinUpEdge(hEffAccPrompt[0].GetNbinsX()), 1.5, \
            'Feed-down;#it{p}_{T} (GeV/#it{c});Acceptance #times efficiency #times 2#it{y}_{fid}')
        hFrameEffFDVsPt.GetYaxis().SetTitleOffset(1.2)
    cEff.cd(2).SetLogy()
    histF.DrawCopy('Esame')
    legEff.Draw()

cSignificanceVsPt.SaveAs('{0}/Significance_vs_Pt.pdf'.format(outputDir))
cSignificanceVsLumi.SaveAs('{0}/Significance_vs_Lumi.pdf'.format(outputDir))
cSoverBVsPt.SaveAs('{0}/SoverB_vs_Pt.pdf'.format(outputDir))
cEff.SaveAs('{0}/Efficiency_vs_Pt.pdf'.format(outputDir))

input('Press enter to exit')