import argparse
import yaml
from ROOT import TFile, TGraphAsymmErrors, TCanvas
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

hSignificanceCent, hSignificanceMin, hSignificanceMax, gSignificance = ([] for iList in range(4))

for indir, col, mark in zip(cfg['inputdirs'], cfg['colors'], cfg['markers']):

    infileSignif = TFile.Open('{0}/Significance.root'.format(indir))
    hSignificanceCent.append(infileSignif.Get('hSignificance'))
    hSignificanceMin.append(infileSignif.Get('hSignificanceMin'))
    hSignificanceMax.append(infileSignif.Get('hSignificanceMax'))
    hSignificanceCent[-1].SetDirectory(0)
    hSignificanceMin[-1].SetDirectory(0)
    hSignificanceMax[-1].SetDirectory(0)

    gSignificance.append(getGraphFromHistos(hSignificanceCent[-1], hSignificanceMin[-1], hSignificanceMax[-1]))
    SetObjectStyle(gSignificance[-1], color=col, fillalpha=0.3)
    SetObjectStyle(hSignificanceCent[-1], color=col, markerstyle=mark)

cSignificanceVsPt = TCanvas('cSignificanceVsPt', '', 500, 500)
hFrameVsPt = cSignificanceVsPt.DrawFrame(hSignificanceCent[0].GetBinLowEdge(1), hSignificanceMin[0].GetMinimum()*0.5, \
    hSignificanceCent[0].GetXaxis().GetBinUpEdge(hSignificanceCent[0].GetNbinsX()), hSignificanceMax[0].GetMaximum()*2, \
        ';#it{p}_{T} (GeV/#it{c});significance')
hFrameVsPt.GetYaxis().SetTitleOffset(1.2)
cSignificanceVsPt.SetLogy()
for hist, gr in zip(hSignificanceCent, gSignificance):
    gr.Draw('2')
    hist.DrawCopy('same')

input('Press enter to exit')