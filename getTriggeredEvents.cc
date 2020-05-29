#if !defined (__CINT__) || defined (__CLING__)

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "yaml-cpp/yaml.h"

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TDatabasePDG.h>

#include "AliAnalysisTaskSECharmTriggerStudy.h"

#endif

using std::vector;
using std::string;
using std::map;
using std::cout;
using std::cerr;
using std::endl;

const double massD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
const double massDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
const double massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
const double massLc = TDatabasePDG::Instance()->GetParticle(4122)->Mass();

//______________________________________________________________________________________________
bool is2ProngSelected(Charm2Prong ch2Prong, map<string, double> cuts)
{
    if(ch2Prong.fd0MinDau < cuts["fd0MinDau"]) return false;
    if(ch2Prong.fNormDecayLengthXY < cuts["fNormDecayLengthXY"]) return false;
    if(ch2Prong.fDecayLength < cuts["fDecayLength"]) return false;
    if(ch2Prong.fImpParProd > cuts["fImpParProd"]) return false;
    if(ch2Prong.fPtMinDau < cuts["fPtMinDau"]) return false;
    if(ch2Prong.fCosP < cuts["fCosP"])  return false;
    if(ch2Prong.fCosPXY < cuts["fCosPXY"]) return false;
    if(TMath::Abs(ch2Prong.fInvMassD0-massD0) > cuts["fDeltaMass"] && TMath::Abs(ch2Prong.fInvMassD0bar-massD0) > cuts["fDeltaMass"]) return false;
    return true;
}

//______________________________________________________________________________________________
bool is3ProngSelected(Charm3Prong ch3Prong, map<string, double> cuts)
{
    if(ch3Prong.fd0MinDau < cuts["fd0MinDau"]) return false;
    if(ch3Prong.fNormDecayLengthXY < cuts["fNormDecayLengthXY"]) return false;
    if(ch3Prong.fDecayLength < cuts["fDecayLength"]) return false;
    if(ch3Prong.fSigmaVtx > cuts["fSigmaVtx"]) return false;
    if(ch3Prong.fPtMinDau < cuts["fPtMinDau"]) return false;
    if(ch3Prong.fCosP < cuts["fCosP"])  return false;
    if(ch3Prong.fCosPXY < cuts["fCosPXY"]) return false;
    if(TMath::Abs(ch3Prong.fInvMassDplus-massDplus) > cuts["fDeltaMass"] &&
       TMath::Abs(ch3Prong.fInvMassDstoKKpi-massDs) > cuts["fDeltaMass"] && TMath::Abs(ch3Prong.fInvMassDstopiKK-massDs) > cuts["fDeltaMass"] &&
       TMath::Abs(ch3Prong.fInvMassLctopKpi-massLc) > cuts["fDeltaMass"] && TMath::Abs(ch3Prong.fInvMassLctopiKp-massLc) > cuts["fDeltaMass"])
        return false;
    return true;
}

//______________________________________________________________________________________________
bool isBeauty3ProngSelected(Beauty3Prong beauty3Prong, map<string, double> cuts)
{
    if(beauty3Prong.fd0MinDauD0 < cuts["fd0MinDauD"]) return false;
    if(beauty3Prong.fNormDecayLengthXYD0 < cuts["fNormDecayLengthXYD"]) return false;
    if(beauty3Prong.fDecayLengthD0 < cuts["fDecayLengthD"]) return false;
    if(beauty3Prong.fImpParProdD0 > cuts["fImpParProdD"]) return false;
    if(beauty3Prong.fPtMinDauD0 < cuts["fPtMinDauD"]) return false;
    if(beauty3Prong.fCosPD0 < cuts["fCosPD"])  return false;
    if(beauty3Prong.fCosPXYD0 < cuts["fCosPXYD"]) return false;
    if(TMath::Abs(beauty3Prong.fInvMassD0-massD0) > cuts["fDeltaMassD"]) return false;
    return true;
}

//______________________________________________________________________________________________
bool isBeauty4ProngSelected(Beauty4Prong beauty4Prong, map<string, double> cuts)
{
    if(beauty4Prong.fd0MinDauD < cuts["fd0MinDauD"]) return false;
    if(beauty4Prong.fNormDecayLengthXYD < cuts["fNormDecayLengthXYD"]) return false;
    if(beauty4Prong.fDecayLengthD < cuts["fDecayLengthD"]) return false;
    if(beauty4Prong.fSigmaVtxD > cuts["fSigmaVtxD"]) return false;
    if(beauty4Prong.fPtMinDauD < cuts["fPtMinDauD"]) return false;
    if(beauty4Prong.fCosPD < cuts["fCosPD"])  return false;
    if(beauty4Prong.fCosPXYD < cuts["fCosPXYD"]) return false;
    if((!(beauty4Prong.fSelBit&BIT(10)) || (beauty4Prong.fSelBit&BIT(10) && TMath::Abs(beauty4Prong.fInvMassDplus-massDplus) > cuts["fDeltaMassD"])) && 
       (!(beauty4Prong.fSelBit&BIT(11)) || (beauty4Prong.fSelBit&BIT(11) && TMath::Abs(beauty4Prong.fInvMassDs-massDs) > cuts["fDeltaMassD"])) &&
       (!(beauty4Prong.fSelBit&BIT(12)) || (beauty4Prong.fSelBit&BIT(12) && TMath::Abs(beauty4Prong.fInvMassLc-massLc) > cuts["fDeltaMassD"])))
       return false;
    return true;
}

//______________________________________________________________________________________________
bool is2ProngSelectedPID(Charm2Prong ch2Prong)
{
    if(ch2Prong.fSelBit&BIT(13) || ch2Prong.fSelBit&BIT(14)) return true;
    return false;
}

//______________________________________________________________________________________________
bool is3ProngSelectedPID(Charm3Prong ch3Prong)
{
    if(ch3Prong.fSelBit&BIT(15) || ch3Prong.fSelBit&BIT(17) || ch3Prong.fSelBit&BIT(18) || ch3Prong.fSelBit&BIT(19) || ch3Prong.fSelBit&BIT(20)) return true;
    return false;
}

//______________________________________________________________________________________________
vector<map<string, double> > buildMapOfCuts(YAML::Node cuts, vector<double> &PtLims)
{
    PtLims = cuts["fPtLims"].as<vector<double> >();
    unsigned int nPtLims = PtLims.size();
    vector<map<string, double> > mapOfCuts;
    for(unsigned int iPt=0; iPt<nPtLims-1; iPt++)
    {
        map<string, double> mapPtBin;
        for(auto var : cuts)
        {
            if(var.first.as<string>() == "fPtLims")
                continue;
            mapPtBin[var.first.as<string>()] = var.second.as<vector<double> >()[iPt];
        }
        mapOfCuts.push_back(mapPtBin);
    }
    return mapOfCuts;
}

//______________________________________________________________________________________________
int findPtBin(double Pt, vector<double> PtLims)
{
    int bin = TMath::BinarySearch(PtLims.size() - 1, &PtLims[0], Pt);
    if (bin < 0) //underflow --> equal to min value
        bin = 0;

    return bin;
}

//______________________________________________________________________________________________
bool AreIdxInCommon2Prong(int idx0first, int idx1first, int idx0second, int idx1second)
{
    if(idx0first == idx0second || idx1first == idx1second || idx0first == idx1second || idx1first == idx0second)
        return true;
    
    return false;
}

//______________________________________________________________________________________________
bool AreIdxInCommon3Prong(int idx0first, int idx1first, int idx2first, int idx0second, int idx1second, int idx2second)
{
    if(idx0first == idx0second || idx1first == idx1second || idx2first == idx2second || 
       idx0first == idx1second || idx0first == idx2second || 
       idx1first == idx0second || idx1first == idx2second ||
       idx2first == idx0second || idx2first == idx1second)
        return true;
    
    return false;
}

//______________________________________________________________________________________________
bool AreIdxInCommon2Prong3Prong(int idx0first, int idx1first, int idx0second, int idx1second, int idx2second)
{
    if(idx0first == idx0second || idx1first == idx1second || 
       idx0first == idx1second || idx0first == idx2second || 
       idx1first == idx0second || idx1first == idx2second)
        return true;
    
    return false;
}

//______________________________________________________________________________________________
void getTriggeredEventsWithD(TString cfgFileName="config.yml")
{

    //Load configs from yaml file
    YAML::Node config = YAML::LoadFile(cfgFileName.Data());
    if (config.IsNull()) {
        cerr << "Yaml config file not found! Exit" << endl;
        return;
    }

    string inFileName = config["infile"]["path"].as<string>();
    string inDirName = config["infile"]["dir"].as<string>();
    string outFileName = config["outfile"].as<string>();

    int nMaxEvents = config["nevents"].as<int>();

    vector<double> PtThresholds = config["ptthresholds"].as<vector<double> >();
    const unsigned int nPtThresholds = PtThresholds.size();
    
    bool applyCuts = static_cast<bool>(config["applycuts"].as<int>());
    bool applyPID = static_cast<bool>(config["applyPID"].as<int>());

    vector<double> PtLims2Prong, PtLims3Prong;
    vector<map<string, double> > mapOfCuts2Prongs = buildMapOfCuts(config["cuts2Prong"], PtLims2Prong);
    vector<map<string, double> > mapOfCuts3Prongs = buildMapOfCuts(config["cuts3Prong"], PtLims3Prong);

    //start analysis
    TFile* infile = TFile::Open(inFileName.data());
    TTree* evTree = dynamic_cast<TTree*>(infile->Get(Form("%s/fRecoTree", inDirName.data())));
    vector<Charm2Prong> *ch2ProngVec = nullptr;
    vector<Charm3Prong> *ch3ProngVec = nullptr;
    float zVtx = -999.;
    evTree->SetBranchAddress("Charm2Prong", &ch2ProngVec);
    evTree->SetBranchAddress("Charm3Prong", &ch3ProngVec);
    evTree->SetBranchAddress("zVtxReco", &zVtx);

    int nTriggeredEvents2Prong[nPtThresholds];
    int nTriggeredEvents3Prong[nPtThresholds];
    int nTriggeredEventsDouble2Prong[nPtThresholds];
    int nTriggeredEventsDouble3Prong[nPtThresholds];
    int nTriggeredEvents2ProngAnd3Prong[nPtThresholds];
    int nTriggeredEvents2ProngOr3Prong[nPtThresholds];

    int nTriggeredSignalEvents2Prong[nPtThresholds];
    int nTriggeredSignalEvents3Prong[nPtThresholds];
    int nTriggeredSignalEventsDouble2Prong[nPtThresholds];
    int nTriggeredSignalEventsDouble3Prong[nPtThresholds];
    int nTriggeredSignalEvents2ProngAnd3Prong[nPtThresholds];
    int nTriggeredSignalEvents2ProngOr3Prong[nPtThresholds];

    for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
    {
        nTriggeredEvents2Prong[iPtThr] = 0;
        nTriggeredEvents3Prong[iPtThr] = 0;
        nTriggeredEventsDouble2Prong[iPtThr] = 0;
        nTriggeredEventsDouble3Prong[iPtThr] = 0;
        nTriggeredEvents2ProngAnd3Prong[iPtThr] = 0;
        nTriggeredEvents2ProngOr3Prong[iPtThr] = 0;

        nTriggeredSignalEvents2Prong[iPtThr] = 0;
        nTriggeredSignalEvents3Prong[iPtThr] = 0;
        nTriggeredSignalEventsDouble2Prong[iPtThr] = 0;
        nTriggeredSignalEventsDouble3Prong[iPtThr] = 0;
        nTriggeredSignalEvents2ProngAnd3Prong[iPtThr] = 0;
        nTriggeredSignalEvents2ProngOr3Prong[iPtThr] = 0;
    }
    int nEvents = (nMaxEvents < evTree->GetEntriesFast()) ? nMaxEvents : evTree->GetEntriesFast();

    for(int iEv=0; iEv<nEvents; iEv++)
    {
        if((iEv+1)%100000 == 0)
            cout << "Processing event " << iEv+1 << endl;

        evTree->GetEvent(iEv);
        if(ch2ProngVec->size() == 0 && ch3ProngVec->size() == 0)
            continue;

        int n2Prong[nPtThresholds];
        int n3Prong[nPtThresholds];
        int n2ProngSignal[nPtThresholds];
        int n3ProngSignal[nPtThresholds];
        for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
        {
            n2Prong[iPtThr] = 0;
            n3Prong[iPtThr] = 0;
            n2ProngSignal[iPtThr] = 0;
            n3ProngSignal[iPtThr] = 0;
        }
        vector<int> idx2Prong0[nPtThresholds];
        vector<int> idx2Prong1[nPtThresholds];
        vector<int> idx3Prong0[nPtThresholds];
        vector<int> idx3Prong1[nPtThresholds];
        vector<int> idx3Prong2[nPtThresholds];

        //loop over 2 prongs
        for(unsigned int i2Prong=0; i2Prong<ch2ProngVec->size(); i2Prong++)
        {
            float pt = ch2ProngVec->at(i2Prong).fPt;
            int ptBin = findPtBin(pt, PtLims2Prong);
            bool isSignal = ch2ProngVec->at(i2Prong).fCandType&BIT(0);
            
            for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
            {
                if(pt > PtThresholds[iPtThr] && (!applyCuts || is2ProngSelected(ch2ProngVec->at(i2Prong), mapOfCuts2Prongs[ptBin])) && (!applyPID || is2ProngSelectedPID(ch2ProngVec->at(i2Prong))))
                {
                    n2Prong[iPtThr]++;
                    if(isSignal)
                        n2ProngSignal[iPtThr]++;
                    idx2Prong0[iPtThr].push_back(ch2ProngVec->at(i2Prong).fProngIdx0);
                    idx2Prong1[iPtThr].push_back(ch2ProngVec->at(i2Prong).fProngIdx1);
                }
            }
        }

        //loop over 3 prongs
        for(unsigned int i3Prong=0; i3Prong<ch3ProngVec->size(); i3Prong++)
        {
            float pt = ch3ProngVec->at(i3Prong).fPt;
            int ptBin = findPtBin(pt, PtLims3Prong);
            bool isSignal = ch3ProngVec->at(i3Prong).fCandType&BIT(0);

            for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
            {
                if(pt > PtThresholds[iPtThr] && (!applyCuts || is3ProngSelected(ch3ProngVec->at(i3Prong), mapOfCuts3Prongs[ptBin])) && (!applyPID || is3ProngSelectedPID(ch3ProngVec->at(i3Prong))))
                {
                    n3Prong[iPtThr]++;
                    if(isSignal)
                        n3ProngSignal[iPtThr]++;
                    idx3Prong0[iPtThr].push_back(ch3ProngVec->at(i3Prong).fProngIdx0);
                    idx3Prong1[iPtThr].push_back(ch3ProngVec->at(i3Prong).fProngIdx1);
                    idx3Prong2[iPtThr].push_back(ch3ProngVec->at(i3Prong).fProngIdx2);
                }
            }
        }

        for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
        {
            if(n2Prong[iPtThr] > 0)
            {
                nTriggeredEvents2Prong[iPtThr]++;
                if(n2ProngSignal[iPtThr] > 0)
                    nTriggeredSignalEvents2Prong[iPtThr]++;

                if(n2Prong[iPtThr] > 1)
                {
                    //check that the pairs of candidates have no daughters in common
                    int nProngsNotInCommon2Prong = 0;

                    for(int i2ProngFirst=0; i2ProngFirst<n2Prong[iPtThr]-1; i2ProngFirst++)
                    {
                        for(int i2ProngSecond=i2ProngFirst+1; i2ProngSecond<n2Prong[iPtThr]; i2ProngSecond++)
                        {
                            if(!AreIdxInCommon2Prong(idx2Prong0[iPtThr][i2ProngFirst], idx2Prong1[iPtThr][i2ProngFirst], idx2Prong0[iPtThr][i2ProngSecond], idx2Prong1[iPtThr][i2ProngSecond]))
                            nProngsNotInCommon2Prong++;
                        }
                    }

                    if(nProngsNotInCommon2Prong > 0)
                    {
                        nTriggeredEventsDouble2Prong[iPtThr]++;
                        if(n2ProngSignal[iPtThr] > 1)
                            nTriggeredSignalEventsDouble2Prong[iPtThr]++;
                    }
                }
            }

            if(n3Prong[iPtThr] > 0)
            {
                //check that the pairs of candidates have no daughters in common
                int nProngsNotInCommon3Prong = 0;

                for(int i3ProngFirst=0; i3ProngFirst<n3Prong[iPtThr]-1; i3ProngFirst++)
                {
                    for(int i3ProngSecond=i3ProngFirst+1; i3ProngSecond<n3Prong[iPtThr]; i3ProngSecond++)
                    {
                        if(!AreIdxInCommon3Prong(idx3Prong0[iPtThr][i3ProngFirst], idx3Prong1[iPtThr][i3ProngFirst], idx3Prong2[iPtThr][i3ProngFirst], idx3Prong0[iPtThr][i3ProngSecond], idx3Prong1[iPtThr][i3ProngSecond], idx3Prong2[iPtThr][i3ProngSecond]))
                           nProngsNotInCommon3Prong++;
                    }
                }

                nTriggeredEvents3Prong[iPtThr]++;
                 if(n3ProngSignal[iPtThr] > 0)
                    nTriggeredSignalEvents3Prong[iPtThr]++;

                if(nProngsNotInCommon3Prong > 0)
                {
                    if(n3Prong[iPtThr] > 1)
                    {
                        nTriggeredEventsDouble3Prong[iPtThr]++;
                        if(n3ProngSignal[iPtThr] > 0)
                            nTriggeredSignalEvents3Prong[iPtThr]++;
                    }
                }
            }

            if(n2Prong[iPtThr] > 0 || n3Prong[iPtThr] > 0)
            {
                nTriggeredEvents2ProngOr3Prong[iPtThr]++;
                 if(n2ProngSignal[iPtThr] > 0 || n3ProngSignal[iPtThr] > 0)
                    nTriggeredSignalEvents2ProngOr3Prong[iPtThr]++;    
            }

            if(n2Prong[iPtThr] > 0 && n3Prong[iPtThr] > 0)
            {
                //check that the pairs of candidates have no daughters in common
                int nProngsNotInCommon2Prong3Prong = 0;

                for(int i2Prong=0; i2Prong<n2Prong[iPtThr]; i2Prong++)
                {
                    for(int i3Prong=0; i3Prong<n3Prong[iPtThr]; i3Prong++)
                    {
                        if(!AreIdxInCommon2Prong3Prong(idx2Prong0[iPtThr][i2Prong], idx2Prong1[iPtThr][i3Prong], idx3Prong0[iPtThr][i3Prong], idx3Prong1[iPtThr][i3Prong], idx3Prong2[iPtThr][i3Prong]))
                            nProngsNotInCommon2Prong3Prong++;
                    }
                }

                if(nProngsNotInCommon2Prong3Prong > 0)
                {
                    nTriggeredEvents2ProngAnd3Prong[iPtThr]++;    
                    if(n2ProngSignal[iPtThr] > 0 && n3ProngSignal[iPtThr] > 0)
                        nTriggeredSignalEvents2ProngAnd3Prong[iPtThr]++;            
                }
            }
        }
    }

    TGraphErrors* gFracTrigEvents2Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEvents2Prong->SetName("gFracTrigEvents2Prong");
    TGraphErrors* gFracTrigEvents3Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEvents3Prong->SetName("gFracTrigEvents3Prong");
    TGraphErrors* gFracTrigEventsDouble2Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEventsDouble2Prong->SetName("gFracTrigEventsDouble2Prong");
    TGraphErrors* gFracTrigEventsDouble3Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEventsDouble3Prong->SetName("gFracTrigEventsDouble3Prong");
    TGraphErrors* gFracTrigEvents2ProngAnd3Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEvents2ProngAnd3Prong->SetName("gFracTrigEvents2ProngAnd3Prong");
    TGraphErrors* gFracTrigEvents2ProngOr3Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEvents2ProngOr3Prong->SetName("gFracTrigEvents2ProngOr3Prong");

    TGraphErrors* gPurityTrigEvents2Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEvents2Prong->SetName("gPurityTrigEvents2Prong");
    TGraphErrors* gPurityTrigEvents3Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEvents3Prong->SetName("gPurityTrigEvents3Prong");
    TGraphErrors* gPurityTrigEventsDouble2Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEventsDouble2Prong->SetName("gPurityTrigEventsDouble2Prong");
    TGraphErrors* gPurityTrigEventsDouble3Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEventsDouble3Prong->SetName("gPurityTrigEventsDouble3Prong");
    TGraphErrors* gPurityTrigEvents2ProngAnd3Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEvents2ProngAnd3Prong->SetName("gPurityTrigEvents2ProngAnd3Prong");
    TGraphErrors* gPurityTrigEvents2ProngOr3Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEvents2ProngOr3Prong->SetName("gPurityTrigEvents2ProngOr3Prong");
    
    for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
    {
        gFracTrigEvents2Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEvents2Prong[iPtThr])/nEvents);
        gFracTrigEvents3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEvents3Prong[iPtThr])/nEvents);
        gFracTrigEventsDouble2Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEventsDouble2Prong[iPtThr])/nEvents);
        gFracTrigEventsDouble3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEventsDouble3Prong[iPtThr])/nEvents);
        gFracTrigEvents2ProngAnd3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEvents2ProngAnd3Prong[iPtThr])/nEvents);
        gFracTrigEvents2ProngOr3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEvents2ProngOr3Prong[iPtThr])/nEvents);

        gPurityTrigEvents2Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEvents2Prong[iPtThr])/static_cast<double>(nTriggeredEvents2Prong[iPtThr]));
        gPurityTrigEvents3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEvents3Prong[iPtThr])/static_cast<double>(nTriggeredEvents3Prong[iPtThr]));
        gPurityTrigEventsDouble2Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEventsDouble2Prong[iPtThr])/static_cast<double>(nTriggeredEventsDouble2Prong[iPtThr]));
        gPurityTrigEventsDouble3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEventsDouble3Prong[iPtThr])/static_cast<double>(nTriggeredEventsDouble3Prong[iPtThr]));
        gPurityTrigEvents2ProngAnd3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEvents2ProngAnd3Prong[iPtThr])/static_cast<double>(nTriggeredEvents2ProngAnd3Prong[iPtThr]));
        gPurityTrigEvents2ProngOr3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEvents2ProngOr3Prong[iPtThr])/static_cast<double>(nTriggeredEvents2ProngOr3Prong[iPtThr]));
    }

    TFile outfile(outFileName.data(), "recreate");
    gFracTrigEvents2Prong->Write();
    gFracTrigEvents3Prong->Write();
    gFracTrigEventsDouble2Prong->Write();
    gFracTrigEventsDouble3Prong->Write();
    gFracTrigEvents2ProngAnd3Prong->Write();
    gFracTrigEvents2ProngOr3Prong->Write();
    gPurityTrigEvents2Prong->Write();
    gPurityTrigEvents3Prong->Write();
    gPurityTrigEventsDouble2Prong->Write();
    gPurityTrigEventsDouble3Prong->Write();
    gPurityTrigEvents2ProngAnd3Prong->Write();
    gPurityTrigEvents2ProngOr3Prong->Write();
    outfile.Close();
}

//______________________________________________________________________________________________
void getTriggeredEventsWithB(TString cfgFileName="config.yml")
{

    //Load configs from yaml file
    YAML::Node config = YAML::LoadFile(cfgFileName.Data());
    if (config.IsNull()) {
        cerr << "Yaml config file not found! Exit" << endl;
        return;
    }

    string inFileName = config["infile"]["path"].as<string>();
    string inDirName = config["infile"]["dir"].as<string>();
    string outFileName = config["outfile"].as<string>();

    int nMaxEvents = config["nevents"].as<int>();

    vector<double> PtThresholds = config["ptthresholds"].as<vector<double> >(); // still use pT of the D
    const unsigned int nPtThresholds = PtThresholds.size();
    
    bool applyCuts = static_cast<bool>(config["applycuts"].as<int>());

    vector<double> PtLimsBeauty3Prong, PtLimsBeauty4Prong;
    vector<map<string, double> > mapOfCutsBeauty3Prongs = buildMapOfCuts(config["cutsBeauty3Prong"], PtLimsBeauty3Prong);
    vector<map<string, double> > mapOfCutsBeauty4Prongs = buildMapOfCuts(config["cutsBeauty4Prong"], PtLimsBeauty4Prong);

    //start analysis
    TFile* infile = TFile::Open(inFileName.data());
    TTree* evTree = dynamic_cast<TTree*>(infile->Get(Form("%s/fRecoTree", inDirName.data())));
    vector<Beauty3Prong> *beauty3ProngVec = nullptr;
    vector<Beauty4Prong> *beauty4ProngVec = nullptr;
    float zVtx = -999.;
    evTree->SetBranchAddress("Beauty3Prong", &beauty3ProngVec);
    evTree->SetBranchAddress("Beauty4Prong", &beauty4ProngVec);
    evTree->SetBranchAddress("zVtxReco", &zVtx);

    int nTriggeredEvents2Prong[nPtThresholds];
    int nTriggeredEvents3Prong[nPtThresholds];
    int nTriggeredEventsDouble2Prong[nPtThresholds];
    int nTriggeredEventsDouble3Prong[nPtThresholds];
    int nTriggeredEvents2ProngAnd3Prong[nPtThresholds];
    int nTriggeredEvents2ProngOr3Prong[nPtThresholds];

    int nTriggeredSignalEvents2Prong[nPtThresholds];
    int nTriggeredSignalEvents3Prong[nPtThresholds];
    int nTriggeredSignalEventsDouble2Prong[nPtThresholds];
    int nTriggeredSignalEventsDouble3Prong[nPtThresholds];
    int nTriggeredSignalEvents2ProngAnd3Prong[nPtThresholds];
    int nTriggeredSignalEvents2ProngOr3Prong[nPtThresholds];

    for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
    {
        nTriggeredEvents2Prong[iPtThr] = 0;
        nTriggeredEvents3Prong[iPtThr] = 0;
        nTriggeredEventsDouble2Prong[iPtThr] = 0;
        nTriggeredEventsDouble3Prong[iPtThr] = 0;
        nTriggeredEvents2ProngAnd3Prong[iPtThr] = 0;
        nTriggeredEvents2ProngOr3Prong[iPtThr] = 0;

        nTriggeredSignalEvents2Prong[iPtThr] = 0;
        nTriggeredSignalEvents3Prong[iPtThr] = 0;
        nTriggeredSignalEventsDouble2Prong[iPtThr] = 0;
        nTriggeredSignalEventsDouble3Prong[iPtThr] = 0;
        nTriggeredSignalEvents2ProngAnd3Prong[iPtThr] = 0;
        nTriggeredSignalEvents2ProngOr3Prong[iPtThr] = 0;
    }
    int nEvents = (nMaxEvents < evTree->GetEntriesFast()) ? nMaxEvents : evTree->GetEntriesFast();

    for(int iEv=0; iEv<nEvents; iEv++)
    {
        if((iEv+1)%100000 == 0)
            cout << "Processing event " << iEv+1 << endl;

        evTree->GetEvent(iEv);

        int n2Prong[nPtThresholds];
        int n3Prong[nPtThresholds];
        int n2ProngSignal[nPtThresholds];
        int n3ProngSignal[nPtThresholds];
        for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
        {
            n2Prong[iPtThr] = 0;
            n3Prong[iPtThr] = 0;
            n2ProngSignal[iPtThr] = 0;
            n3ProngSignal[iPtThr] = 0;
        }

        //loop over beauty 3 prongs (D0 + pion)
        for(unsigned int i2Prong=0; i2Prong<beauty3ProngVec->size(); i2Prong++)
        {
            float pt = beauty3ProngVec->at(i2Prong).fPtD0; // still used pT D
            int ptBin = findPtBin(pt, PtLimsBeauty3Prong);
            bool isSignalD = beauty3ProngVec->at(i2Prong).fCandType&BIT(7);
            bool isSignalB = beauty3ProngVec->at(i2Prong).fCandType&BIT(0);

            for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
            {
                if(pt > PtThresholds[iPtThr] && (!applyCuts || isBeauty3ProngSelected(beauty3ProngVec->at(i2Prong), mapOfCutsBeauty3Prongs[ptBin])))
                {
                    n2Prong[iPtThr]++;
                    if(isSignalD)
                        n2ProngSignal[iPtThr]++;
                }
            }
        }

        //loop over beauty 4 prongs (charm 3 prong + pion)
        for(unsigned int i3Prong=0; i3Prong<beauty4ProngVec->size(); i3Prong++)
        {
            float pt = beauty4ProngVec->at(i3Prong).fPtD; // still used pT D
            int ptBin = findPtBin(pt, PtLimsBeauty4Prong);
            bool isSignalD = beauty4ProngVec->at(i3Prong).fCandType&BIT(7);
            bool isSignalB = beauty4ProngVec->at(i3Prong).fCandType&BIT(0);

            for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
            {
                if(pt > PtThresholds[iPtThr] && (!applyCuts || isBeauty4ProngSelected(beauty4ProngVec->at(i3Prong), mapOfCutsBeauty4Prongs[ptBin])))
                {
                    n3Prong[iPtThr]++;
                    if(isSignalD)
                        n3ProngSignal[iPtThr]++;
                }
            }
        }

        for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
        {
            if(n2Prong[iPtThr] > 0)
            {
                nTriggeredEvents2Prong[iPtThr]++;
                if(n2ProngSignal[iPtThr] > 0)
                    nTriggeredSignalEvents2Prong[iPtThr]++;

                if(n2Prong[iPtThr] > 1)
                {
                    nTriggeredEventsDouble2Prong[iPtThr]++;
                    if(n2ProngSignal[iPtThr] > 1)
                        nTriggeredSignalEventsDouble2Prong[iPtThr]++;
                }
            }

            if(n3Prong[iPtThr] > 0)
            {

                nTriggeredEvents3Prong[iPtThr]++;
                 if(n3ProngSignal[iPtThr] > 0)
                    nTriggeredSignalEvents3Prong[iPtThr]++;

                if(n3Prong[iPtThr] > 1)
                {
                    nTriggeredEventsDouble3Prong[iPtThr]++;
                    if(n3ProngSignal[iPtThr] > 0)
                        nTriggeredSignalEvents3Prong[iPtThr]++;
                }
            }

            if(n2Prong[iPtThr] > 0 || n3Prong[iPtThr] > 0)
            {
                nTriggeredEvents2ProngOr3Prong[iPtThr]++;
                 if(n2ProngSignal[iPtThr] > 0 || n3ProngSignal[iPtThr] > 0)
                    nTriggeredSignalEvents2ProngOr3Prong[iPtThr]++;    
            }

            if(n2Prong[iPtThr] > 0 && n3Prong[iPtThr] > 0)
            {
                nTriggeredEvents2ProngAnd3Prong[iPtThr]++;    
                if(n2ProngSignal[iPtThr] > 0 && n3ProngSignal[iPtThr] > 0)
                    nTriggeredSignalEvents2ProngAnd3Prong[iPtThr]++;            
            }
        }
    }

    TGraphErrors* gFracTrigEvents2Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEvents2Prong->SetName("gFracTrigEvents2Prong");
    TGraphErrors* gFracTrigEvents3Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEvents3Prong->SetName("gFracTrigEvents3Prong");
    TGraphErrors* gFracTrigEventsDouble2Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEventsDouble2Prong->SetName("gFracTrigEventsDouble2Prong");
    TGraphErrors* gFracTrigEventsDouble3Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEventsDouble3Prong->SetName("gFracTrigEventsDouble3Prong");
    TGraphErrors* gFracTrigEvents2ProngAnd3Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEvents2ProngAnd3Prong->SetName("gFracTrigEvents2ProngAnd3Prong");
    TGraphErrors* gFracTrigEvents2ProngOr3Prong = new TGraphErrors(nPtThresholds);
    gFracTrigEvents2ProngOr3Prong->SetName("gFracTrigEvents2ProngOr3Prong");

    TGraphErrors* gPurityTrigEvents2Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEvents2Prong->SetName("gPurityTrigEvents2Prong");
    TGraphErrors* gPurityTrigEvents3Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEvents3Prong->SetName("gPurityTrigEvents3Prong");
    TGraphErrors* gPurityTrigEventsDouble2Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEventsDouble2Prong->SetName("gPurityTrigEventsDouble2Prong");
    TGraphErrors* gPurityTrigEventsDouble3Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEventsDouble3Prong->SetName("gPurityTrigEventsDouble3Prong");
    TGraphErrors* gPurityTrigEvents2ProngAnd3Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEvents2ProngAnd3Prong->SetName("gPurityTrigEvents2ProngAnd3Prong");
    TGraphErrors* gPurityTrigEvents2ProngOr3Prong = new TGraphErrors(nPtThresholds);
    gPurityTrigEvents2ProngOr3Prong->SetName("gPurityTrigEvents2ProngOr3Prong");
    
    for(unsigned int iPtThr=0; iPtThr<nPtThresholds; iPtThr++)
    {
        gFracTrigEvents2Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEvents2Prong[iPtThr])/nEvents);
        gFracTrigEvents3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEvents3Prong[iPtThr])/nEvents);
        gFracTrigEventsDouble2Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEventsDouble2Prong[iPtThr])/nEvents);
        gFracTrigEventsDouble3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEventsDouble3Prong[iPtThr])/nEvents);
        gFracTrigEvents2ProngAnd3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEvents2ProngAnd3Prong[iPtThr])/nEvents);
        gFracTrigEvents2ProngOr3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredEvents2ProngOr3Prong[iPtThr])/nEvents);

        gPurityTrigEvents2Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEvents2Prong[iPtThr])/static_cast<double>(nTriggeredEvents2Prong[iPtThr]));
        gPurityTrigEvents3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEvents3Prong[iPtThr])/static_cast<double>(nTriggeredEvents3Prong[iPtThr]));
        gPurityTrigEventsDouble2Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEventsDouble2Prong[iPtThr])/static_cast<double>(nTriggeredEventsDouble2Prong[iPtThr]));
        gPurityTrigEventsDouble3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEventsDouble3Prong[iPtThr])/static_cast<double>(nTriggeredEventsDouble3Prong[iPtThr]));
        gPurityTrigEvents2ProngAnd3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEvents2ProngAnd3Prong[iPtThr])/static_cast<double>(nTriggeredEvents2ProngAnd3Prong[iPtThr]));
        gPurityTrigEvents2ProngOr3Prong->SetPoint(iPtThr, PtThresholds[iPtThr], static_cast<double>(nTriggeredSignalEvents2ProngOr3Prong[iPtThr])/static_cast<double>(nTriggeredEvents2ProngOr3Prong[iPtThr]));
    }

    TFile outfile(outFileName.data(), "recreate");
    gFracTrigEvents2Prong->Write();
    gFracTrigEvents3Prong->Write();
    gFracTrigEventsDouble2Prong->Write();
    gFracTrigEventsDouble3Prong->Write();
    gFracTrigEvents2ProngAnd3Prong->Write();
    gFracTrigEvents2ProngOr3Prong->Write();
    gPurityTrigEvents2Prong->Write();
    gPurityTrigEvents3Prong->Write();
    gPurityTrigEventsDouble2Prong->Write();
    gPurityTrigEventsDouble3Prong->Write();
    gPurityTrigEvents2ProngAnd3Prong->Write();
    gPurityTrigEvents2ProngOr3Prong->Write();
    outfile.Close();
}
