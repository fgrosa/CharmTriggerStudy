#if !defined (__CINT__) || defined (__CLING__)

#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include <TChain.h>
#include <TGrid.h>

#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"

#include "AliPhysicsSelectionTask.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAnalysisTaskSECleanupVertexingHF.h"

#include "AliAnalysisTaskSECharmTriggerStudy.h"

#endif

using namespace std;

//______________________________________________
void RunAnalysisCharmTrigger(TString configfilename, TString runMode="full", bool mergeviajdl=true)
{
    //load config
    YAML::Node config = YAML::LoadFile(configfilename.Data());
    string aliPhysVersion = config["AliPhysicVersion"].as<string>();
    string gridDataDir = config["datadir"].as<string>();
    string gridDataPattern = config["datapattern"].as<string>();
    string gridWorkingDir = config["gridworkdir"].as<string>();

    string sSystem = config["system"].as<string>();
    int System=-1;
    if(sSystem=="pp")
        System = AliAnalysisTaskSECharmTriggerStudy::kpp;
    else if(sSystem=="PbPb")
        System = AliAnalysisTaskSECharmTriggerStudy::kPbPb;
    else
    {
        cerr << "ERROR: Only pp and PbPb are supported, please check your config file. Exit" << endl;
        return;
    }

    vector<int> runlist = config["runs"].as<vector<int>>();
    const int nruns = runlist.size();
    bool isRunOnMC = strstr(gridDataDir.c_str(),"sim");
    if(!isRunOnMC)
    {
        cerr << "ERROR: AliAnalysisTaskSECharmTriggerStudy can only be run on a MC production!" << endl;
        return;
    }

    bool local = false;
    bool gridTest = false;
    string pathToLocalAODfiles = "";
    if(config["runtype"].as<string>() == "local")
    {
        local = true;
        pathToLocalAODfiles = config["pathtolocalAOD"].as<string>();
    }
    else
    {
        TGrid::Connect("alien://");
        if(config["runtype"].as<string>() == "test")
            gridTest = true;
    }

    // since we will compile a class, tell root where to look for headers
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    AliPhysicsSelectionTask *physseltask = reinterpret_cast<AliPhysicsSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s (%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"),isRunOnMC)));
    AliAnalysisTaskPIDResponse *pidResp = reinterpret_cast<AliAnalysisTaskPIDResponse *>(gInterpreter->ProcessLine(Form(".x %s (%d)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"),isRunOnMC)));

    gInterpreter->ProcessLine(".L AliAnalysisTaskSECharmTriggerStudy.cxx+g");
    AliAnalysisTaskSECharmTriggerStudy *tasktrigger = reinterpret_cast<AliAnalysisTaskSECharmTriggerStudy*>(gInterpreter->ProcessLine(Form(".x %s(%d,%d,%d,%d,%d,\"%s\")", gSystem->ExpandPathName("AddTaskCharmTriggerStudy.C"), System, true, true, true, true, "Test")));

    if(System==AliAnalysisTaskSECharmTriggerStudy::kPbPb) {
        AliAnalysisTaskSECleanupVertexingHF *taskclean =reinterpret_cast<AliAnalysisTaskSECleanupVertexingHF *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskCleanupVertexingHF.C"))));
    }

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chainAOD = new TChain("aodTree");
        TChain *chainAODfriend = new TChain("aodTree");

        // add a few files to the chain (change this so that your local files are added)
        chainAOD->Add(Form("%s/AliAOD.root",pathToLocalAODfiles.data()));
        chainAODfriend->Add(Form("%s/AliAOD.VertexingHF.root",pathToLocalAODfiles.data()));

        chainAOD->AddFriend(chainAODfriend);


        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chainAOD);

    }
    else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();

        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

        //make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskSECharmTriggerStudy.h AliAnalysisTaskSECharmTriggerStudy.cxx");
        alienHandler->SetAnalysisSource("AliAnalysisTaskSECharmTriggerStudy.cxx");

        // select the aliphysics version.
        alienHandler->SetAliPhysicsVersion(aliPhysVersion.data());

        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");

        // select the input data
        alienHandler->SetGridDataDir(gridDataDir.data());
        alienHandler->SetDataPattern(Form("%s/*AliAOD.root",gridDataPattern.data()));
        alienHandler->SetFriendChainName("AliAOD.VertexingHF.root");

        // MC has no prefix, data has prefix 000
        if(!isRunOnMC)
            alienHandler->SetRunPrefix("000");

        for(int k=0; k<nruns; k++)
            alienHandler->AddRunNumber(runlist[k]);
        alienHandler->SetNrunsPerMaster(nruns);

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(1);
        alienHandler->SetExecutable("myTask.sh");

        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(true);
        alienHandler->SetKeepLogs(true);

        // merging: run with true to merge on grid
        // after re-running the jobs in SetRunMode("terminate")
        // (see below) mode, set SetMergeViaJDL(false)
        // to collect final results
        alienHandler->SetMaxMergeStages(3); //2, 3
        alienHandler->SetMergeViaJDL(mergeviajdl);

        // define the output folders
        alienHandler->SetGridWorkingDir(gridWorkingDir.data());
        alienHandler->SetGridOutputDir("output");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);

        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        }
        else {
            // else launch the full grid analysis
            alienHandler->SetRunMode(runMode.Data()); //terminate
            mgr->StartAnalysis("grid");
        }
    }
}
