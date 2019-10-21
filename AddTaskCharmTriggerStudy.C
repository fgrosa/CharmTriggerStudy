AliAnalysisTaskSECharmTriggerStudy *AddTaskCharmTriggerStudy(int system = AliAnalysisTaskSECharmTriggerStudy::kpp,
                                                             bool enable2prongs = true,
                                                             bool enable3prongs = true,
                                                             bool enableDstars = false,
                                                             bool enableCascades = false,
                                                             bool fillOnlySignal = false,
                                                             TString suffix = "")
{
    //
    // Test macro for the AliAnalysisTaskSE for charm-trigger studies

    // Get the pointer to the existing analysis manager via the static access method.
    //============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AliAnalysisTaskSECharmTriggerStudy", "No analysis manager to connect to.");
        return NULL;
    }
    if (!mgr->GetInputEventHandler())
    {
        ::Error("AliAnalysisTaskSECharmTriggerStudy", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if (type.Contains("ESD"))
    {
        ::Error("AliAnalysisTaskSECharmTriggerStudy", "This task requires to run on AOD");
        return NULL;
    }

    // Analysis task
    AliAnalysisTaskSECharmTriggerStudy *chTask = new AliAnalysisTaskSECharmTriggerStudy("CharmTriggerStudyTask");
    chTask->Enable2Prongs(enable2prongs);
    chTask->Enable3Prongs(enable3prongs);
    chTask->EnableDstars(enableDstars);
    chTask->EnableCascades(enableCascades);
    chTask->SetSystem(system);
    mgr->AddTask(chTask);

    // Create containers for input/output
    TString contname = Form("cinputChTrigger%s", suffix.Data());
    AliAnalysisDataContainer *cinputcont = mgr->CreateContainer(contname.Data(), TChain::Class(), AliAnalysisManager::kInputContainer);

    TString outputfile = AliAnalysisManager::GetCommonFileName();
    TString outputdirname = Form("%s:PWGHF_D2H_CharmTrigger_%s", outputfile.Data(), suffix.Data());

    contname = Form("coutputChTrigger%s", suffix.Data());
    AliAnalysisDataContainer *coutputlist = mgr->CreateContainer(contname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outputdirname.Data());

    contname = Form("coutputChTriggerRecoTree%s", suffix.Data());
    AliAnalysisDataContainer *coutputrecotree = mgr->CreateContainer(contname.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer, outputdirname.Data());

    contname = Form("coutputChTriggerGenTree%s", suffix.Data());
    AliAnalysisDataContainer *coutputgentree = mgr->CreateContainer(contname.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer, outputdirname.Data());

    mgr->ConnectInput(chTask, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(chTask, 1, coutputlist);
    mgr->ConnectOutput(chTask, 2, coutputrecotree);
    mgr->ConnectOutput(chTask, 3, coutputgentree);

    return chTask;
}