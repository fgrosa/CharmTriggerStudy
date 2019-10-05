// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliAnalysisTaskSECharmTriggerStudy
// \brief task that produces an output tree for the charm trigger studies
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
// M. Puccio, maximiliano.puccio@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TList.h>
#include <TString.h>
#include <TH1F.h>
#include <TDatabasePDG.h>

#include "AliAnalysisTaskSECharmTriggerStudy.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliRDHFCuts.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSECharmTriggerStudy);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSECharmTriggerStudy::AliAnalysisTaskSECharmTriggerStudy() : AliAnalysisTaskSE(),
    fOutput(nullptr),
    fHistNEvents(nullptr),
    fRecoTree(nullptr),
    fGenTree(nullptr),
    fEventCuts{},
    fSystem(kpp),
    fAOD(nullptr),
    fAODProtection(1),
    fMCArray(nullptr),
    fRecoZvtx(-999.),
    fGenZvtx(-999.),
    fCharm2Prong{},
    fCharm3Prong{},
    fDstar{},
    fCharmCascade{},
    fGenCharmHadron{},
    fEnable2Prongs(true),
    fEnable3Prongs(true),
    fEnableDstars(false),
    fEnableCascades(false)
{
    /// Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSECharmTriggerStudy::AliAnalysisTaskSECharmTriggerStudy(const char *name) : AliAnalysisTaskSE(name),
    fOutput(nullptr),
    fHistNEvents(nullptr),
    fRecoTree(nullptr),
    fGenTree(nullptr),
    fEventCuts{},
    fSystem(kpp),
    fAOD(nullptr),
    fAODProtection(1),
    fMCArray(nullptr),
    fRecoZvtx(-999.),
    fGenZvtx(-999.),
    fCharm2Prong{},
    fCharm3Prong{},
    fDstar{},
    fCharmCascade{},
    fGenCharmHadron{},
    fEnable2Prongs(true),
    fEnable3Prongs(true),
    fEnableDstars(false),
    fEnableCascades(false)
{
    /// Default constructor

    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskSECharmTriggerStudy::~AliAnalysisTaskSECharmTriggerStudy()
{
    // Destructor
    if(fOutput) delete fOutput;
    if(fRecoTree) delete fRecoTree;
    if(fGenTree) delete fGenTree;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::UserCreateOutputObjects()
{
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");
    fEventCuts.AddQAplotsToList(fOutput);

    fHistNEvents = new TH1F("hNEvents", "number of events ", 13, 0.5, 12.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1, "nEventsRead");
    fHistNEvents->GetXaxis()->SetBinLabel(2, "nEvents Matched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(3, "nEvents Mismatched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(4, "nEventsAnal");
    fHistNEvents->GetXaxis()->SetBinLabel(5, "n. rejected due to not reco vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(6, "n. passing IsEvSelected");
    fHistNEvents->GetXaxis()->SetBinLabel(7, "no. of 2 prong candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(8, "no. of 3 prong candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(9, "no. of Dstar candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(10, "no. of cascade candidates");

    fHistNEvents->GetXaxis()->SetNdivisions(1, false);
    fHistNEvents->SetMinimum(0);

    fOutput->Add(fHistNEvents);

    fRecoTree = new TTree("fRecoTree", "Reconstructed charm hadron candidates");
    fRecoTree->Branch("zVtxReco", &fRecoZvtx);
    if(fEnable2Prongs)
        fRecoTree->Branch("Charm2Prong", &fCharm2Prong);
    if(fEnable3Prongs)
        fRecoTree->Branch("Charm3Prong", &fCharm3Prong);
    if(fEnableDstars)
        fRecoTree->Branch("Dstar", &fDstar);
    if(fEnableCascades)
        fRecoTree->Branch("CharmCascade", &fCharmCascade);

    fGenTree = new TTree("fGenTree", "Generate charm hadrons");
    fGenTree->Branch("zVtxGen", &fGenZvtx);
    fGenTree->Branch("GenCharmHadron", &fGenCharmHadron);

    PostData(1, fOutput);
    PostData(2, fRecoTree);
    PostData(3, fGenTree);

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::UserExec(Option_t * /*option*/)
{

    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());

    fHistNEvents->Fill(1); // all events
    if (fAODProtection >= 0)
    {
        int matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
        if (matchingAODdeltaAODlevel < 0 || (matchingAODdeltaAODlevel == 0 && fAODProtection == 1))
        {
            fHistNEvents->Fill(3);
            PostData(1, fOutput);
            return;
        }
        fHistNEvents->Fill(2);
    }

    TClonesArray *array3Prong = nullptr, *array2Prong = nullptr, *arrayCasc = nullptr, *arrayDstar = nullptr;

    if (!fAOD && AODEvent() && IsStandardAOD())
    {
        fAOD = dynamic_cast<AliAODEvent *>(AODEvent());
        AliAODHandler *aodHandler = dynamic_cast<AliAODHandler *>((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
        if (aodHandler->GetExtensions())
        {
            AliAODExtension *ext = dynamic_cast<AliAODExtension *>(aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root"));
            AliAODEvent *aodFromExt = ext->GetAOD();
            array2Prong = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("D0toKpi"));
            array3Prong = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Charm3Prong"));
            arrayDstar = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("Dstar"));
            arrayCasc = dynamic_cast<TClonesArray *>(aodFromExt->GetList()->FindObject("CascadesHF"));
        }
    }
    else if (fAOD)
    {
        array2Prong = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("D0toKpi"));
        array3Prong = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Charm3Prong"));
        arrayDstar = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("Dstar"));
        arrayCasc = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject("CascadesHF"));
    }

    if (!fAOD || (fEnable3Prongs && !array3Prong) || (fEnable2Prongs && !array2Prong) || (fEnableCascades && !arrayCasc) || (fEnableDstars && !arrayDstar))
    {
        AliWarning("Candidate branch not found!");
        return;
    }

    if (!fAOD->GetPrimaryVertex() || TMath::Abs(fAOD->GetMagneticField()) < 0.001)
        return;

    fHistNEvents->Fill(4); // count event

    fEventCuts.AcceptEvent(fAOD); // do not return yet (no physics selection applied for upgrade MC)

    if(!fEventCuts.PassedCut(AliEventCuts::kVertex) || !fEventCuts.PassedCut(AliEventCuts::kVertexPosition) || !fEventCuts.PassedCut(AliEventCuts::kVertexQuality))
    {
        fHistNEvents->Fill(5); // rejected for primary vtx
        PostData(1, fOutput);
        return;
    }

    fHistNEvents->Fill(6); // selected event

    fMCArray = dynamic_cast<TClonesArray*>(fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
    if(!fMCArray)
    {
      AliWarning("MC particles branch not found!");
      return;
    }

    // load MC header
    AliAODMCHeader* mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if(!mcHeader)
    {
      AliWarning("MC header branch not found!");
      return;
    }

    fGenZvtx = mcHeader->GetVtxZ();

    AliAODVertex *primVtx = dynamic_cast<AliAODVertex*>(fAOD->GetPrimaryVertex());
    fRecoZvtx = primVtx->GetZ();

    //loop on generated particles
    for (int iPart = 0; iPart < fMCArray->GetEntriesFast(); iPart++)
    {
        AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(fMCArray->UncheckedAt(iPart));
        if(!part)
            continue;

        int pdgCode = TMath::Abs(part->GetPdgCode());
        if(pdgCode == 411 || pdgCode == 421 || pdgCode == 431 || pdgCode == 413 || pdgCode == 4122)
        {
            int origin = AliVertexingHFUtils::CheckOrigin(fMCArray, part, true);
            if(origin != 4 && origin != 5)
                continue; //keep only prompt or feed-down

            int labDau[3] = {-1, -1, -1};
            int pdgCodeDau0 = -1;
            int decay = -1;
            if(pdgCode == 421 && fEnable2Prongs) //Dzero
            {
                if(part->GetNDaughters()==2)
                    decay = AliVertexingHFUtils::CheckD0Decay(fMCArray, part, labDau);
                if(decay!=1 || labDau[0]<0 || labDau[1]<0)
                    continue;
                pdgCodeDau0 = (dynamic_cast<AliAODMCParticle*>(fMCArray->UncheckedAt(labDau[0])))->GetPdgCode();
                if(TMath::Abs(pdgCodeDau0) == 321)
                    fGenCharmHadron.push_back(FillCharmGen(part, origin, kDzerotoKpi));
                else
                    fGenCharmHadron.push_back(FillCharmGen(part, origin, kDzerotopiK));
            }
            else if(pdgCode == 411 && fEnable3Prongs) //Dplus
            {
                decay = AliVertexingHFUtils::CheckDplusDecay(fMCArray, part, labDau);
                if(decay>=1 && labDau[0]>=0 && labDau[1]>=0) {
                    fGenCharmHadron.push_back(FillCharmGen(part, origin, kDplustoKpipi));
                    continue;
                }
                decay = AliVertexingHFUtils::CheckDplusKKpiDecay(fMCArray, part, labDau);
                if(decay!=1 && labDau[0]<0 && labDau[1]<0)
                    continue;

                pdgCodeDau0 = (dynamic_cast<AliAODMCParticle*>(fMCArray->UncheckedAt(labDau[0])))->GetPdgCode();
                if(TMath::Abs(pdgCodeDau0) == 321)
                    fGenCharmHadron.push_back(FillCharmGen(part, origin, kDplustoKKpi));
                else
                    fGenCharmHadron.push_back(FillCharmGen(part, origin, kDplustopiKK));
            }
            else if(pdgCode == 431 && fEnable3Prongs) //Ds
            {
                decay = AliVertexingHFUtils::CheckDsDecay(fMCArray, part, labDau);
                if(decay!=1 || labDau[0]<0 || labDau[1]<0) //keep only Ds -> phipi --> KKpi (to be discussed)
                    continue;
                pdgCodeDau0 = (dynamic_cast<AliAODMCParticle*>(fMCArray->UncheckedAt(labDau[0])))->GetPdgCode();
                if(TMath::Abs(pdgCodeDau0) == 321)
                    fGenCharmHadron.push_back(FillCharmGen(part, origin, kDstoKKpi));
                else
                    fGenCharmHadron.push_back(FillCharmGen(part, origin, kDstopiKK));
            }
            else if(pdgCode == 413 && fEnableDstars) //Dstar
            {
                decay = AliVertexingHFUtils::CheckDstarDecay(fMCArray, part, labDau);
                if(decay!=1 || labDau[0]<0 || labDau[1]<0)
                    continue;
                fGenCharmHadron.push_back(FillCharmGen(part, origin, kDstartoKpipi));
            }
            else if(pdgCode == 4122) //Lc
            {
                if(fEnable3Prongs)
                {
                    decay = AliVertexingHFUtils::CheckLcpKpiDecay(fMCArray, part, labDau);
                    if(decay>=1 && labDau[0]>=0 && labDau[1]>=0) {
                        pdgCodeDau0 = (dynamic_cast<AliAODMCParticle*>(fMCArray->UncheckedAt(labDau[0])))->GetPdgCode();
                        if(TMath::Abs(pdgCodeDau0) == 2212)
                            fGenCharmHadron.push_back(FillCharmGen(part, origin, kLctopKpi));
                        else
                            fGenCharmHadron.push_back(FillCharmGen(part, origin, kLctopiKp));
                        continue;
                    }
                }
                if(fEnableCascades)
                {
                    decay = AliVertexingHFUtils::CheckLcV0bachelorDecay(fMCArray, part, labDau);
                    if(labDau[0]>=0 && labDau[1]>=0) {
                        if(decay == 1)
                            fGenCharmHadron.push_back(FillCharmGen(part, origin, kLctopiLambda));
                        else if(decay == 2)
                            fGenCharmHadron.push_back(FillCharmGen(part, origin, kLctopK0s));
                    }
                }
            }
        }
    }

    //loop on 2 prongs
    if(fEnable2Prongs)
    {
        for (int i2Prong = 0; i2Prong < array2Prong->GetEntriesFast(); i2Prong++)
        {
            AliAODRecoDecayHF2Prong *d = dynamic_cast<AliAODRecoDecayHF2Prong *>(array2Prong->UncheckedAt(i2Prong));
            if (!d || (d->GetSelectionMap() && !d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)))
                continue;
            fHistNEvents->Fill(9);

            //check if primary vtx is set
            bool unsetvtx = false;
            if(!d->GetOwnPrimaryVtx())
            {
                if(!d->GetOwnPrimaryVtx())
                {
                    d->SetOwnPrimaryVtx(primVtx);
                    unsetvtx = true;
                }
            }

            //if pp recompute primary vertex without daughters
            bool isvtxrecalc = false;
            AliAODVertex* origownvtx = nullptr;
            if(fSystem == kpp) {
                origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
                isvtxrecalc = RecalcOwnPrimaryVertex(d);
                if(!isvtxrecalc)
                    CleanOwnPrimaryVertex(d, origownvtx);
            }

            //fill vector of 2prongs
            fCharm2Prong.push_back(FillCharm2Prong(d));

            if(isvtxrecalc)
                CleanOwnPrimaryVertex(d, origownvtx);
            if(unsetvtx)
                d->UnsetOwnPrimaryVtx();
        }
    }

    //loop on 3 prongs
    if(fEnable3Prongs)
    {
        for (int i3Prong = 0; i3Prong < array3Prong->GetEntriesFast(); i3Prong++)
        {
            AliAODRecoDecayHF3Prong *d = dynamic_cast<AliAODRecoDecayHF3Prong *>(array3Prong->UncheckedAt(i3Prong));
            if (!d || (d->GetSelectionMap() && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts) && !d->HasSelectionBit(AliRDHFCuts::kDsCuts) && !d->HasSelectionBit(AliRDHFCuts::kLcCuts)))
                continue;

            fHistNEvents->Fill(10);

            //check if primary vtx is set
            bool unsetvtx = false;
            if(!d->GetOwnPrimaryVtx())
            {
                if(!d->GetOwnPrimaryVtx())
                {
                    d->SetOwnPrimaryVtx(primVtx);
                    unsetvtx = true;
                }
            }

            //if pp recompute primary vertex without daughters
            bool isvtxrecalc = false;
            AliAODVertex* origownvtx = nullptr;
            if(fSystem == kpp) {
                origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
                isvtxrecalc = RecalcOwnPrimaryVertex(d);
                if(!isvtxrecalc)
                    CleanOwnPrimaryVertex(d, origownvtx);
            }

            //fill vector of 3prongs
            fCharm3Prong.push_back(FillCharm3Prong(d));

            if(isvtxrecalc)
                CleanOwnPrimaryVertex(d, origownvtx);
            if(unsetvtx)
                d->UnsetOwnPrimaryVtx();
        }
    }

    //loop on Dstars
    if(fEnableDstars)
    {
        for (int iDstar = 0; iDstar < arrayDstar->GetEntriesFast(); iDstar++)
        {
            AliAODRecoCascadeHF *d = dynamic_cast<AliAODRecoCascadeHF *>(arrayDstar->UncheckedAt(iDstar));
            if (!d || (d->GetSelectionMap() && !d->HasSelectionBit(AliRDHFCuts::kDstarCuts)))
                continue;
            AliAODRecoDecayHF2Prong *d0 = d->Get2Prong();
            if (!d0)
                continue;

            fHistNEvents->Fill(11);

            //check if primary vtx is set
            bool unsetvtx = false;
            if(!d->GetOwnPrimaryVtx())
            {
                if(!d->GetOwnPrimaryVtx())
                {
                    d->SetOwnPrimaryVtx(primVtx);
                    unsetvtx = true;
                }
            }

            //if pp recompute primary vertex without daughters
            bool isvtxrecalc = false;
            AliAODVertex* origownvtx = nullptr;
            if(fSystem == kpp) {
                origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
                isvtxrecalc = RecalcOwnPrimaryVertex(d);
                if(!isvtxrecalc)
                    CleanOwnPrimaryVertex(d, origownvtx);
            }

            //fill vector of dstars
            fDstar.push_back(FillDstar(d, d0));

            if(isvtxrecalc)
                CleanOwnPrimaryVertex(d, origownvtx);
            if(unsetvtx)
                d->UnsetOwnPrimaryVtx();
        }
    }

    //loop on cascades
    if(fEnableCascades)
    {
        for (int iCasc = 0; iCasc < arrayCasc->GetEntriesFast(); iCasc++)
        {
            AliAODRecoCascadeHF *lc = dynamic_cast<AliAODRecoCascadeHF *>(arrayCasc->UncheckedAt(iCasc));
            if (!lc)
                continue;
            AliAODv0 *v0part = lc->Getv0();
            if (!v0part)
                continue;
            fHistNEvents->Fill(12);

            //check if primary vtx is set
            bool unsetvtx = false;
            if(!lc->GetOwnPrimaryVtx())
            {
                if(!lc->GetOwnPrimaryVtx())
                {
                    lc->SetOwnPrimaryVtx(primVtx);
                    unsetvtx = true;
                }
            }

            //if pp recompute primary vertex without daughters
            bool isvtxrecalc = false;
            AliAODVertex* origownvtx = nullptr;
            if(fSystem == kpp) {
                origownvtx = new AliAODVertex(*lc->GetOwnPrimaryVtx());
                isvtxrecalc = RecalcOwnPrimaryVertex(lc);
                if(!isvtxrecalc)
                    CleanOwnPrimaryVertex(lc, origownvtx);
            }

            //fill vector of cascades
            fCharmCascade.push_back(FillCharmCascade(lc, v0part));

            if(isvtxrecalc)
                CleanOwnPrimaryVertex(lc, origownvtx);
            if(unsetvtx)
                lc->UnsetOwnPrimaryVtx();
        }
    }

    fRecoTree->Fill();
    fGenTree->Fill();

    PostData(1, fOutput);
    PostData(2, fRecoTree);
    PostData(3, fGenTree);

    return;
}

//________________________________________________________________________
Charm2Prong AliAnalysisTaskSECharmTriggerStudy::FillCharm2Prong(AliAODRecoDecayHF2Prong* cand)
{
    Charm2Prong ch2Prong;
    ch2Prong.fPt = cand->Pt();
    ch2Prong.fY = cand->Y(421);
    ch2Prong.fInvMassD0 = cand->InvMassD0();
    ch2Prong.fInvMassD0bar = cand->InvMassD0bar();
    ch2Prong.fCosP = cand->CosPointingAngle();
    ch2Prong.fCosPXY = cand->CosPointingAngleXY();
    ch2Prong.fDecayLength = cand->DecayLength();
    ch2Prong.fDecayLengthXY = cand->DecayLengthXY();
    ch2Prong.fNormDecayLength = cand->NormalizedDecayLength();
    ch2Prong.fNormDecayLengthXY = cand->NormalizedDecayLengthXY();
    ch2Prong.fImpParProd = cand->Getd0Prong(0)*cand->Getd0Prong(1);

    int pdgDgD0toKpi[2] = {321, 211};
    ch2Prong.fGenLabel = cand->MatchToMC(421, fMCArray, 2, pdgDgD0toKpi);
    ch2Prong.fDecay = kNone;
    ch2Prong.fCandType = 0;
    int origin = -1;
    if(ch2Prong.fGenLabel >= 0){
        AliAODMCParticle* partD0 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(ch2Prong.fGenLabel));
        origin = AliVertexingHFUtils::CheckOrigin(fMCArray, partD0, true);
        if(origin == 4)
        {
            ch2Prong.fCandType |= kSignal;
            ch2Prong.fCandType |= kPrompt;
        }
        else if(origin == 5)
        {
            ch2Prong.fCandType |= kSignal;
            ch2Prong.fCandType |= kFeedDown;
        }
        else
        {
            ch2Prong.fCandType |= kSignal; // no prompt, no feed-down --> weird stuff
        }

        int labDau0 = dynamic_cast<AliAODTrack*>(cand->GetDaughter(0))->GetLabel();
        AliAODMCParticle* dauPart0 = dynamic_cast<AliAODMCParticle*>(fMCArray->UncheckedAt(TMath::Abs(labDau0)));
        int pdgCode0 = TMath::Abs(dauPart0->GetPdgCode());
        if(pdgCode0 == 321)
            ch2Prong.fDecay = kDzerotoKpi;
        else if(pdgCode0 == 211)
            ch2Prong.fDecay = kDzerotopiK;
        else
            ch2Prong.fDecay = kNone;
    }
    else {
        ch2Prong.fCandType |= kBackground;
    }

    return ch2Prong;
}

//________________________________________________________________________
Charm3Prong AliAnalysisTaskSECharmTriggerStudy::FillCharm3Prong(AliAODRecoDecayHF3Prong* cand)
{
    Charm3Prong ch3Prong;
    ch3Prong.fPt = cand->Pt();
    ch3Prong.fInvMassDplus = cand->InvMassDplus();
    ch3Prong.fInvMassDstoKKpi = cand->InvMassDsKKpi();
    ch3Prong.fInvMassDstopiKK = cand->InvMassDspiKK();
    ch3Prong.fInvMassLctopKpi = cand->InvMassLcpKpi();
    ch3Prong.fInvMassLctopiKp = cand->InvMassLcpiKp();
    ch3Prong.fYDplus = cand->YDplus();
    ch3Prong.fYDs = cand->YDs();
    ch3Prong.fYLc = cand->YLc();
    ch3Prong.fInvMassPhiKKpi = cand->InvMass2Prongs(0, 1, 321, 321);
    ch3Prong.fInvMassPhipiKK = cand->InvMass2Prongs(1, 2, 321, 321);
    ch3Prong.fCosP = cand->CosPointingAngle();
    ch3Prong.fCosPXY = cand->CosPointingAngleXY();
    ch3Prong.fDecayLength = cand->DecayLength();
    ch3Prong.fDecayLengthXY = cand->DecayLengthXY();
    ch3Prong.fNormDecayLength = cand->NormalizedDecayLength();
    ch3Prong.fNormDecayLengthXY = cand->NormalizedDecayLengthXY();
    ch3Prong.fSigmaVtx = cand->GetSigmaVert();

    int pdgDgDplustoKpipi[3] = {321, 211, 211};
    int pdgDgDstoKKpi[3] = {321, 321, 211};
    int pdgDgLctopKpi[3] = {2122, 321, 211};
    int origin = -1;
    ch3Prong.fGenLabel = -1;
    ch3Prong.fDecay = kNone;

    AliAODMCParticle* part3prong = dynamic_cast<AliAODMCParticle*>(fMCArray->At(ch3Prong.fGenLabel));
    int labDplus = cand->MatchToMC(411, fMCArray, 3, pdgDgDplustoKpipi);
    int labDs = -1;
    int labDplustoKKpi = -1;
    int labLc = -1;
    if (labDplus < 0)
    {
        labDs = cand->MatchToMC(431, fMCArray, 3, pdgDgDstoKKpi);
        if (labDs < 0)
        {
            labDplustoKKpi = cand->MatchToMC(411, fMCArray, 3, pdgDgDstoKKpi);
            if (labDplustoKKpi < 0)
            {
                labLc = cand->MatchToMC(4122, fMCArray, 3, pdgDgDstoKKpi);
                if (labLc >= 0)
                {
                    ch3Prong.fGenLabel = labLc;

                    int labDau0 = dynamic_cast<AliAODTrack*>(cand->GetDaughter(0))->GetLabel();
                    AliAODMCParticle* dauPart0 = dynamic_cast<AliAODMCParticle*>(fMCArray->UncheckedAt(TMath::Abs(labDau0)));
                    int pdgCode0 = TMath::Abs(dauPart0->GetPdgCode());
                    if(pdgCode0 == 2122)
                        ch3Prong.fDecay = kLctopKpi;
                    else if(pdgCode0 == 211)
                        ch3Prong.fDecay = kLctopiKp;
                }
            }
            else {
                ch3Prong.fGenLabel = labDplustoKKpi;
                int labDau0 = dynamic_cast<AliAODTrack*>(cand->GetDaughter(0))->GetLabel();
                AliAODMCParticle* dauPart0 = dynamic_cast<AliAODMCParticle*>(fMCArray->UncheckedAt(TMath::Abs(labDau0)));
                int pdgCode0 = TMath::Abs(dauPart0->GetPdgCode());
                if(pdgCode0 == 321)
                    ch3Prong.fDecay = kDplustoKKpi;
                else if(pdgCode0 == 211)
                    ch3Prong.fDecay = kDplustopiKK;
            }
        }
        else
        {
            ch3Prong.fGenLabel = labDs;
            int labDau0 = dynamic_cast<AliAODTrack*>(cand->GetDaughter(0))->GetLabel();
            AliAODMCParticle* dauPart0 = dynamic_cast<AliAODMCParticle*>(fMCArray->UncheckedAt(TMath::Abs(labDau0)));
            int pdgCode0 = TMath::Abs(dauPart0->GetPdgCode());
            if(pdgCode0 == 321)
                ch3Prong.fDecay = kDstoKKpi;
            else if(pdgCode0 == 211)
                ch3Prong.fDecay = kDstopiKK;
        }
    }
    else
    {
        ch3Prong.fGenLabel = labDplus;
        ch3Prong.fDecay = kDplustoKpipi;
    }

    if(ch3Prong.fGenLabel >= 0){
        origin = AliVertexingHFUtils::CheckOrigin(fMCArray, part3prong, true);
        if(origin == 4)
        {
            ch3Prong.fCandType |= kSignal;
            ch3Prong.fCandType |= kPrompt;
        }
        else if(origin == 5)
        {
            ch3Prong.fCandType |= kSignal;
            ch3Prong.fCandType |= kFeedDown;
        }
        else
        {
            ch3Prong.fCandType |= kSignal; // no prompt, no feed-down --> weird stuff
        }
    }
    else {
        ch3Prong.fCandType |= kBackground;
    }

    return ch3Prong;
}

//________________________________________________________________________
Dstar AliAnalysisTaskSECharmTriggerStudy::FillDstar(AliAODRecoCascadeHF* cand, AliAODRecoDecayHF2Prong* dau)
{
    Dstar dstar;
    dstar.fInvMass = cand->InvMassDstarKpipi();
    dstar.fInvMassD0 = cand->InvMassD0();
    dstar.fPt = cand->Pt();
    dstar.fY = cand->Y(431);
    dstar.fCosPD0 = dau->CosPointingAngle();
    dstar.fCosPXYD0 = dau->CosPointingAngleXY();
    dstar.fDecayLengthD0 = dau->DecayLength();
    dstar.fDecayLengthXYD0 = dau->DecayLengthXY();
    dstar.fNormDecayLengthD0 = dau->NormalizedDecayLength();
    dstar.fNormDecayLengthXYD0 = dau->NormalizedDecayLengthXY();

    int pdgDgDstartoKpipi[3] = {321, 211, 211};
    int pdgDgD0toKpi[2] = {321, 211};
    dstar.fGenLabel = cand->MatchToMC(413, 421, pdgDgDstartoKpipi, pdgDgD0toKpi, fMCArray);
    dstar.fDecay = kNone;
    dstar.fCandType = 0;
    int origin = -1;
    if(dstar.fGenLabel >= 0){
        AliAODMCParticle* partD0 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(dstar.fGenLabel));
        origin = AliVertexingHFUtils::CheckOrigin(fMCArray, partD0, true);
        if(origin == 4)
        {
            dstar.fCandType |= kSignal;
            dstar.fCandType |= kPrompt;
        }
        else if(origin == 5)
        {
            dstar.fCandType |= kSignal;
            dstar.fCandType |= kFeedDown;
        }
        else
        {
            dstar.fCandType |= kSignal; // no prompt, no feed-down --> weird stuff
        }
        dstar.fDecay = kDstartoKpipi;
    }
    else {
        dstar.fCandType |= kBackground;
    }

    return dstar;
}

//________________________________________________________________________
CharmCascade AliAnalysisTaskSECharmTriggerStudy::FillCharmCascade(AliAODRecoCascadeHF* cand, AliAODv0* dau)
{
    CharmCascade chCasc;

    chCasc.fInvMassLctopK0s = cand->InvMassLctoK0sP();
    chCasc.fInvMassLctopiLambda = cand->InvMassLctoLambdaPi();
    chCasc.fInvMassK0s = dau->MassK0Short();
    if (cand->Charge() > 0)
        chCasc.fInvMassLambda = dau->MassLambda();
    else
        chCasc.fInvMassLambda = dau->MassAntiLambda();
    chCasc.fPt = cand->Pt();
    chCasc.fY = cand->Y(4122);
    chCasc.fCosPV0 = cand->CosV0PointingAngle();
    chCasc.fCosPXYV0 = cand->CosV0PointingAngleXY();

    int pdgDgLctopK0s[3] = {2122, 211, 211};
    int pdgDgLctopiLambda[3] = {2122, 211, 211};
    int pdgDgK0s[2] = {211, 211};
    int pdgDgLambda[2] = {2122, 211};

    chCasc.fGenLabel = -1;
    chCasc.fDecay = kNone;
    int labtoK0s = cand->MatchToMC(4122, 310, pdgDgLctopK0s, pdgDgK0s, fMCArray, true);
    int labtoLambda = -1;
    if (labtoK0s < 0)
    {
        labtoLambda = cand->MatchToMC(4122, 3122, pdgDgLctopiLambda, pdgDgLambda, fMCArray, true);
        if(labtoLambda >= 0)
        {
            chCasc.fGenLabel = labtoLambda;
            chCasc.fDecay = kLctopK0s;
        }
    }
    else
    {
        chCasc.fGenLabel = labtoK0s;
        chCasc.fDecay = kLctopiLambda;
    }

    chCasc.fCandType = 0;
    int origin = -1;
    if(chCasc.fGenLabel >= 0){
        AliAODMCParticle* partD0 = dynamic_cast<AliAODMCParticle*>(fMCArray->At(chCasc.fGenLabel));
        origin = AliVertexingHFUtils::CheckOrigin(fMCArray, partD0, true);
        if(origin == 4)
        {
            chCasc.fCandType |= kSignal;
            chCasc.fCandType |= kPrompt;
        }
        else if(origin == 5)
        {
            chCasc.fCandType |= kSignal;
            chCasc.fCandType |= kFeedDown;
        }
        else
        {
            chCasc.fCandType |= kSignal; // no prompt, no feed-down --> weird stuff
        }
    }
    else {
        chCasc.fCandType |= kBackground;
    }

    return chCasc;
}

//________________________________________________________________________
GenCharmHadron AliAnalysisTaskSECharmTriggerStudy::FillCharmGen(AliAODMCParticle* part, int origin, int decay)
{
    GenCharmHadron genCharm;
    genCharm.fPt = part->Pt();
    genCharm.fY = part->Y();
    genCharm.fGenLabel = part->GetLabel();
    if(origin == 4)
        genCharm.fCandType |= kPrompt;
    else if(origin == 5)
        genCharm.fCandType |= kFeedDown;
    genCharm.fDecay = decay;

    return genCharm;
}

//________________________________________________________________________
bool AliAnalysisTaskSECharmTriggerStudy::RecalcOwnPrimaryVertex(AliAODRecoDecayHF* cand)
{
    AliAODVertex *recvtx = cand->RemoveDaughtersFromPrimaryVtx(fAOD);
    if(!recvtx){
        AliDebug(2,"Removal of daughter tracks failed");
        return false;
    }

    cand->SetOwnPrimaryVtx(recvtx);
    cand->RecalculateImpPars(recvtx, fAOD);
    delete recvtx;

    return true;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmTriggerStudy::CleanOwnPrimaryVertex(AliAODRecoDecayHF* cand, AliAODVertex* origvtx)
{
    cand->UnsetOwnPrimaryVtx();
    if(origvtx) {
      cand->SetOwnPrimaryVtx(origvtx);
      delete origvtx;
      origvtx = NULL;
    }
    cand->RecalculateImpPars(cand->GetPrimaryVtx(), fAOD);
}
