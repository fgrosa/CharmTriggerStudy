#ifndef ALIANALYSISTASKSECHARMTRIGGERSTUDY_H
#define ALIANALYSISTASKSECHARMTRIGGERSTUDY_H

//**************************************************************************************
// \class AliAnalysisTaskSECharmTriggerStudy
// \brief task that produces an output tree for the charm trigger studies
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
// M. Puccio, maximiliano.puccio@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include <TList.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include "AliAnalysisTaskSE.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODv0.h"
#include "AliEventCuts.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsLctoV0.h"

using std::vector;

struct Charm2Prong
{
    float fInvMassD0;                   /// inv mass of D0 hypothesis
    float fInvMassD0bar;                /// inv mass of D0bar hypothesis
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fY;                      //[-1.023,1.023,11]
    Double32_t fCosP;                   //[0.67233000,1.,15]
    Double32_t fCosPXY;                 //[0.67233000,1.,15]
    Double32_t fDecayLength;            //[0.0,6.5536,16]
    Double32_t fNormDecayLengthXY;      //[0.0,102.4,10]
    Double32_t fImpParProd;             //[-0.32768,0.32768,16]
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct Charm3Prong
{
    float fInvMassDplus;
    float fInvMassDstoKKpi;
    float fInvMassDstopiKK;
    float fInvMassLctopKpi;
    float fInvMassLctopiKp;
    float fInvMassPhiKKpi;
    float fInvMassPhipiKK;
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fYDplus;                 //[-1.023,1.023,11]
    Double32_t fYDs;                    //[-1.023,1.023,11]
    Double32_t fYLc;                    //[-1.023,1.023,11]
    Double32_t fCosP;                   //[0.67233000,1.,15]
    Double32_t fCosPXY;                 //[0.67233000,1.,15]
    Double32_t fDecayLength;            //[0.0,6.5535,16]
    Double32_t fNormDecayLengthXY;      //[0.0,102.3,10]
    Double32_t fSigmaVtx;               //[0.0,0.08190,12]
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct Dstar
{
    float fInvMass;
    float fInvMassD0;
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fY;                      //[-1.023,1.023,11]
    Double32_t fCosPD0;                 //[0.67233000,1.,15]
    Double32_t fCosPXYD0;               //[0.67233000,1.,15]
    Double32_t fDecayLengthD0;          //[0.0,6.5535,16]
    Double32_t fNormDecayLengthXYD0;    //[0.0,102.3,10]
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct CharmCascade
{
    float fInvMassLctopK0s;
    float fInvMassLctopiLambda;
    float fInvMassK0s;
    float fInvMassLambda;
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fY;                      //[-1.023,1.023,11]
    Double32_t fCosPV0;                 //[0.67233000,1.,15]
    Double32_t fCosPXYV0;               //[0.67233000,1.,15]
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct Beauty2Prong
{
    float fInvMassBplustoD0pi;
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fY;                      //[-1.023,1.023,11]
    Double32_t fCosP;                   //[0.67233000,1.,15]
    Double32_t fCosPXY;                 //[0.67233000,1.,15]
    Double32_t fDecayLength;            //[0.0,6.5536,16]
    Double32_t fNormDecayLengthXY;      //[0.0,102.4,10]
    Double32_t fImpParProd;             //[-0.32768,0.32768,16]
    Double32_t fPtD0;                   //[0.0,65.535,16]
    Double32_t fCosPD0;                 //[0.67233000,1.,15]
    Double32_t fCosPXYD0;               //[0.67233000,1.,15]
    Double32_t fDecayLengthD0;          //[0.0,6.5536,16]
    Double32_t fNormDecayLengthXYD0;    //[0.0,102.4,10]
    Double32_t fImpParProdD0;           //[-0.32768,0.32768,16]
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct GenHadron
{
    float fPt;
    float fY;
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
};

class AliAnalysisTaskSECharmTriggerStudy : public AliAnalysisTaskSE
{
public:
    enum kCandType
    {
        kSignal      = BIT(0),
        kBackground  = BIT(1),
        kPrompt      = BIT(2),
        kFeedDown    = BIT(3),
        kHasDauInAcc = BIT(4),
        kIsInFidAcc  = BIT(5)
    };

    enum kDecay
    {
        kNone,
        kDzerotoKpi,
        kDzerotopiK,
        kDplustoKpipi,
        kDstoKKpi,
        kDstopiKK,
        kDplustoKKpi,
        kDplustopiKK,
        kDstartoKpipi,
        kLctopKpi,
        kLctopiKp,
        kLctopK0s,
        kLctopiLambda,
        kBplustoD0pi
    };

    enum kSelBit
    {
        kDzerotoKpiCuts      = BIT(0),
        kDzerotopiKCuts      = BIT(1),
        kDplustoKpipiCuts    = BIT(2),
        kDstartoKpipiCuts    = BIT(3),
        kDstoKKpiCuts        = BIT(4),
        kDstopiKKCuts        = BIT(5),
        kLctopKpiCuts        = BIT(6),
        kLctopiKpCuts        = BIT(7),
        kLctoV0bachCuts      = BIT(8),
        kBplustoD0piCuts     = BIT(9),
        kDzerotoKpiCutsPID   = BIT(10),
        kDzerotopiKCutsPID   = BIT(11),
        kDplustoKpipiCutsPID = BIT(12),
        kDstartoKpipiCutsPID = BIT(13),
        kDstoKKpiCutsPID     = BIT(14),
        kDstopiKKCutsPID     = BIT(15),
        kLctopKpiCutsPID     = BIT(16),
        kLctopiKpCutsPID     = BIT(17),
        kLctoV0bachCutsPID   = BIT(18),
        kDzerotoKpiFidAcc    = BIT(19),
        kDplustoKpipiFidAcc  = BIT(20),
        kDstartoKpipiFidAcc  = BIT(21),
        kDstoKKpiFidAcc      = BIT(22),
        kLctopKpiFidAcc      = BIT(23),
        kLctoV0bachFidAcc    = BIT(24),
        kBplustoD0piFidAcc   = BIT(25)
    };

    enum kSystem
    {
        kpp,
        kPbPb
    };

    AliAnalysisTaskSECharmTriggerStudy();
    AliAnalysisTaskSECharmTriggerStudy(const char *name, TList *cutlist);
    virtual ~AliAnalysisTaskSECharmTriggerStudy();

    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);

    void Enable2Prongs(bool enable = true)               {fEnable2Prongs = enable;}
    void Enable3Prongs(bool enableDplus = true,
                       bool enableDs = true,
                       bool enableLc = true)             {fEnable3Prongs = 0; if(enableDplus) fEnable3Prongs |= BIT(0); if(enableDs) fEnable3Prongs |= BIT(1); if(enableLc) fEnable3Prongs |= BIT(2);}
    void EnableDstars(bool enable = true)                {fEnableDstars = enable;}
    void EnableCascades(bool enable = true)              {fEnableCascades = enable;}
    void EnableBplus(bool enable = true)                 {fEnableBplus = enable;}
    void SetFillOnlySignal(bool fillonlysignal = true)   {fFillOnlySignal = fillonlysignal;}

    void SetSystem(int system = kpp)                     {fSystem = system;}
    void ApplyCuts(bool applycuts = true)                {fApplyCuts = applycuts;}

private:

    void FillCharm2Prong(AliAODRecoDecayHF2Prong* cand, int issel);
    void FillCharm3Prong(AliAODRecoDecayHF3Prong* cand, bool isselDplus, int isselDs, int isselLc);
    void FillDstar(AliAODRecoCascadeHF* cand, AliAODRecoDecayHF2Prong* dau, bool issel);
    void FillCharmCascade(AliAODRecoCascadeHF* cand, AliAODv0* dau, int issel);
    void FillBeauty2Prong(AliAODRecoDecayHF2Prong* cand, AliAODRecoDecayHF2Prong* dau, bool issel);
    void FillGenerated(AliAODMCParticle* part, int origin, int decay, bool aredauinacc);
    bool RecalcOwnPrimaryVertex(AliAODRecoDecayHF* cand);
    void CleanOwnPrimaryVertex(AliAODRecoDecayHF* cand, AliAODVertex* origvtx);
    bool AreDauInAcc(int nProng, int *labDau);
    bool IsInFiducialAcceptance(double pt, double y);
    AliAODVertex* ReconstructDisplVertex(const AliVVertex *primary, TObjArray *tracks, double bField, double dispersion);

    TList *fOutput;                             //!<! List of output histograms
    TH1F* fHistNEvents;                         //!<! Histogram for event info
    TTree *fRecoTree;                           //!<! Output tree with reco candidates
    TTree *fGenTree;                            //!<! Output tree with generated particles

    AliEventCuts fEventCuts;                    /// object for event selection
    int fSystem;                                /// system (pp or PbPb)
    AliAODEvent *fAOD;                          //!<! AOD event
    int fAODProtection;                         /// protection for delta AOD mismatch
    TClonesArray* fMCArray;                     //!<! MC array
    float fRecoZvtx;                            /// Z of the reconstructed primary vertex
    float fGenZvtx;                             /// Z of the generated primary vertex
    int fNtracklets;                            /// number of tracklets in |eta| < 1

    vector<Charm2Prong> fCharm2Prong;           /// vector of charm 2 prongs
    vector<Charm3Prong> fCharm3Prong;           /// vector of charm 3 prongs
    vector<Dstar> fDstar;                       /// vector of Dstar
    vector<CharmCascade> fCharmCascade;         /// vector of charm cascades
    vector<Beauty2Prong> fBeauty2Prong;         /// vector of beauty 2 prongs
    vector<GenHadron> fGenHadron;               /// vector of generated charm/beauty hadrons

    bool fEnable2Prongs;                        /// flag to enable 2-prong branch
    int fEnable3Prongs;                         /// flag to enable 3-prong branch (with D+ and/or Ds+ and/or Lc)
    bool fEnableDstars;                         /// flag to enable Dstar branch
    bool fEnableCascades;                       /// flag to enable cascade branch
    bool fEnableBplus;                          /// flag to enable B+
    bool fFillOnlySignal;                       /// flag to fill only signal

    AliRDHFCutsD0toKpi* fCutsD0toKpi;           /// cut object for D0->Kpi
    AliRDHFCutsDplustoKpipi* fCutsDplustoKpipi; /// cut object for D+->Kpipi
    AliRDHFCutsDStartoKpipi* fCutsDstartoKpipi; /// cut object for D*+->D0pi->Kpipi
    AliRDHFCutsDstoKKpi* fCutsDstoKKpi;         /// cut object for Ds+->phipi->KKpi
    AliRDHFCutsLctopKpi* fCutsLctopKpi;         /// cut object for Lc->pKpi
    AliRDHFCutsLctoV0* fCutsLctoV0bach;         /// cut object for Lc->V0bachelor

    bool fApplyCuts;                            /// flag to enable cuts application
    TList* fListCuts;                           /// list of cut objects

    ClassDef(AliAnalysisTaskSECharmTriggerStudy, 1);
};

#endif