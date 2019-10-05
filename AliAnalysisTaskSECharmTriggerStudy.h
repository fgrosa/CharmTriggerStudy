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
    Double32_t fDecayLengthXY;          //[0.0,6.5536,16]
    Double32_t fNormDecayLength;        //[0.0,102.4.,10]
    Double32_t fNormDecayLengthXY;      //[0.0,102.4,10]
    Double32_t fImpParProd;             //[-3276.8,3276.8,16]
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
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
    Double32_t fDecayLengthXY;          //[0.0,6.5535,16]
    Double32_t fNormDecayLength;        //[0.0,102.3.,10]
    Double32_t fNormDecayLengthXY;      //[0.0,102.3,10]
    Double32_t fSigmaVtx;               //[0.0,0.08190,12]
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
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
    Double32_t fDecayLengthXYD0;        //[0.0,6.5535,16]
    Double32_t fNormDecayLengthD0;      //[0.0,102.3.,10]
    Double32_t fNormDecayLengthXYD0;    //[0.0,102.3,10]
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
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
};

struct GenCharmHadron
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
        kSignal     = BIT(0),
        kBackground = BIT(1),
        kPrompt     = BIT(2),
        kFeedDown   = BIT(3)
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
        kLctopiLambda
    };

    enum kSystem
    {
        kpp,
        kPbPb
    };

    AliAnalysisTaskSECharmTriggerStudy();
    AliAnalysisTaskSECharmTriggerStudy(const char *name = "CharmTriggerStudy");
    virtual ~AliAnalysisTaskSECharmTriggerStudy();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);

    void Enable2Prongs(bool enable = true)  {fEnable2Prongs = enable;}
    void Enable3Prongs(bool enable = true)  {fEnable3Prongs = enable;}
    void EnableDstars(bool enable = true)   {fEnableDstars = enable;}
    void EnableCascades(bool enable = true) {fEnableCascades = enable;}

    void SetSystem(int system = kpp) {fSystem = system;}

private:

    Charm2Prong FillCharm2Prong(AliAODRecoDecayHF2Prong* cand);
    Charm3Prong FillCharm3Prong(AliAODRecoDecayHF3Prong* cand);
    Dstar FillDstar(AliAODRecoCascadeHF* cand, AliAODRecoDecayHF2Prong* dau);
    CharmCascade FillCharmCascade(AliAODRecoCascadeHF* cand, AliAODv0* dau);
    GenCharmHadron FillCharmGen(AliAODMCParticle* part, int origin, int decay);
    bool RecalcOwnPrimaryVertex(AliAODRecoDecayHF* cand);
    void CleanOwnPrimaryVertex(AliAODRecoDecayHF* cand, AliAODVertex* origvtx);

    TList *fOutput;                         //!<! List of output histograms
    TH1F* fHistNEvents;                     //!<! Histogram for event info
    TTree *fRecoTree;                       //!<! Output tree with reco candidates
    TTree *fGenTree;                        //!<! Output tree with generated particles

    AliEventCuts fEventCuts;                /// object for event selection
    int fSystem;                            /// system (pp or PbPb)
    AliAODEvent *fAOD;                      //!<! AOD event
    int fAODProtection;                     ///protection for delta AOD mismatch
    TClonesArray* fMCArray;                 //!<! MC array
    float fRecoZvtx;                        /// Z of the reconstructed primary vertex
    float fGenZvtx;                         /// Z of the generated primary vertex

    vector<Charm2Prong> fCharm2Prong;       /// vector of 2 prongs
    vector<Charm3Prong> fCharm3Prong;       /// vector of 3 prongs
    vector<Dstar> fDstar;                   /// vector of Dstar
    vector<CharmCascade> fCharmCascade;     /// vector of cascades
    vector<GenCharmHadron> fGenCharmHadron; /// vector of generated charm hadrons

    bool fEnable2Prongs;                    /// flag to enable 2-prong branch
    bool fEnable3Prongs;                    /// flag to enable 3-prong branch
    bool fEnableDstars;                     /// flag to enable Dstar branch
    bool fEnableCascades;                   /// flag to enable cascade branch

    ClassDef(AliAnalysisTaskSECharmTriggerStudy, 1);
};

#endif