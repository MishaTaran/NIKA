#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include <G4UserEventAction.hh>
#include <G4ThreeVector.hh>
#include <G4String.hh>
#include <G4Event.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4SDManager.hh>
#include <G4HCofThisEvent.hh>
#include <G4SystemOfUnits.hh>
#include <cfloat>
#include <vector>
#include <map>

#include "Geometry.hh"
#include "RunAction.hh"
#include "Sizes.hh"
#include "Configuration.hh"
#include "AnalysisManager.hh"
#include "SDHit.hh"
#include "SiPMOpticalSD.hh"

class G4Event;
class Geometry;

struct PrimaryRec {
    int index = 0;
    int pdg = 0;
    G4String name;
    double E_MeV = 0.0;
    G4ThreeVector dir;
    G4ThreeVector pos_mm;
    double t0_ns = 0.0;
};

struct InteractionRec {
    int trackID = -1;
    int parentID = -1;
    G4String process;
    G4String volumeName;
    int volumeID = -1;
    G4ThreeVector pos_mm;
    double t_ns = 0.0;

    int secIndex = -1;
    int secPDG = 0;
    G4String secName;
    double secE_MeV = 0.0;
    G4ThreeVector secDir;
};

struct PhotonRec {
    G4int photonID = -1;
    G4String detName;
    G4int detCh;
    G4double energy;
    G4ThreeVector pos_mm;
};

class EventAction : public G4UserEventAction {
public:
    std::vector<PrimaryRec> primBuf;
    std::vector<InteractionRec> interBuf;
    std::vector<G4int> photonCountBuf{0, 0, 0};
    std::vector<PhotonRec> photonBuf;

    EventAction(AnalysisManager *, RunAction *);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event *) override;
    void EndOfEventAction(const G4Event *) override;

    [[nodiscard]] bool HasTOFAndNoAC() const {
        int nTrig = 0;
        for (int i = 0; i < Sizes::Trigger::nLayers; ++i)
            if (hasTriggerLayer[i]) nTrig++;

        return (nTrig >= 2) && !hasSideVeto && !hasBottomVeto && !hasUpperVeto;
    }

private:
    void WritePrimaries_(int eventID);
    void WriteSiPMFromSD_(int eventID);
    int WriteInteractions_(int eventID);
    int WritePhotonsCount_(int eventID);
    int WritePhotons_(int eventID);
    int WriteEdepFromSD_(const G4Event *evt, int eventID);

    // Методы для отметки детекторов
    void MarkTrigger() { hasTrigger = true; }
    void MarkVeto() {hasSideVeto = hasBottomVeto = hasUpperVeto = true; }
    void MarkSideVeto() { hasSideVeto = true; }
    void MarkBottomVeto() { hasBottomVeto = true; }
    void MarkUpperVeto() { hasUpperVeto = true; }
    void MarkTriggerOpt() { hasTriggerOpt = true; }
    void MarkVetoOpt() { hasSideVetoOpt = hasUpperVetoOpt = hasBottomVetoOpt = true; }
    void MarkSideVetoOpt() { hasSideVetoOpt = true; }
    void MarkUpperVetoOpt() { hasUpperVetoOpt = true; }
    void MarkBottomVetoOpt() { hasBottomVetoOpt = true; }

    AnalysisManager *analysisManager = nullptr;

    std::vector<std::tuple<G4String, int, G4String>> detMap;
    std::vector<int> HCIDs;

    int nPrimaries = 0;
    int nInteractions = 0;
    int nPhotons = 0;
    int nEdepHits = 0;

    RunAction* run = nullptr;

    // Флаги для энергетических детекторов
    bool hasTrigger = false;
    bool hasSideVeto = false;
    bool hasBottomVeto = false;
    bool hasUpperVeto = false;
    
    // Флаги для оптических детекторов
    bool hasTriggerOpt = false;
    bool hasSideVetoOpt = false;
    bool hasBottomVetoOpt = false;
    bool hasUpperVetoOpt = false;
    bool hasTriggerLayer[3] = {false};

    double stripE[3] = {0.0};
    std::map<G4String, double> detE;

    const double kVetoThreshold  = 0.2 * MeV;
    const double kPlaneThreshold = 0.15 * MeV;
};

#endif