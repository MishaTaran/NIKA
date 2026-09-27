#include "SteppingAction.hh"
#include "G4ios.hh"

SteppingAction::SteppingAction() {
    debugEventID = -1;
    debugEnabled = false;
    gammaLogEnabled = false;
    gammaLogEventID = -1;
}

SteppingAction::~SteppingAction() {
    // Пустой
}

G4String SteppingAction::HumanReadableProcess(const G4String& procName) {
    if (procName == "phot")          return "Photoelectric effect";
    if (procName == "compt")         return "Compton scattering";
    if (procName == "conv")          return "Pair production";
    if (procName == "Rayl")          return "Rayleigh scattering";
    if (procName == "gammaNuclear")  return "Gamma nuclear";
    if (procName == "eBrem")         return "Bremsstrahlung";
    if (procName == "msc")           return "Multiple scattering";
    if (procName == "Transportation")return "Transportation";
    if (procName == "none" || procName.empty()) return "none";
    return procName;
}

void SteppingAction::LogGammaInteraction(const G4Step* step,
                                         const G4String& volName,
                                         G4int trackID, G4int parentID,
                                         G4double edep, G4double t,
                                         const G4ThreeVector& x) {
    const G4int evtID = G4EventManager::GetEventManager()
                            ->GetConstCurrentEvent()->GetEventID();

    if (gammaLogEventID >= 0 && gammaLogEventID != evtID) return;

    const auto* post = step->GetPostStepPoint();
    const auto* pre  = step->GetPreStepPoint();
    const auto* proc = post->GetProcessDefinedStep();
    const G4String procName = proc ? proc->GetProcessName() : "none";

    const bool isGammaInteraction =
        procName == "phot"         ||
        procName == "compt"        ||
        procName == "conv"         ||
        procName == "Rayl"         ||
        procName == "gammaNuclear";

    if (!isGammaInteraction) return;

    const G4double eGamma = pre->GetKineticEnergy() / MeV;

    G4cout << "========================================\n"
           << " [Gamma interaction] Event " << evtID << "\n"
           << "   Process : " << HumanReadableProcess(procName)
           << " (" << procName << ")\n"
           << "   Volume  : " << volName << "\n"
           << "   TrackID : " << trackID
           << "   ParentID: " << parentID << "\n"
           << "   E_gamma : " << eGamma << " MeV\n"
           << "   Edep    : " << edep << " MeV\n"
           << "   Position: (" << x.x() << ", " << x.y() << ", " << x.z() << ") mm\n"
           << "   Time    : " << t << " ns\n"
           << "========================================"
           << G4endl;
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    const auto* post = step->GetPostStepPoint();
    const auto* pre  = step->GetPreStepPoint();
    const auto* postProc = post->GetProcessDefinedStep();
    const auto* track = step->GetTrack();

    const G4ThreeVector x = post->GetPosition() / mm;
    const G4double t = post->GetGlobalTime() / ns;
    const G4double edep = step->GetTotalEnergyDeposit() / MeV;
    const G4double stepLen = step->GetStepLength() / mm;

    const auto* touch = post->GetTouchable();
    G4String volName = "World";
    int copyNo = -1;
    if (touch && touch->GetVolume()) {
        volName = touch->GetVolume()->GetName();
        copyNo = touch->GetVolume()->GetCopyNo();
    }

    const G4int trackID = track->GetTrackID();
    const G4int parentID = track->GetParentID();
    const G4double trackE = track->GetKineticEnergy() / MeV;
    const G4String particleName = track->GetDefinition()->GetParticleName();
    const G4int pdg = track->GetDefinition()->GetPDGEncoding();

    // === логирование гамма-взаимодействий в терминал ===
    if (gammaLogEnabled && track->GetDefinition() == G4Gamma::Gamma()) {
        LogGammaInteraction(step, volName, trackID, parentID, edep, t, x);
    }

    auto* ea = static_cast<EventAction*>(
        G4EventManager::GetEventManager()->GetUserEventAction());
    if (!ea) return;

    auto secs = step->GetSecondaryInCurrentStep();
    if (secs && !secs->empty()) {
        for (size_t i = 0; i < secs->size(); ++i) {
            const auto* sc = (*secs)[i];
            const auto* cp = sc->GetCreatorProcess();

            G4String secParticleName = sc->GetDefinition()->GetParticleName();
            G4int secPdg = sc->GetDefinition()->GetPDGEncoding();
            G4String processName = cp ? cp->GetProcessName() : "unknown";

            std::string procName = processName;
            G4bool isDeltaElectron = (secPdg == 11 &&
                                      (processName == "eIoni" ||
                                       processName == "ionIoni" ||
                                       processName == "hIoni" ||
                                       processName == "muIoni" ||
                                       procName.find("Ioni") != std::string::npos));

            if (isDeltaElectron) {
                G4double secE_MeV = sc->GetKineticEnergy() / MeV;
                G4ThreeVector secPos = sc->GetPosition() / mm;
                G4ThreeVector secDir = sc->GetMomentumDirection();

                G4double primaryE = trackE + secE_MeV + edep;

                InteractionRec rec;
                rec.trackID = trackID;
                rec.parentID = parentID;
                rec.process = processName;
                rec.volumeName = volName;
                rec.volumeID = copyNo;
                rec.pos_mm = x;
                rec.t_ns = t;
                rec.secIndex = static_cast<G4int>(i);
                rec.secPDG = secPdg;
                rec.secName = secParticleName;
                rec.secE_MeV = secE_MeV;
                rec.secDir = secDir;
                ea->interBuf.emplace_back(std::move(rec));
            }
        }
    }
}