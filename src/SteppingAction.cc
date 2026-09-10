#include "SteppingAction.hh"
#include "G4ios.hh"

SteppingAction::SteppingAction() {
    debugEventID = -1;
    debugEnabled = false;
}

SteppingAction::~SteppingAction() {
    // Пустой
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    const auto* post = step->GetPostStepPoint();
    const auto* pre = step->GetPreStepPoint();
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
    
    auto* ea = static_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction());
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