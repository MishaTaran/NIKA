#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH

#include <G4Step.hh>
#include <G4VProcess.hh>
#include <G4EventManager.hh>
#include <G4TrackingManager.hh>
#include <G4RunManager.hh>
#include <G4ParticleDefinition.hh>
#include <G4Gamma.hh>
#include <G4TouchableHistory.hh>
#include <G4SystemOfUnits.hh>
#include <G4UserSteppingAction.hh>
#include <G4ios.hh>
#include <vector>
#include <fstream>
#include <map>
#include <string>

#include "EventAction.hh"


class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction();
    ~SteppingAction();

    void UserSteppingAction(const G4Step* step) override;

    void SetDebugEventID(G4int id) { debugEventID = id; }
    void EnableDebug(bool enable) { debugEnabled = enable; }

    // === логирование гамма-взаимодействий в терминал ===
    void EnableGammaLog(bool enable) { gammaLogEnabled = enable; }
    void SetGammaLogEventID(G4int id) { gammaLogEventID = id; }

private:
    G4int debugEventID = -1;
    G4bool debugEnabled = false;

    G4bool gammaLogEnabled = false;
    G4int  gammaLogEventID = -1;   // -1 = все события

    static G4String HumanReadableProcess(const G4String& procName);

    void LogGammaInteraction(const G4Step* step,
                             const G4String& volName,
                             G4int trackID, G4int parentID,
                             G4double edep, G4double t,
                             const G4ThreeVector& x);
};

#endif // STEPPINGACTION_HH