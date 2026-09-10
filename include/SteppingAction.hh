#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH

#include <G4Step.hh>
#include <G4VProcess.hh>
#include <G4EventManager.hh>
#include <G4TrackingManager.hh>
#include <G4RunManager.hh>
#include <G4ParticleDefinition.hh>
#include <G4TouchableHistory.hh>
#include <G4SystemOfUnits.hh>
#include <G4UserSteppingAction.hh>
#include <vector>
#include <fstream>
#include <map>

#include "EventAction.hh"


class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction();
    ~SteppingAction();
    
    void UserSteppingAction(const G4Step* step) override;
    
    void SetDebugEventID(G4int id) { debugEventID = id; }
    void EnableDebug(bool enable) { debugEnabled = enable; }

private:
    G4int debugEventID = -1;
    G4bool debugEnabled = false;
};

#endif