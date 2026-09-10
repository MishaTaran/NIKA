#ifndef SENSITIVEDETECTOR_HH
#define SENSITIVEDETECTOR_HH

#include <G4VSensitiveDetector.hh>
#include <G4RunManager.hh>
#include <G4AnalysisManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <G4Track.hh>
#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4ParticleDefinition.hh>
#include <SteppingAction.hh>
#include <unordered_map>
#include <G4Step.hh>
#include <G4SDManager.hh>
#include <G4TouchableHistory.hh>
#include <G4OpticalPhoton.hh>
#include <G4VPhysicalVolume.hh>
#include <G4AffineTransform.hh>
#include <G4ios.hh>

#include "SDHit.hh"

class SensitiveDetector : public G4VSensitiveDetector {
public:
    SensitiveDetector(const G4String &sdName, G4int detID, G4String detName);
    ~SensitiveDetector() override = default;

    void Initialize(G4HCofThisEvent *hce) override;
    G4bool ProcessHits(G4Step *step, G4TouchableHistory *) override;

    G4int GetHCID() const { return HCID; }
    SDHitCollection *GetHits() { return hits; }
    G4int GetDetID() const { return detID; }
    const G4String &GetDetName() const { return detName; }

private:
    SDHit *FindOrCreateHit(G4int volumeID);

    SDHitCollection *hits = nullptr;
    SDHitCollection *optHC  = nullptr;
    G4int HCID = -1;

    std::unordered_map<int, int> indexByVol;

    G4int detID = -1;
    G4String detName;

    G4bool isLED = false;
};


#endif //SENSITIVEDETECTOR_HH