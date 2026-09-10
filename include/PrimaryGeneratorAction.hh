#ifndef PRMIARYGENERATIONACTION_HH
#define PRMIARYGENERATIONACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ParticleGun.hh>
#include <G4ThreeVector.hh>
#include <G4String.hh>
#include <utility>
#include <G4VVisManager.hh>
#include <G4Event.hh>
#include <G4Circle.hh>
#include <G4EventManager.hh>
#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <cmath>
#include <fstream>
#include <numeric>
#include <utility>
#include <sstream>
#include <filesystem>

#include "EventAction.hh"
#include "Geometry.hh"
#include "Flux/Flux.hh"
#include "Flux/UniformFlux.hh"
#include "Flux/PLAWFlux.hh"
#include "Flux/COMPFlux.hh"
#include "Flux/SEPFlux.hh"
#include "Flux/TableFlux.hh"
#include "Flux/GalacticFlux.hh"
#include "Configuration.hh"

namespace fs = std::filesystem;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction(G4String , const G4String &, G4double cThreshold, bool saveSeeds = false);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event *evt) override;
    
    void SetSaveSeeds(bool save) { saveSeeds = save; }

private:
    G4ParticleGun *particleGun = nullptr;

    G4double radius;
    G4ThreeVector center;
    G4ThreeVector detectorHalfSize;

    G4String fluxDirection;
    ParticleInfo pInfo{};

    Flux *flux;

    G4double eTriggerThreshold;
    
    bool saveSeeds = false;
    
    const std::string eventsDir = "events/";

    void GenerateOnSphere(G4ThreeVector &pos, G4ThreeVector &dir) const;
};

#endif