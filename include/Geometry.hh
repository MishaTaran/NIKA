#ifndef GEOMETRY_HH
#define GEOMETRY_HH

#include <G4VUserDetectorConstruction.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4VisAttributes.hh>
#include <G4SystemOfUnits.hh>
#include <algorithm>
#include <string>

#include "Detector.hh"
#include "Sizes.hh"

class Geometry : public G4VUserDetectorConstruction {
public:
    Geometry();
    ~Geometry() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

private:
    G4NistManager* nist{};
    G4Material* worldMat{};

    G4LogicalVolume* worldLV{};
    G4VPhysicalVolume* worldPV{};

    Detector* detector{};
};

#endif // GEOMETRY_HH
