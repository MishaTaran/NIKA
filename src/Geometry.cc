#include "Geometry.hh"
#include "SensitiveDetector.hh"
#include "SiPMOpticalSD.hh"
#include "Configuration.hh"

using namespace Sizes;

Geometry::Geometry() {
    nist = G4NistManager::Instance();
}

G4VPhysicalVolume* Geometry::Construct() {
    worldMat = nist->FindOrBuildMaterial("G4_Galactic");

    const G4double halfWorld = 5.0 * 100;

    auto* worldBox = new G4Box("World",
                               halfWorld,
                               halfWorld,
                               halfWorld);

    worldLV = new G4LogicalVolume(worldBox, worldMat, "WorldLV");
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    worldPV = new G4PVPlacement(nullptr,
                                G4ThreeVector(0, 0, 0),
                                worldLV,
                                "WorldPV",
                                nullptr,
                                false,
                                0,
                                false);
    // auto* temp = new G4Box("hehe", 10, 10, 10);
    // auto* tempLV = new G4LogicalVolume(temp, worldMat, "heheLV");
    // new G4PVPlacement(nullptr, G4ThreeVector(), tempLV, "hehePVP", worldLV, false, 0, true);

    detector = new Detector(worldLV, nist);
    detector->Construct();

    return worldPV;
}

void Geometry::ConstructSDandField() {
    G4cout << "\n========== CONSTRUCTING SD AND FIELD ==========" << G4endl;
    
    auto* sdManager = G4SDManager::GetSDMpointer();

    auto* triggerLV = detector->GetTriggerLV();
    if (triggerLV) {
        auto* triggerSD = new SensitiveDetector("TriggerSD", 0, "Trigger");
        sdManager->AddNewDetector(triggerSD);
        triggerLV->SetSensitiveDetector(triggerSD);
        G4cout << "[Geometry] ✓ Registered Trigger SD" << G4endl;
    }
    auto* sideVetoLV = detector->GetVetoLV();
    if (sideVetoLV) {
        auto* sideVetoSD = new SensitiveDetector("SideVetoSD", 2, "SideVeto");
        sdManager->AddNewDetector(sideVetoSD);
        sideVetoLV->SetSensitiveDetector(sideVetoSD);
        G4cout << "[Geometry] ✓ Registered SideVeto SD" << G4endl;
    }
    auto* upperVetoLV = detector->GetUpperVetoLV();
    if (upperVetoLV) {
        auto* upperVetoSD = new SensitiveDetector("UpperVetoSD", 3, "UpperVeto");
        sdManager->AddNewDetector(upperVetoSD);
        upperVetoLV->SetSensitiveDetector(upperVetoSD);
        G4cout << "[Geometry] ✓ Registered UpperVeto SD" << G4endl;
    }
    auto* bottomVetoLV = detector->GetBottomVetoLV();
    if (bottomVetoLV) {
        auto* bottomVetoSD = new SensitiveDetector("BottomVetoSD", 4, "BottomVeto");
        sdManager->AddNewDetector(bottomVetoSD);
        bottomVetoLV->SetSensitiveDetector(bottomVetoSD);
        G4cout << "[Geometry] ✓ Registered BottomVeto SD" << G4endl;
    }
    if (Configuration::useOptics) {
        auto* sipmSD = new SiPMOpticalSD("SiPMOpticalSD");
        auto* SipmLV = detector->GetSiPMLV();
        sipmSD->SetSiPMLV(SipmLV);
        sdManager->AddNewDetector(sipmSD);
        SipmLV->SetSensitiveDetector(sipmSD);
    }
}