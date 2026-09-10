#ifndef DETECTOR_HH
#define DETECTOR_HH

#include <G4VisAttributes.hh>
#include <G4SubtractionSolid.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4Element.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4LogicalVolume.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4MultiUnion.hh>

#include "Sizes.hh"

class Detector {
public:
    Detector(G4LogicalVolume* worldLV, G4NistManager* nistMan);
    ~Detector() = default;

    void Construct();
    G4LogicalVolume* GetTriggerLV() const { return triggerLV; }
    G4LogicalVolume* GetUpperVetoLV() const { return upperVetoLV; }
    G4LogicalVolume* GetBottomVetoLV() const { return bottomVetoLV; }
    G4LogicalVolume* GetVetoLV() const { return vetoLV; }
    G4LogicalVolume* GetSiPMLV() const { return SiPMWindowLV; }
    G4VPhysicalVolume* GetCubeOuterPV() const { return cubeOuterPV; }

private:
    G4LogicalVolume* worldLV{};
    G4NistManager* nist{};

    // ===== Материалы =====
    G4Material* galacticMat{};
    G4Material* SiPMMat{};
    G4Material* SiPMEncapsulantMat{};
    G4Material* vetoMat{};
    G4Material* alMat{};
    G4Material* tyvekMat{};
    G4Material* rubberMat{};
    G4Material* SiO2{};
    G4Material* Epoxy{};
    G4Material* texMat{};
    G4Material* csIMat{};
    G4Material* glassMat{};

    // ===== Логические объемы =====
    G4LogicalVolume* cubeOuterLV{};
    G4VPhysicalVolume* cubeOuterPV{};
    G4LogicalVolume* bottomVetoLV{};
    G4LogicalVolume* triggerLV{};
    G4LogicalVolume* vetoLV{};
    G4LogicalVolume* upperVetoLV{};
    G4LogicalVolume* SiPMBodyLV{};
    G4LogicalVolume* SiPMWindowLV{};
    G4LogicalVolume* SiPMFrameLV{};
    G4LogicalVolume* SiPMLV{};

    G4VPhysicalVolume* triggerPV{};
    G4VPhysicalVolume* tyvekPV{};
    G4VPhysicalVolume* bottomVetoPV{};
    G4VPhysicalVolume* tyvekBottomPV{};
    G4VPhysicalVolume* tyvekSidePV{};
    G4VPhysicalVolume* tyvekUpperPV{};
    G4VPhysicalVolume* upperVetoPV{};
    G4VPhysicalVolume* sideVetoPV{};
    G4VPhysicalVolume* SiPMPVPV{};

    // ===== Tyvek поверхность =====
    G4OpticalSurface* tyvekSurf{nullptr};

    // ===== Алюминиевая поверхность =====
    G4OpticalSurface* aluminumSurf{nullptr};

    // ===== Визуализация =====
    G4VisAttributes* visGalactic{};
    G4VisAttributes* visVeto{};
    G4VisAttributes* visAl{};
    G4VisAttributes* visTyvek{};
    G4VisAttributes* visRubber{};
    G4VisAttributes* visTex{};
    G4VisAttributes* visCsI{};
    G4VisAttributes* visSiPM{};
    G4VisAttributes* visGlass{};

    // ===== Методы =====
    void DefineMaterials();
    void DefineVisual();
    void ConstructVeto();
    void ConstructTOF();
    void ConstructSideVeto();
    void AddBorderSurface();
    void AddBorderSurface(const G4String& name,
                          G4VPhysicalVolume* pvFrom,
                          G4VPhysicalVolume* pvTo,
                          G4OpticalSurface* surf);
    void AddBidirectionalBorder(const G4String& nameAToB,
                                const G4String& nameBToA,
                                G4VPhysicalVolume* pvA,
                                G4VPhysicalVolume* pvB,
                                G4OpticalSurface* surf);
    void ConstructOpticalSurfaces();
    void ConstructShellAndContainer();
    
    // ===== SiPM методы =====
    void ConstructSiPM();
    void PlaceAllSiPMs();
    void CreateAluminumOpticalSurface();
    G4OpticalSurface* SiPMPhotocathodeSurf = nullptr;
};

#endif