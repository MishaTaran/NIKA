#include "Detector.hh"
#include "Utils.hh"
#include "Configuration.hh"
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4MultiUnion.hh>
#include "G4UnionSolid.hh"

namespace
{
    inline G4double ZInInstrument(G4double zWorld) {
        return zWorld - Sizes::Instrument::sizeZ / 2.0;
    }
}

using namespace Sizes;

// Векторы для хранения указателей на Tyvek и соседние объемы
static std::vector<G4VPhysicalVolume*> g_tyvekPVs;
static std::vector<G4VPhysicalVolume*> g_triggerPVs;
static std::vector<G4VPhysicalVolume*> g_sideVetoPVs;
static std::vector<G4VPhysicalVolume*> g_tyvekBottomPVs;
static std::vector<G4VPhysicalVolume*> g_tyvekSidePVs;
static std::vector<G4VPhysicalVolume*> g_bottomVetoPVs;

Detector::Detector(G4LogicalVolume* wLV, G4NistManager* nistMan)
    : worldLV(wLV), nist(nistMan) {
    DefineMaterials();
    DefineVisual();
}

void Detector::DefineMaterials() {
    G4Element* elC = nist->FindOrBuildElement("C");
    G4Element* elH = nist->FindOrBuildElement("H");
    G4Element* elCs = nist->FindOrBuildElement("Cs");
    G4Element* elI = nist->FindOrBuildElement("I");
    G4Element* elTl = nist->FindOrBuildElement("Tl");
    G4Element* elO = nist->FindOrBuildElement("O");
    G4Element* elN = nist->FindOrBuildElement("N");

    const std::string base = "../OpticalParameters/";
    auto loadRIndex = [&](const std::string& prefix) -> Utils::Table {
        // rindex is dimensionless
        return Utils::ReadCSV(base + prefix + "_refractive_index.csv", /*valueScale=*/1.0, /*clampNonNegative=*/false);
    };
    auto loadAbsLengthMM = [&](const std::string& prefix) -> Utils::Table {
        // ABSLENGTH values are stored in mm in your files
        return Utils::ReadCSV(base + prefix + "_absorption_length.csv", /*valueScale=*/mm, /*clampNonNegative=*/true);
    };
    auto loadEmission2 = [&](const std::string& prefix,
                             const std::string& suffix = "_normalised_emission_intensity.csv")
        -> Utils::EmissionTables {
        auto e = Utils::ReadEmissionCSV(base + prefix + suffix, /*valueScale=*/1.0, /*clampNonNegative=*/true);
        Utils::NormalizeMaxToOne(e.c1);
        Utils::NormalizeMaxToOne(e.c2);
        return e;
    };
    auto loadConsts = [&](const std::string& prefix) -> Utils::ConstMap {
        return Utils::ReadConstFile(base + prefix + "_optical_consts.txt");
    };
    auto applyYieldScaleIfPresent = [&](Utils::ConstMap& c) {
        auto it = c.find("SCINTILLATIONYIELD");
        if (it != c.end() && Configuration::yieldScale > 0) {
            it->second /= static_cast<G4double>(Configuration::yieldScale);
        }
    };
    {
    SiPMMat = nist->FindOrBuildMaterial("G4_Si");
    vetoMat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    const std::string p = "Veto";
    const auto rindex = loadRIndex(p);
    const auto absl = loadAbsLengthMM(p);
    auto emission = loadEmission2(p);
    auto c = loadConsts(p);
    applyYieldScaleIfPresent(c);
    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", rindex.E, rindex.V, rindex.E.size());
    mpt->AddProperty("ABSLENGTH", absl.E, absl.V, absl.E.size());
    Utils::ApplyScintillation(vetoMat, mpt, c, emission.c1, emission.c2, true);
    Utils::ApplyBirksIfPresent(vetoMat, c);
    }
    alMat = nist->FindOrBuildMaterial("G4_Al");
    SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    Epoxy = nist->FindOrBuildMaterial("G4_POLYSTYRENE");

    tyvekMat = new G4Material("Tyvek", 0.38 * g / cm3, 2, kStateSolid);
    tyvekMat->AddElement(elC, 2);
    tyvekMat->AddElement(elH, 4);

    rubberMat = new G4Material("Rubber", 0.92 * g / cm3, 2, kStateSolid);
    rubberMat->AddElement(elC, 5);
    rubberMat->AddElement(elH, 8); 

    texMat = new G4Material("Tex", 1.85 * g/cm3, 2, kStateSolid);
    texMat->AddMaterial(SiO2, 0.675); 
    texMat->AddMaterial(Epoxy, 0.325);
    {
    csIMat = new G4Material("CsI", 4.51 * g/cm3, 3, kStateSolid);
    const G4double wTlcsI = 0.0008;
    const G4double scalecsI = 1.0 - wTlcsI;
    csIMat->AddElement(elCs, 0.511549 * scalecsI);
    csIMat->AddElement(elI, 0.488451 * scalecsI);
    csIMat->AddElement(elTl, wTlcsI);
    auto rindex = loadRIndex("CsI");
    auto absl = loadAbsLengthMM("CsI");
    auto emission = loadEmission2("CsI");
    auto c = loadConsts("CsI");
    applyYieldScaleIfPresent(c);  
    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", rindex.E, rindex.V, rindex.E.size());
    mpt->AddProperty("ABSLENGTH", absl.E, absl.V, absl.E.size());  
    auto it = c.find("SCINTILLATIONYIELD");
    Utils::ApplyScintillation(csIMat, mpt, c, emission.c1, emission.c2, true);
    Utils::ApplyBirksIfPresent(csIMat, c);
    }
    {
    SiPMEncapsulantMat = new G4Material("SiPMEncapsulant_DGEBA", 1.16*g/cm3, 3, kStateSolid);
    // Stoichiometry from molecular formula C21H24O4
    SiPMEncapsulantMat->AddElement(elC, 21);
    SiPMEncapsulantMat->AddElement(elH, 24);
    SiPMEncapsulantMat->AddElement(elO, 4);
    const auto rindex = loadRIndex("SiPM_Encapsulant");
    const auto absl = loadAbsLengthMM("SiPM_Encapsulant");
    Utils::ApplyMaterialTable(SiPMEncapsulantMat, rindex, &absl);
    }
    {
    glassMat = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
    const auto rindex = loadRIndex("Glass");
    const auto absl = loadAbsLengthMM("Glass");
    Utils::ApplyMaterialTable(glassMat, rindex, &absl);
    }
    {
        galacticMat = nist->FindOrBuildMaterial("G4_Galactic");
        auto tN = Utils::MakeConstantTable(1.0 * eV, 4.0 * eV, 1.0);
        auto tA = Utils::MakeConstantTable(1.0 * eV, 4.0 * eV, 1e6 * m);
        Utils::ApplyMaterialTable(galacticMat, tN, &tA);
    }
    // ============================================================
    // СОЗДАНИЕ ПОВЕРХНОСТИ TYVEK
    // ============================================================
        Utils::Table tyvekRefl;
        try {
            tyvekRefl = Utils::ReadCSV("../OpticalParameters/Tyvek_reflectivity.csv", 1.0, true);
        } catch (const std::exception& e) {
            tyvekRefl.E = {1.0*eV, 4.0*eV};
            tyvekRefl.V = {0.95, 0.95};
        }

        auto* tyvekMPT = new G4MaterialPropertiesTable();
        tyvekMPT->AddProperty("REFLECTIVITY", tyvekRefl.E.data(), tyvekRefl.V.data(), tyvekRefl.E.size());

        tyvekSurf = new G4OpticalSurface("TyvekSurface");
        tyvekSurf->SetModel(unified);

        if (Configuration::polishedTyvek) {
            tyvekSurf->SetType(dielectric_metal);
            tyvekSurf->SetFinish(polished);
            tyvekSurf->SetSigmaAlpha(0.0);
        } else {
            tyvekSurf->SetType(dielectric_dielectric);
            tyvekSurf->SetFinish(groundfrontpainted);
            tyvekSurf->SetSigmaAlpha(0.2);

            std::vector<G4double> E = {1.0 * eV, 4.0 * eV};
            std::vector<G4double> spike = {0.0, 0.0};
            std::vector<G4double> lobe = {0.02, 0.02};
            std::vector<G4double> back = {0.0, 0.0};
            tyvekMPT->AddProperty("SPECULARSPIKECONSTANT", E.data(), spike.data(), 2, true);
            tyvekMPT->AddProperty("SPECULARLOBECONSTANT", E.data(), lobe.data(), 2, true);
            tyvekMPT->AddProperty("BACKSCATTERCONSTANT", E.data(), back.data(), 2, true);
            G4cout << "[Detector] Tyvek surface: diffuse" << G4endl;
        }
        tyvekSurf->SetMaterialPropertiesTable(tyvekMPT);
    }

void Detector::DefineVisual() {
    visVeto = new G4VisAttributes(G4Color(0.0, 0.7, 0.0, 0.4));
    visVeto->SetForceSolid(true);
    visAl = new G4VisAttributes(G4Color(0.8, 0.8, 0.8, 1.0));
    visAl->SetForceSolid(true);
    visTyvek = new G4VisAttributes(G4Color(0.0, 0.0, 1.0));
    visTyvek->SetForceSolid(true);
    visCsI = new G4VisAttributes(G4Color(0.7, 0.0, 0.0, 0.3));
    visCsI->SetForceSolid(true);
    visTex = new G4VisAttributes(G4Color(1.0, 1.0, 0.0, 1.0));
    visTex->SetForceSolid(true);
    visRubber = new G4VisAttributes(G4Color(1.0, 0.5, 0.0, 1.0));
    visRubber->SetForceSolid(true);   
    visSiPM = new G4VisAttributes(G4Color(1.0, 0.5, 1.0, 0.8));
    visSiPM->SetForceSolid(true);
    visGlass = new G4VisAttributes(G4Color(0.0, 0.876, 0.96, 0.5));
    visGlass->SetForceSolid(true);
    visGalactic = new G4VisAttributes(G4Color(1.0, 1.0, 10, 0.0));
    visGalactic->SetForceSolid(true);
}

void Detector::CreateAluminumOpticalSurface() {
    if (!Configuration::useOptics) return;
    if (aluminumSurf) return;
    
    Utils::Table alRefl;
    try {
        alRefl = Utils::ReadCSV("../OpticalParameters/Aluminum_reflectivity.csv", 1.0, true);
    } catch (const std::exception& e) {
        alRefl.E = {1.0*eV, 4.0*eV};
        alRefl.V = {0.9, 0.9};
    }
    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("REFLECTIVITY", alRefl.E.data(), alRefl.V.data(), alRefl.E.size());
    aluminumSurf = new G4OpticalSurface("AluminumSurface");
    aluminumSurf->SetModel(unified);
    aluminumSurf->SetType(dielectric_metal);
    aluminumSurf->SetFinish(polished);
    aluminumSurf->SetSigmaAlpha(0.0);
    aluminumSurf->SetMaterialPropertiesTable(mpt);
}

void Detector::Construct() {
    ConstructShellAndContainer();
    ConstructSiPM();
    CreateAluminumOpticalSurface();
    ConstructSideVeto();
    ConstructVeto();
    ConstructTOF();
    PlaceAllSiPMs();
    ConstructOpticalSurfaces();
}
// ============================================================
// СОЗДАНИЕ SiPM
// ============================================================
void Detector::ConstructSiPM() {
    auto* SiPMFrameBase = new G4Box("SiPMFrameBase",
                            Trigger::spmXY / 2.0,
                            Trigger::spmXY / 2.0,
                            Trigger::spmZ / 2.0);
    auto* SiPMHole = new G4Box("SiPMHole",
                            Trigger::spmXY / 2.0 - Trigger::SiPMFrame,
                            Trigger::spmXY / 2.0 - Trigger::SiPMFrame,
                            Trigger::spmZ / 2.0 + 5 * mm);
    auto* SiPMFrame = new G4SubtractionSolid("SiPMFrame", SiPMFrameBase, SiPMHole);
    SiPMFrameLV = new G4LogicalVolume(SiPMFrame, galacticMat, "SiPMFrameLV");
    SiPMFrameLV->SetVisAttributes(visGalactic);
    auto* SiPMBody = new G4Box("SiPMBody",
                            Trigger::spmXY / 2.0 - Trigger::SiPMFrame, 
                            Trigger::spmXY / 2.0 - Trigger::SiPMFrame,
                            (Trigger::spmZ - Trigger::SiPMWindowThick) / 2.0);
    SiPMBodyLV = new G4LogicalVolume(SiPMBody, SiPMMat, "SiPMBodyLV");
    SiPMBodyLV->SetVisAttributes(visSiPM);
    auto* SiPMWindow = new G4Box("SiPMWindow",
                            Trigger::spmXY / 2.0 - Trigger::SiPMFrame,
                            Trigger::spmXY / 2.0 - Trigger::SiPMFrame,
                            Trigger::SiPMWindowThick / 2.0);
    SiPMWindowLV = new G4LogicalVolume(SiPMWindow, SiPMEncapsulantMat, "SiPMWindowLV");
    SiPMWindowLV->SetVisAttributes(visGlass);
    auto* SiPMPV = new G4Box("SiPMWindow",
                            Trigger::spmXY / 2.0,
                            Trigger::spmXY / 2.0,
                            Trigger::spmZ / 2.0);                        
    SiPMLV = new G4LogicalVolume(SiPMPV, galacticMat, "SiPMLV");
    SiPMLV->SetVisAttributes(visGalactic);
    new G4PVPlacement(nullptr,
                  G4ThreeVector(0, 0, 0),
                  SiPMFrameLV,
                  "SiPMFramePV",
                  SiPMLV,
                  false,
                  0, true);
    auto* bodyPVP = new G4PVPlacement(nullptr,
                  G4ThreeVector(0, 0, -Trigger::SiPMWindowThick / 2.0),
                  SiPMBodyLV,
                  "SiPMBodyPV",
                  SiPMLV,
                  false,
                  0, true);
    auto* windowPVP = new G4PVPlacement(nullptr,
                  G4ThreeVector(0, 0, (Trigger::spmZ - Trigger::SiPMWindowThick) / 2.0),
                  SiPMWindowLV,
                  "SiPMWindowPV",
                  SiPMLV,
                  false,
                  0, true);
    if (SiPMPhotocathodeSurf == nullptr) {
        G4cout << "[Detector] Creating SiPMPhotocathode surface..." << G4endl;
        SiPMPhotocathodeSurf = new G4OpticalSurface("SiPMPhotocathode");
        SiPMPhotocathodeSurf->SetModel(unified);
        SiPMPhotocathodeSurf->SetType(dielectric_metal);
        SiPMPhotocathodeSurf->SetFinish(polished);
        auto* mpt = new G4MaterialPropertiesTable();
        auto pde = Utils::ReadCSV("../OpticalParameters/SiPM_PDE.csv", 1.0, true);
        Utils::NormalizeMaxToOne(pde);
        mpt->AddProperty("EFFICIENCY", pde.E, pde.V, pde.E.size());
        std::vector refl(pde.E.size(), 0.02);
        mpt->AddProperty("REFLECTIVITY", pde.E, refl, pde.E.size());
        SiPMPhotocathodeSurf->SetMaterialPropertiesTable(mpt);
    }
    new G4LogicalBorderSurface("SiPM_Photocathode", windowPVP, bodyPVP, SiPMPhotocathodeSurf);
}
// ============================================================
// ПРЯМОЕ РАЗМЕЩЕНИЕ SiPM
// ============================================================
void Detector::PlaceAllSiPMs() {
    G4RotationMatrix* rotX = new G4RotationMatrix();
    rotX->rotateX(90 * deg);
    G4RotationMatrix* rotXNeg = new G4RotationMatrix();
    rotXNeg->rotateX(-90 * deg);
    G4RotationMatrix* rotY = new G4RotationMatrix();
    rotY->rotateY(90 * deg);
    G4RotationMatrix* rotYNeg = new G4RotationMatrix();
    rotYNeg->rotateY(-90 * deg);
    struct SiPMPosition {
        G4double x, y, z;
        G4RotationMatrix* rot;
        const char* name;
    };
    std::vector<SiPMPosition> allPositions; 
    // ============================================================
    // 1) TRIGGER SiPM (24 штуки: 3 слоя по 8 SiPM)
    // ============================================================
    {
        G4double zStart = VetoAC::DownVeloAll + Trigger::polkaSmallZ + Trigger::standZ / 2.0 - Trigger::spmZ / 2.0;
        
        for (G4int layer = 0; layer < 3; ++layer) {
            G4double zPos = zStart + layer * (2.0 * (Trigger::tyvek + Trigger::polkaSmallZ) + Trigger::thicknessTrig) + Trigger::spmZ / 2.0;
            
            G4double xOuter = Trigger::halfX + Trigger::spmZ / 2.0;
            G4double yInner = (Trigger::halfY - Trigger::spmZ) / 3.0 + Trigger::spmZ / 2.0;
            
            allPositions.push_back({ xOuter,  yInner, zPos, rotY, "Trigger"});
            allPositions.push_back({ xOuter, -yInner, zPos, rotY, "Trigger"});
            allPositions.push_back({-xOuter,  yInner, zPos, rotYNeg, "Trigger"});
            allPositions.push_back({-xOuter, -yInner, zPos, rotYNeg, "Trigger"});
            
            G4double xInner = (Trigger::halfX - Trigger::spmZ) / 3.0 + Trigger::spmZ / 2.0;
            G4double yOuter = Trigger::halfY + Trigger::spmZ / 2.0;
            
            allPositions.push_back({ xInner,  yOuter, zPos, rotXNeg, "Trigger"});
            allPositions.push_back({ xInner, -yOuter, zPos, rotX, "Trigger"});
            allPositions.push_back({-xInner,  yOuter, zPos, rotXNeg, "Trigger"});
            allPositions.push_back({-xInner, -yOuter, zPos, rotX, "Trigger"});
        }
    }
    // ============================================================
    // 2) SIDE VETO SiPM (8 штук)
    // ============================================================
    {
        G4double zPos = VetoAC::alTSideDown + VetoAC::PlataSide + Trigger::spmZ / 2.0;
        G4double pos = (VetoAC::sideVetoXY - Trigger::spmXY) / 3.0 + Trigger::spmXY / 2.0;
        G4double edgeX = Instrument::halfX - VetoAC::alTSide - Trigger::tyvek - VetoAC::thickness / 2.0;
        G4double edgeY = Instrument::halfY - VetoAC::alTSide - Trigger::tyvek - VetoAC::thickness / 2.0;
        
        std::vector<std::pair<G4double, G4double>> coords = {
            { pos,  edgeY}, {-pos,  edgeY},
            { pos, -edgeY}, {-pos, -edgeY},
            { edgeX,  pos}, { edgeX, -pos},
            {-edgeX,  pos}, {-edgeX, -pos}
        };
        
        for (const auto& coord : coords) {
            allPositions.push_back({coord.first, coord.second, zPos, nullptr, "SideVeto"});
        }
    }
    // ============================================================
    // 3) UPPER VETO SiPM (4 штуки)
    // ============================================================
    {
    G4double zPos = VetoAC::DownVeloAll + 6.0 * (Trigger::tyvek + Trigger::polkaSmallZ) + 
                    3.0 * Trigger::thicknessTrig + VetoAC::alTDown + Trigger::tyvek + 
                    VetoAC::thickness + Trigger::spmZ / 2.0;
    G4double half = VetoAC::sideVetoXY / 2.0;
    G4RotationMatrix* rotY180 = new G4RotationMatrix();
    rotY180->rotateX(180 * deg);
    std::vector<std::pair<G4double, G4double>> coords = {
        { half,  half}, {-half,  half},
        { half, -half}, {-half, -half}
    };
    for (const auto& coord : coords) {
        allPositions.push_back({coord.first, coord.second, zPos, rotY180, "UpperVeto"});
    }
}
    // ============================================================
    // 4) BOTTOM VETO SiPM (4 штуки)
    // ============================================================
    {
        G4double zPos = VetoAC::alTDownZ + VetoAC::alTDown + VetoAC::Plata + Trigger::tyvek / 2.0 - 0.1 * mm;
        std::vector<std::pair<G4double, G4double>> coords = {
            {VetoAC::sideVetoXY/2.0, VetoAC::sideVetoXY/2.0}, {-VetoAC::sideVetoXY/2.0, VetoAC::sideVetoXY/2.0},
            {VetoAC::sideVetoXY/2.0, -VetoAC::sideVetoXY/2.0}, {-VetoAC::sideVetoXY/2.0, -VetoAC::sideVetoXY/2.0}
        }; 
        for (const auto& coord : coords) {
            allPositions.push_back({coord.first, coord.second, zPos, nullptr, "BottomVeto"});
        }
    }
    // ============================================================
    // РАЗМЕЩАЕМ SiPM НАПРЯМУЮ В cubeOuterLV
    // ============================================================
    std::vector<G4VPhysicalVolume*> siPMPVs;
    for (size_t i = 0; i < allPositions.size(); ++i) {
        const auto& p = allPositions[i];
        G4ThreeVector sipmPos(p.x * mm, p.y * mm, ZInInstrument(p.z));
        G4String name = G4String(p.name) + "SiPM_PVP_" + std::to_string(i);
        G4VPhysicalVolume* siPMPV = new G4PVPlacement(
            p.rot,
            sipmPos,
            SiPMLV,
            name,
            cubeOuterLV,
            false,
            static_cast<G4int>(i),
            true
        );
        siPMPVs.push_back(siPMPV);
    }
    // ============================================================
    // ПРИМЕНЯЕМ ОПТИЧЕСКУЮ ПОВЕРХНОСТЬ К КАЖДОМУ SiPM
    // ============================================================
    if (Configuration::useOptics && SiPMPhotocathodeSurf) {
        int nSurfaces = 0;
        for (size_t i = 0; i < siPMPVs.size(); ++i) {
            G4VPhysicalVolume* siPMPV = siPMPVs[i];  
            new G4LogicalBorderSurface(
                "SiPMToWorld_" + std::to_string(i),
                siPMPV,
                cubeOuterPV,
                SiPMPhotocathodeSurf
            );
            nSurfaces++;
        }
    }
}

void Detector::ConstructTOF() {
    const G4bool checkOverlaps = true;
    G4double z = VetoAC::DownVeloAll;
    // ============================================================
    // 1 Выступ маленький
    // ============================================================
    auto* outerpolkasmall = new G4Box("PolkaOuter",
                            Trigger::halfX + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY + Trigger::standrubXY + Trigger::standXY,
                            Trigger::halfY + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY + Trigger::standrubXY + Trigger::standXY,
                            Trigger::polkaSmallZ / 2.0);
    auto* innerpolkasmall = new G4Box("PolkaInner",
                            Trigger::halfX + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY + Trigger::standrubXY + Trigger::standXY - Trigger::polkaSmallXY,
                            Trigger::halfY + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY + Trigger::standrubXY + Trigger::standXY - Trigger::polkaSmallXY,
                            Trigger::polkaSmallZ);
    auto* polkaSolidSmall = new G4SubtractionSolid("PolkaSmall", outerpolkasmall, innerpolkasmall);
    auto* polkaSmallLV = new G4LogicalVolume(polkaSolidSmall, alMat, "PolkaSmallLV");
    polkaSmallLV->SetVisAttributes(visAl);
    // ============================================================
    // 2 Cтойка
    // ============================================================
    auto* outerstand = new G4Box("StandOuter",
                            Trigger::halfX + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY + Trigger::standrubXY + Trigger::standXY,
                            Trigger::halfY + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY + Trigger::standrubXY + Trigger::standXY,
                            Trigger::standZ / 2.0);
    auto* innerstand = new G4Box("StandInner",
                            Trigger::halfX + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY + Trigger::standrubXY,
                            Trigger::halfY + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY + Trigger::standrubXY,
                            Trigger::standZ);
    auto* standSolid = new G4SubtractionSolid("Stand", outerstand, innerstand);
    auto* standLV = new G4LogicalVolume(standSolid, alMat, "StandLV");
    standLV->SetVisAttributes(visAl);
    // ============================================================
    // 3 Cтойка резиновая
    // ============================================================
    auto* outerstandrub = new G4Box("StandrubOuter",
                            Trigger::halfX + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY + Trigger::standrubXY,
                            Trigger::halfY + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY + Trigger::standrubXY,
                            Trigger::standZ / 2.0);
    auto* innerstandrub = new G4Box("StandrubInner",
                            Trigger::halfX + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY,
                            Trigger::halfY + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY,
                            Trigger::standZ);
    auto* standrubSolid = new G4SubtractionSolid("Standrub", outerstandrub, innerstandrub);
    auto* standrubLV = new G4LogicalVolume(standrubSolid, rubberMat, "StandrubLV");
    standrubLV->SetVisAttributes(visRubber); 
    // ============================================================
    // 4 Cтойка с электроникой
    // ============================================================
    auto* outerstandel = new G4Box("StandelOuterEl",
                            Trigger::halfX + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY,
                            Trigger::halfY + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ + Trigger::plataTOFXY, 
                            Trigger::standZ / 2.0);
    auto* innerstandel = new G4Box("StandelInnerEl",
                            Trigger::halfX + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ,
                            Trigger::halfY + Trigger::tyvek + Trigger::standAlXY + Trigger::spmZ,
                            Trigger::standZ);
    auto* standElSolid = new G4SubtractionSolid("Standel", outerstandel, innerstandel);
    auto* standElLV = new G4LogicalVolume(standElSolid, texMat, "StandelLV");
    standElLV->SetVisAttributes(visTex);   
    // ============================================================
    // 6 Пластина алюминия
    // ============================================================
    auto* outerstandal = new G4Box("StandalOuter",
                            Trigger::halfX + Trigger::standAlXY + Trigger::tyvek,
                            Trigger::halfY + Trigger::standAlXY + Trigger::tyvek,
                            Trigger::standZ / 2.0);
    auto* innerstandal = new G4Box("StandalInner",
                            Trigger::halfX + Trigger::tyvek,
                            Trigger::halfY + Trigger::tyvek,
                            Trigger::standZ);
    auto* standAlSolidBase = new G4SubtractionSolid("StandAl", outerstandal, innerstandal);
    auto* squareHoleY = new G4Box("SquareHole",
                             Trigger::spmXY / 2.0,
                             Instrument::halfY,
                             Trigger::spmXY / 2.0);
    auto* squareHoleX = new G4Box("SquareHole",
                             Instrument::halfX,
                             Trigger::spmXY / 2.0,
                             Trigger::spmXY / 2.0); 
    G4MultiUnion* holesUnion = new G4MultiUnion("HolesUnion");
    holesUnion->AddNode(squareHoleY, G4Translate3D((Trigger::halfX - Trigger::spmZ) / 3.0 + Trigger::spmZ / 2.0, 0, 0));
    holesUnion->AddNode(squareHoleY, G4Translate3D(-((Trigger::halfX - Trigger::spmZ) / 3.0 + Trigger::spmZ / 2.0), 0, 0));
    holesUnion->AddNode(squareHoleX, G4Translate3D(0, ((Trigger::halfY - Trigger::spmZ) / 3.0 + Trigger::spmZ / 2.0), 0));
    holesUnion->AddNode(squareHoleX, G4Translate3D(0, -((Trigger::halfY - Trigger::spmZ) / 3.0 + Trigger::spmZ / 2.0), 0));
    holesUnion->Voxelize();
    G4VSolid* standAlSolid = new G4SubtractionSolid("standAlWithAllHoles",
                                               standAlSolidBase,
                                               holesUnion,
                                               nullptr,
                                               G4ThreeVector(0, 0, 0));
    auto* standAlLV = new G4LogicalVolume(standAlSolid, alMat, "StandAlLV");
    standAlLV->SetVisAttributes(visAl);
    // ============================================================
    // 7 СЦИНТИЛЛЯТОР 5 мм + ТАЙВИК
    // ============================================================
    auto* outerBox = new G4Box("TyvekOuter",
                           Trigger::halfX + Trigger::tyvek,
                           Trigger::halfY + Trigger::tyvek,
                           (Trigger::thicknessTrig / 2.0) + Trigger::tyvek);
    auto* innerBox = new G4Box("TyvekInner",
                           Trigger::halfX,
                           Trigger::halfY,
                           Trigger::thicknessTrig / 2.0);
    auto* tyvekSolid = new G4SubtractionSolid("TyvekSolid", outerBox, innerBox);
    auto* tyvekSolidWithHoles = new G4SubtractionSolid("TyvekWithHoles",
                                                tyvekSolid,
                                                holesUnion,
                                                nullptr,
                                                G4ThreeVector(0, 0, 0));
    auto* tyvekLV = new G4LogicalVolume(tyvekSolidWithHoles, tyvekMat, "TyvekLV");
    tyvekLV->SetVisAttributes(visTyvek);

    auto* scintTrigger = new G4Box("Trigger",
                            Trigger::halfX,
                            Trigger::halfY,
                            Trigger::thicknessTrig / 2.0);
    triggerLV = new G4LogicalVolume(scintTrigger, csIMat, "TriggerLV");
    triggerLV->SetVisAttributes(visCsI);

    for (G4int i = 0; i <= 2; i++) {
        // ============================================================
        // ЦИКЛ Выступ маленький
        // ============================================================
        new G4PVPlacement(nullptr,
                          G4ThreeVector(0, 0, ZInInstrument(z + Trigger::polkaSmallZ / 2.0)),
                          polkaSmallLV,
                          "PolkaSmall",
                          cubeOuterLV,
                          false,
                          i,
                          checkOverlaps);
        z += Trigger::polkaSmallZ;
        
        // ============================================================
        // ЦИКЛ Cтойка
        // ============================================================
        new G4PVPlacement(nullptr,
                          G4ThreeVector(0, 0, ZInInstrument(z + Trigger::standZ / 2.0)),
                          standLV,
                          "Stand",
                          cubeOuterLV,
                          false,
                          i,
                          checkOverlaps);
        new G4PVPlacement(nullptr,
                          G4ThreeVector(0, 0, ZInInstrument(z + Trigger::standZ / 2.0)),
                          standrubLV,
                          "Standrub",
                          cubeOuterLV,
                          false,
                          i,
                          checkOverlaps);
        new G4PVPlacement(nullptr,
                          G4ThreeVector(0, 0, ZInInstrument(z + Trigger::standZ / 2.0)),
                          standElLV,
                          "StandEl",
                          cubeOuterLV,
                          false,
                          i,
                          checkOverlaps);
        new G4PVPlacement(nullptr,
                          G4ThreeVector(0, 0, ZInInstrument(z + Trigger::standZ / 2.0)),
                          standAlLV,
                          "StandAl",
                          cubeOuterLV,
                          false,
                          i,
                          checkOverlaps);
        
        // ============================================================
        // ЦИКЛ СЦИНТИЛЛЯТОР 5 мм + ТАЙВИК
        // ============================================================
        triggerPV = new G4PVPlacement(nullptr,
                              G4ThreeVector(0, 0, ZInInstrument(z + Trigger::thicknessTrig/2.0 + Trigger::tyvek)),
                              triggerLV,
                              "TriggerPV",
                              cubeOuterLV,
                              false,
                              i,
                              checkOverlaps);
        g_triggerPVs.push_back(triggerPV);

        tyvekPV = new G4PVPlacement(nullptr,
                              G4ThreeVector(0, 0, ZInInstrument(z + Trigger::tyvek + Trigger::thicknessTrig / 2.0)),
                              tyvekLV,
                              "Tyvek",
                              cubeOuterLV,
                              false,
                              i,
                              checkOverlaps);
        g_tyvekPVs.push_back(tyvekPV); 
        z += Trigger::thicknessTrig + 2 * Trigger::tyvek;
        // ============================================================
        // ЦИКЛ Выступ маленький (закрывающий)
        // ============================================================
        new G4PVPlacement(nullptr,
                          G4ThreeVector(0, 0, ZInInstrument(z + Trigger::polkaSmallZ / 2.0)),
                          polkaSmallLV,
                          "PolkaSmall",
                          cubeOuterLV,
                          false,
                          i,
                          checkOverlaps);
        z += Trigger::polkaSmallZ;
    }
}

void Detector::ConstructSideVeto() {
    const G4bool checkOverlaps = true;
    // ============================================================
    // 1 Наружняя пластина из алюминия
    // ============================================================
    auto* outeral = new G4Box("SidePLate1Outer",
                            Instrument::halfX,
                            Instrument::halfY,
                            (Instrument::sizeZ - VetoAC::alTSideDown) / 2.0);
    auto* inneral = new G4Box("SidePLate1Inner",
                            Instrument::halfX - VetoAC::alTSide,
                            Instrument::halfY - VetoAC::alTSide,
                            Instrument::sizeZ - VetoAC::alTSideDown);
    auto* plateSideSolid = new G4SubtractionSolid("SidePLate1", outeral, inneral);
    auto* plateSideLV = new G4LogicalVolume(plateSideSolid, alMat, "SidePLateLV");
    plateSideLV->SetVisAttributes(visAl);
    G4VPhysicalVolume* sidePlateAlPVP = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument((Instrument::sizeZ + VetoAC::alTSideDown) / 2.0)),
                      plateSideLV,
                      "SidePLate",
                      cubeOuterLV,
                      false,
                      0,
                      checkOverlaps);
    // ============================================================
    // 2 Нижняя пластина
    // ============================================================
    auto* outerAlDown = new G4Box("SideVetoOuterDown",
                            Instrument::halfX,
                            Instrument::halfY,
                            VetoAC::alTSideDown / 2.0);
    auto* innerAlDown = new G4Box("SideVetoInnerDown",
                            Instrument::halfX - VetoAC::alTDownXY,
                            Instrument::halfY - VetoAC::alTDownXY,
                            VetoAC::alTSideDown);
    auto* plateSideDownSolid = new G4SubtractionSolid("SidePLateDown", outerAlDown, innerAlDown);
    auto* plateSideDownLV = new G4LogicalVolume(plateSideDownSolid, alMat, "SidePLateDownLV");
    plateSideDownLV->SetVisAttributes(visAl);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(VetoAC::alTSideDown / 2.0)),
                      plateSideDownLV,
                      "SidePLateDownAl",
                      cubeOuterLV,
                      false,
                      0,
                      checkOverlaps);
    // ============================================================
    // 3 Нижняя плата
    // ============================================================
    auto* outerDownPlata = new G4Box("SideVetoOuterDown",
                            Instrument::halfX - VetoAC::alTSide,
                            Instrument::halfY - VetoAC::alTSide,
                            VetoAC::PlataSide / 2.0);
    auto* innerDownPlata = new G4Box("SideVetoInnerDown",
                            Instrument::halfX - VetoAC::alTSide - 2.0 * Trigger::tyvek - VetoAC::thickness,
                            Instrument::halfY - VetoAC::alTSide - 2.0 * Trigger::tyvek - VetoAC::thickness,
                            VetoAC::PlataSide);
    auto* electroSideDownSolid = new G4SubtractionSolid("SideElectroSideDown", outerDownPlata, innerDownPlata);
    auto* electroSideDownLV = new G4LogicalVolume(electroSideDownSolid, texMat, "SideElectroSideDownLV");
    electroSideDownLV->SetVisAttributes(visTex);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(VetoAC::alTSideDown + VetoAC::PlataSide / 2.0)),
                      electroSideDownLV,
                      "SideElectro",
                      cubeOuterLV,
                      false,
                      0,
                      checkOverlaps);
    // ============================================================
    // 4 Боковое вето
    // ============================================================
    auto* outer = new G4Box("SideVetoOuter",
                            Instrument::halfX - VetoAC::alTSide - Trigger::tyvek,
                            Instrument::halfY - VetoAC::alTSide - Trigger::tyvek,
                            VetoAC::VetoHeight / 2.0);
    auto* inner = new G4Box("SideVetoInner",
                            Instrument::halfX - VetoAC::alTSide - Trigger::tyvek - VetoAC::thickness,
                            Instrument::halfY - VetoAC::alTSide - Trigger::tyvek - VetoAC::thickness,
                            VetoAC::VetoHeight);
    auto* shell = new G4SubtractionSolid("SideVetoShell", outer, inner);
    vetoLV = new G4LogicalVolume(shell, vetoMat, "SideVetoLV");
    vetoLV->SetVisAttributes(visVeto);
    sideVetoPV = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(VetoAC::VetoHeight / 2.0 + VetoAC::alTSideDown + VetoAC::PlataSide + Trigger::spmZ + Trigger::tyvek)),
                      vetoLV,
                      "SideVeto",
                      cubeOuterLV,
                      false,
                      0,
                      checkOverlaps);
    g_sideVetoPVs.push_back(sideVetoPV);
    // ============================================================
    // 5 Тайвик
    // ============================================================
    auto* outerBox = new G4Box("TyvekOuter",
                            Instrument::halfX - VetoAC::alTSide,
                            Instrument::halfY - VetoAC::alTSide,
                            VetoAC::VetoHeight/2.0 + Trigger::tyvek);
    auto* innerBox = new G4Box("TyvekOuter",
                            Instrument::halfX - VetoAC::alTSide - 2.0 * Trigger::tyvek - VetoAC::thickness,
                            Instrument::halfY - VetoAC::alTSide - 2.0 * Trigger::tyvek - VetoAC::thickness,
                            Instrument::sizeZ);
    auto* tyvekSolidBase1 = new G4SubtractionSolid("TyvekSolidBase", outerBox, innerBox);
    auto* tyvekSolidBase = new G4SubtractionSolid("TyvekSolidBase", tyvekSolidBase1, shell);
    auto* squareHole = new G4Box("SquareHole",
                            Trigger::spmXY / 2.0,
                            Trigger::spmXY / 2.0,
                            10 * mm);
    G4double pos = (VetoAC::sideVetoXY - Trigger::spmXY) / 3.0 + Trigger::spmXY / 2.0;
    G4double edgeX = Instrument::halfX - VetoAC::alTSide - Trigger::tyvek - VetoAC::thickness / 2.0;
    G4double edgeY = Instrument::halfY - VetoAC::alTSide - Trigger::tyvek - VetoAC::thickness / 2.0;
    G4MultiUnion* holesUnion = new G4MultiUnion("HolesUnion");
    holesUnion->AddNode(squareHole, G4Translate3D(pos, edgeY, -40 * mm));
    holesUnion->AddNode(squareHole, G4Translate3D(pos, -edgeY, -40 * mm));
    holesUnion->AddNode(squareHole, G4Translate3D(-pos, edgeY, -40 * mm));
    holesUnion->AddNode(squareHole, G4Translate3D(-pos, -edgeY, -40 * mm));
    holesUnion->AddNode(squareHole, G4Translate3D(edgeX, pos, -40 * mm));
    holesUnion->AddNode(squareHole, G4Translate3D(edgeX, -pos, -40 * mm));
    holesUnion->AddNode(squareHole, G4Translate3D(-edgeX, pos, -40 * mm));
    holesUnion->AddNode(squareHole, G4Translate3D(-edgeX, -pos, -40 * mm));
    holesUnion->Voxelize();
    G4VSolid* tyvekSolid = new G4SubtractionSolid("TyvekSolid",
                                               tyvekSolidBase,
                                               holesUnion,
                                               nullptr,
                                               G4ThreeVector(0, 0, 0));
    auto* tyvekLV = new G4LogicalVolume(tyvekSolid, tyvekMat, "Tyvek");
    tyvekLV->SetVisAttributes(visTyvek);
    tyvekSidePV = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(VetoAC::VetoHeight / 2.0 + VetoAC::alTSideDown + VetoAC::PlataSide + Trigger::spmZ + Trigger::tyvek)),
                      tyvekLV,
                      "Tyvek",
                      cubeOuterLV,
                      false,
                      3,
                      checkOverlaps);
    g_tyvekSidePVs.push_back(tyvekSidePV);
    // ============================================================
    // 6 Верхняя пластина
    // ============================================================
    auto* outerAlUp = new G4Box("SideVetoOuterUp",
                            Instrument::halfX - VetoAC::alTSide,
                            Instrument::halfY - VetoAC::alTSide,
                            VetoAC::alTSideUp / 2.0);
    auto* innerAlUp = new G4Box("SideVetoInnerUp",
                            Instrument::halfX - VetoAC::alTSide - 2.0 * Trigger::tyvek - VetoAC::thickness,
                            Instrument::halfY - VetoAC::alTSide - 2.0 * Trigger::tyvek - VetoAC::thickness,
                            VetoAC::alTSideUp);
    auto* plateSideUpSolid = new G4SubtractionSolid("SidePLateUp", outerAlUp, innerAlUp);
    auto* plateSideUpLV = new G4LogicalVolume(plateSideUpSolid, alMat, "SidePLateUpLV");
    plateSideUpLV->SetVisAttributes(visAl);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(VetoAC::alTSideDown + VetoAC::PlataSide + VetoAC::VetoHeight + VetoAC::alTSideUp / 2.0 + Trigger::spmZ + 2.0 * Trigger::tyvek)),
                      plateSideUpLV,
                      "SidePLatealTSideUp",
                      cubeOuterLV,
                      false,
                      0,
                      checkOverlaps);
    // ============================================================
    // 7 Внутренняя пластина
    // ============================================================
    auto* outeral2 = new G4Box("SidePLate2Outer",
                            Instrument::halfX - VetoAC::alTSide - 2.0 * Trigger::tyvek - VetoAC::thickness,
                            Instrument::halfY - VetoAC::alTSide - 2.0 * Trigger::tyvek - VetoAC::thickness,
                            (Instrument::sizeZ - VetoAC::alTSideDown) / 2.0);
    auto* inneral2 = new G4Box("SidePLate2Inner",
                            Instrument::halfX - 2.0 * (VetoAC::alTSide + Trigger::tyvek) - VetoAC::thickness,
                            Instrument::halfY - 2.0 * (VetoAC::alTSide + Trigger::tyvek) - VetoAC::thickness,
                            Instrument::sizeZ - VetoAC::alTSideDown);
    auto* plateSide2Solid = new G4SubtractionSolid("SidePLate2", outeral2, inneral2);
    auto* plateSide2LV = new G4LogicalVolume(plateSide2Solid, alMat, "SidePLate2");
    plateSide2LV->SetVisAttributes(visAl);
    
    G4VPhysicalVolume* sidePlate2PVP = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument((Instrument::sizeZ + VetoAC::alTSideDown) / 2.0)),
                      plateSide2LV,
                      "SidePLate",
                      cubeOuterLV,
                      false,
                      1,
                      checkOverlaps);
}

void Detector::ConstructVeto() {
    const G4bool checkOverlaps = true;
    G4double z = VetoAC::alTDownZ - 0.05 * mm;
    // ============================================================
    // 1 Алюминьевая пластина нижняя
    // ============================================================
    auto* alDownSolid = new G4Box("Al1",
                               VetoAC::downVetoXY + Trigger::tyvek,
                               VetoAC::downVetoXY + Trigger::tyvek,
                               VetoAC::alTDown / 2.0);
    auto* AlPlateLV = new G4LogicalVolume(alDownSolid, alMat, "AlPlateLV");
    AlPlateLV->SetVisAttributes(visAl);
    G4VPhysicalVolume* alPlatePVP0 = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z + VetoAC::alTDown / 2.0)),
                      AlPlateLV,
                      "AlPlate",
                      cubeOuterLV,
                      false,
                      0,
                      checkOverlaps);
    z += VetoAC::alTDown;
    // ============================================================
    // 2 Выступ
    // ============================================================
    auto* outerpolka = new G4Box("SidePolkaOuter",
                            VetoAC::downVetoXY + Trigger::tyvek + VetoAC::polkaXY,
                            VetoAC::downVetoXY + Trigger::tyvek + VetoAC::polkaXY,
                            VetoAC::polkaZ / 2.0);
    auto* innerpolka = new G4Box("SidePolkaInner",
                            VetoAC::downVetoXY + Trigger::tyvek,
                            VetoAC::downVetoXY + Trigger::tyvek,
                            VetoAC::polkaZ);
    auto* polkaSideSolid = new G4SubtractionSolid("SidePolka", outerpolka, innerpolka);
    auto* polkaSideLV = new G4LogicalVolume(polkaSideSolid, alMat, "SidePolkaLV");
    polkaSideLV->SetVisAttributes(visAl);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z + VetoAC::polkaZ / 2.0)),
                      polkaSideLV,
                      "SidePolka",
                      cubeOuterLV,
                      false,
                      0,
                      checkOverlaps);
    // ============================================================
    // 3 Электроника
    // ============================================================
    auto* outerElectroDown = new G4Box("ElectroDownOuter",
                            VetoAC::downVetoXY,
                            VetoAC::downVetoXY,
                            VetoAC::Plata / 2.0);
    auto* innerElectroDown = new G4Box("ElectroDownInner",
                            VetoAC::downVetoXY - 1 * mm,
                            VetoAC::downVetoXY - 1 * mm,
                            VetoAC::Plata);
    auto* electroDownSolid = new G4SubtractionSolid("ElectroDown", outerElectroDown, innerElectroDown);
    auto* electroDownLV = new G4LogicalVolume(electroDownSolid, texMat, "ElectroDown");
    electroDownLV->SetVisAttributes(visTex);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z + VetoAC::Plata / 2.0)),
                      electroDownLV,
                      "ElectroDown",
                      cubeOuterLV,
                      false,
                      0,
                      checkOverlaps);
    z += VetoAC::Plata + VetoAC::thickness/2.0 + Trigger::tyvek;
    // ============================================================
    // 4 Тайвик
    // ============================================================
    auto* outerBottomBox = new G4Box("TyvekOuter",
                            VetoAC::downVetoXY + Trigger::tyvek,
                            VetoAC::downVetoXY + Trigger::tyvek,
                            VetoAC::thickness/2.0 + Trigger::tyvek);
    auto* innerBottomBox = new G4Box("TyvekOuter",
                            VetoAC::downVetoXY,
                            VetoAC::downVetoXY,
                            VetoAC::thickness/2.0);
    auto* tyvekSolidBase = new G4SubtractionSolid("TyvekSolidBase", outerBottomBox, innerBottomBox);
    auto* squareHole = new G4Box("SquareHole",
                            Trigger::spmXY / 2.0,
                            Trigger::spmXY / 2.0,
                            5 * mm);
    G4MultiUnion* holesUnionB = new G4MultiUnion("HolesUnion");
    holesUnionB->AddNode(squareHole, G4Translate3D(VetoAC::sideVetoXY/2.0, VetoAC::sideVetoXY / 2.0, -5 * mm));
    holesUnionB->AddNode(squareHole, G4Translate3D(-VetoAC::sideVetoXY/2.0, VetoAC::sideVetoXY/2.0, -5 * mm));
    holesUnionB->AddNode(squareHole, G4Translate3D(VetoAC::sideVetoXY/2.0, -VetoAC::sideVetoXY/2.0, -5 * mm));
    holesUnionB->AddNode(squareHole, G4Translate3D(-VetoAC::sideVetoXY/2.0, -VetoAC::sideVetoXY/2.0, -5 * mm));
    holesUnionB->Voxelize();
    G4MultiUnion* holesUnionU = new G4MultiUnion("HolesUnion");
    holesUnionU->AddNode(squareHole, G4Translate3D(VetoAC::sideVetoXY/2.0, VetoAC::sideVetoXY / 2.0, 5 * mm));
    holesUnionU->AddNode(squareHole, G4Translate3D(-VetoAC::sideVetoXY/2.0, VetoAC::sideVetoXY/2.0, 5 * mm));
    holesUnionU->AddNode(squareHole, G4Translate3D(VetoAC::sideVetoXY/2.0, -VetoAC::sideVetoXY/2.0, 5 * mm));
    holesUnionU->AddNode(squareHole, G4Translate3D(-VetoAC::sideVetoXY/2.0, -VetoAC::sideVetoXY/2.0, 5 * mm));
    holesUnionU->Voxelize();
    G4VSolid* tyvekBottomSolid = new G4SubtractionSolid("TyvekSolid",
                                               tyvekSolidBase,
                                               holesUnionB,
                                               nullptr,
                                               G4ThreeVector(0, 0, 0));
    G4VSolid* tyvekUpperSolid = new G4SubtractionSolid("TyvekSolid",
                                               tyvekSolidBase,
                                               holesUnionU,
                                               nullptr,
                                               G4ThreeVector(0, 0, 0));
    auto* tyvekBottomLV = new G4LogicalVolume(tyvekBottomSolid, tyvekMat, "Tyvek");
    tyvekBottomLV->SetVisAttributes(visTyvek);
    auto* tyvekUpperLV = new G4LogicalVolume(tyvekUpperSolid, tyvekMat, "Tyvek");
    tyvekUpperLV->SetVisAttributes(visTyvek);
    tyvekBottomPV = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z)),
                      tyvekBottomLV,
                      "Tyvek",
                      cubeOuterLV,
                      false,
                      4,
                      checkOverlaps);
    g_tyvekBottomPVs.push_back(tyvekBottomPV);
    // ============================================================
    // 5 Нижнее антисовпадение
    // ============================================================
    auto* scintBottomVeto = new G4Box("ScintBottomVeto",
                                  VetoAC::downVetoXY,
                                  VetoAC::downVetoXY,
                                  VetoAC::thickness / 2.0);
    bottomVetoLV = new G4LogicalVolume(scintBottomVeto, vetoMat, "BottomVetoLV");
    bottomVetoLV->SetVisAttributes(visVeto);
    bottomVetoPV = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z)),
                      bottomVetoLV,
                      "BottomVeto",
                      cubeOuterLV,
                      false,
                      0,
                      checkOverlaps);
    g_bottomVetoPVs.push_back(bottomVetoPV);
    z += VetoAC::thickness/2.0 + Trigger::tyvek;
    // ============================================================
    // 7 Алюминьевая пластина
    // ============================================================
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z + VetoAC::alTDown / 2.0)),
                      AlPlateLV,
                      "AlPlate",
                      cubeOuterLV,
                      false,
                      1,
                      checkOverlaps);
    z += VetoAC::alTDown + 3.0 * Trigger::thicknessTrig + 6.0 * (Trigger::tyvek + Trigger::polkaSmallZ);
    // ============================================================
    // 8 Алюминьевая пластина
    // ============================================================
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z + VetoAC::alTDown / 2.0)),
                      AlPlateLV,
                      "AlPlate",
                      cubeOuterLV,
                      false,
                      2,
                      checkOverlaps);
    z += VetoAC::alTDown;
    // ============================================================
    // 9 Выступ
    // ============================================================
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z + VetoAC::polkaZ / 2.0)),
                      polkaSideLV,
                      "SidePolka",
                      cubeOuterLV,
                      false,
                      1,
                      checkOverlaps);
    z += Trigger::tyvek + VetoAC::thickness/2.0;
    // ===========================================================
    // 11 Тайвик
    // ============================================================
    tyvekUpperPV = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z)),
                      tyvekUpperLV,
                      "Tyvek",
                      cubeOuterLV,
                      false,
                      5,
                      checkOverlaps);
    g_tyvekBottomPVs.push_back(tyvekUpperPV);
    // ============================================================
    // 12 Верхнее антисовпадение
    // ============================================================
    auto* scintUpperVeto = new G4Box("ScintUpperVeto",
                                  VetoAC::downVetoXY,
                                  VetoAC::downVetoXY,
                                  VetoAC::thickness / 2.0);
    upperVetoLV = new G4LogicalVolume(scintUpperVeto, vetoMat, "UpperVetoLV");
    upperVetoLV->SetVisAttributes(visVeto);
    upperVetoPV = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z)),
                      upperVetoLV,
                      "UpperVeto",
                      cubeOuterLV,
                      false,
                      0,
                      checkOverlaps);
    g_bottomVetoPVs.push_back(upperVetoPV);
    z += Trigger::tyvek + VetoAC::thickness/2.0;
    // ============================================================
    // 13 Электроника
    // ============================================================
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z + VetoAC::Plata/2.0)),
                      electroDownLV,
                      "ElectroDown",
                      cubeOuterLV,
                      false,
                      1,
                      checkOverlaps);
    z += VetoAC::Plata;
    // ============================================================
    // 14 Алюминьевая пластина
    // ============================================================
    G4VPhysicalVolume* alPlatePVP3 = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, ZInInstrument(z + VetoAC::alTDown / 2.0)),
                      AlPlateLV,
                      "AlPlate",
                      cubeOuterLV,
                      false,
                      3,
                      checkOverlaps);
}
void Detector::ConstructShellAndContainer() {
    const G4bool checkOverlaps = true;
    auto* solid = new G4Box("InstrumentSolid",
                            Instrument::halfX,
                            Instrument::halfY,
                            Instrument::halfZ);
    cubeOuterLV = new G4LogicalVolume(solid, galacticMat, "InstrumentLV");
    cubeOuterLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    cubeOuterPV = new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, Instrument::sizeZ / 2.0),
                      cubeOuterLV,
                      "InstrumentPV",
                      worldLV,
                      false,
                      0,
                      checkOverlaps);
}
void Detector::AddBorderSurface(const G4String& name,
                                G4VPhysicalVolume* pvFrom,
                                G4VPhysicalVolume* pvTo,
                                G4OpticalSurface* surf) {
    if (!pvFrom || !pvTo || !surf) return;
    new G4LogicalBorderSurface(name, pvFrom, pvTo, surf);
}
void Detector::AddBidirectionalBorder(const G4String& nameAToB,
                                      const G4String& nameBToA,
                                      G4VPhysicalVolume* pvA,
                                      G4VPhysicalVolume* pvB,
                                      G4OpticalSurface* surf) {
    if (!pvA || !pvB || !surf) return;
    new G4LogicalBorderSurface(nameAToB, pvA, pvB, surf);
    new G4LogicalBorderSurface(nameBToA, pvB, pvA, surf);
}
void Detector::ConstructOpticalSurfaces() {
    const auto refl = Utils::ReadCSV("../OpticalParameters/Tyvek_reflectivity.csv", 1.0, true);
    auto* mpt = new G4MaterialPropertiesTable();

    auto* tyvekSurf = new G4OpticalSurface("TyvekSurface");
    tyvekSurf->SetModel(unified);

    if (Configuration::polishedTyvek) {
        tyvekSurf->SetType(dielectric_metal);
        tyvekSurf->SetFinish(polished);
        tyvekSurf->SetSigmaAlpha(0.0);
    } else {
        tyvekSurf->SetType(dielectric_dielectric);
        tyvekSurf->SetFinish(groundfrontpainted);
        tyvekSurf->SetSigmaAlpha(0.2);
        std::vector E = {1.0 * eV, 4.0 * eV};
        std::vector spike = {0.0, 0.0};
        std::vector lobe = {0.02, 0.02};
        std::vector back = {0.0, 0.0};
        mpt->AddProperty("SPECULARSPIKECONSTANT", E.data(), spike.data(), static_cast<G4int>(E.size()), true);
        mpt->AddProperty("SPECULARLOBECONSTANT", E.data(), lobe.data(), static_cast<G4int>(E.size()), true);
        mpt->AddProperty("BACKSCATTERCONSTANT", E.data(), back.data(), static_cast<G4int>(E.size()), true);
    }
    mpt->AddProperty("REFLECTIVITY", refl.E, refl.V, refl.E.size());
    tyvekSurf->SetMaterialPropertiesTable(mpt);
    // 1) Trigger <-> Tyvek
    AddBidirectionalBorder("TriggerToTyvek", "TyvekToTrigger", triggerPV, tyvekPV, tyvekSurf);
    // 2) BottomVeto <-> Tyvek
    AddBidirectionalBorder("BottomVetoToTyvek", "TyvekToBottomVeto", bottomVetoPV, tyvekBottomPV, tyvekSurf);
    // 3) UpperVeto <-> Tyvek
    AddBidirectionalBorder("BottomVetoToTyvek", "TyvekToBottomVeto", upperVetoPV, tyvekUpperPV, tyvekSurf);
    // 4) SideVeto <-> Tyvek
    AddBidirectionalBorder("SideVetoTyvek", "TyvekToSideVeto", sideVetoPV, tyvekSidePV, tyvekSurf);
    // 5) Trigger <-> SiPM
    AddBidirectionalBorder("TriggerToSiPM", "SiPMToTrigger", triggerPV, SiPMPVPV, tyvekSurf);
    // 6) BottomVeto <-> SiPM
    AddBidirectionalBorder("BottomVetoToSiPM", "SiPMToBottomVeto", bottomVetoPV, SiPMPVPV, tyvekSurf);
    // 7) UpperVeto <-> SiPM
    AddBidirectionalBorder("UpperVetoToSiPM", "SiPMToVetoSiPM", upperVetoPV, SiPMPVPV, tyvekSurf);
}