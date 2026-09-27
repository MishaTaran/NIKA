#include "PrimaryGeneratorAction.hh"
#include <G4RunManager.hh>
#include <cstdlib>
#include <filesystem>

using namespace Configuration;

namespace fs = std::filesystem;

PrimaryGeneratorAction::PrimaryGeneratorAction(G4String fDir, const G4String& fluxType, const G4double cThreshold, bool saveSeeds)
    : particleGun(new G4ParticleGun(1)),
      center(G4ThreeVector(0, 0, Sizes::Instrument::halfZ)),
      detectorHalfSize(G4ThreeVector(Sizes::Instrument::halfX,
                                     Sizes::Instrument::halfY,
                                     Sizes::Instrument::halfZ)),
      fluxDirection(std::move(fDir)),
      eTriggerThreshold(cThreshold),
      saveSeeds(saveSeeds) {
    radius = sqrt(detectorHalfSize.x() * detectorHalfSize.x() + detectorHalfSize.y() * detectorHalfSize.y() + detectorHalfSize.z() * detectorHalfSize.z()) + 5 * mm;

    std::vector<G4String> fluxDirList = {
        "isotropic", "isotropic_up", "isotropic_down", "vertical_up", "vertical_down", "horizontal"
    };
    if (std::find(fluxDirList.begin(), fluxDirList.end(), fluxDirection) == fluxDirList.end()) {
        G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "FluxDirection", FatalException,
                    ("Flux direction is not implemented: " + fluxDirection +
                        ".\nAvailable flux directions: isotropic, isotropic_up, isotropic_down, vertical_up," +
                        " vertical_down, horizontal").c_str());
    }
    std::vector<G4String> fluxTypeList = {"Uniform", "PLAW", "COMP", "SEP", "Galactic", "Table"};
    if (std::find(fluxTypeList.begin(), fluxTypeList.end(), fluxType) == fluxTypeList.end()) {
        G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "FluxType", FatalException,
                    ("Flux type not found: " + fluxType + ".\nAvailable flux types: Uniform, PLAW, SEP, Galactic")
                    .
                    c_str());
    }

    if (fluxType == "Uniform") {
        flux = new UniformFlux(eTriggerThreshold);
    } else if (fluxType == "PLAW") {
        flux = new PLAWFlux(eTriggerThreshold);
    } else if (fluxType == "COMP") {
        flux = new COMPFlux(eTriggerThreshold);
    } else if (fluxType == "SEP") {
        flux = new SEPFlux(eTriggerThreshold);
    } else if (fluxType == "Galactic") {
        flux = new GalacticFlux(eTriggerThreshold);
    } else if (fluxType == "Table") {
        flux = new TableFlux(eTriggerThreshold);
    }

    if (saveSeeds) {
        fs::create_directories(eventsDir);
    }
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
    delete flux;
}

void PrimaryGeneratorAction::GenerateOnSphere(G4ThreeVector& pos, G4ThreeVector& dir) const {
    G4double u = 0;
    if (fluxDirection == "isotropic") {
        u = 2.0 * G4UniformRand() - 1.0;
    } else if (fluxDirection == "isotropic_up") {
        u = G4UniformRand();
    } else if (fluxDirection == "isotropic_down") {
        u = -G4UniformRand();
    }
    const G4double phi = 2.0 * M_PI * G4UniformRand();
    const G4double l = std::sqrt(std::max(0.0, 1.0 - u * u));
    const G4ThreeVector rhat(l * std::cos(phi), l * std::sin(phi), u);

    pos = center + radius * rhat;

    const G4ThreeVector z = rhat.unit();
    const G4ThreeVector a = std::fabs(z.z()) < 0.999 ? G4ThreeVector(0, 0, 1) : G4ThreeVector(1, 0, 0);
    const G4ThreeVector x = z.cross(a).unit();
    const G4ThreeVector y = z.cross(x).unit();

    const G4double ksi = G4UniformRand();
    const G4double sinTh = std::sqrt(ksi);
    const G4double cosTh = std::sqrt(1.0 - ksi);
    const G4double phi2 = 2.0 * M_PI * G4UniformRand();

    const G4ThreeVector v_local(sinTh * std::cos(phi2), sinTh * std::sin(phi2), cosTh);

    dir = -(v_local.x() * x + v_local.y() * y + v_local.z() * z);
    dir = dir.unit();
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* evt) {
    G4ThreeVector x, v;

    if (fluxDirection == "vertical_up") {
        v = G4ThreeVector(0., 0., 1.);
        const G4double halfSize = 2 * mm;   
        const G4double height   = -10.0 * cm;  
        const G4double x_ = (2 * G4UniformRand() - 1) * halfSize;
        const G4double y_ = (2 * G4UniformRand() - 1) * halfSize;
        const G4double z_ = -height;
        x = G4ThreeVector(x_, y_, z_);
    } else if (fluxDirection == "vertical_down") {
        v = G4ThreeVector (0., 0., -1.0);
        const G4double halfSize = 2 * mm;
        const G4double height = 10.0 * cm;  
        const G4double x_ = (2 * G4UniformRand() - 1) * halfSize;
        const G4double y_ = (2 * G4UniformRand() - 1) * halfSize;
        const G4double z_ = height;
        x = G4ThreeVector(x_, y_, z_);
    } else if (fluxDirection == "horizontal") {
        v = G4ThreeVector(-1., 0., 0.);
        const G4double halfSize = 2.0 * mm;
        const G4double y_ = (2 * G4UniformRand() - 1) * halfSize;
        const G4double z_ = (2 * G4UniformRand() - 1) * halfSize + 39.75;
        const G4double x_ = 2 * Sizes::Instrument::halfX;
        
        x = G4ThreeVector(x_, y_, z_);
    } else {
        GenerateOnSphere(x, v);
    }
    
    ParticleInfo info = flux->GenerateParticle();

    if (saveSeeds) {
        G4int eventID = evt->GetEventID();
        G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
        
        std::ostringstream fileNameStream;
        fileNameStream << eventsDir << "run" << runID << "_evt" << eventID << ".rndm";
        std::string fileName = fileNameStream.str();
        
        G4Random::saveEngineStatus(fileName.c_str());
        
        std::ofstream file(fileName, std::ios::app);
        file << "\n# ===== METADATA =====\n";
        file << "# fluxDirection: " << std::string(fluxDirection) << "\n";
        file << "# fluxType: " << std::string(fluxType) << "\n";
        file << "# particle: " << std::string(info.name) << "\n";
        file << "# pdg: " << info.pdg << "\n";
        file << "# energy_MeV: " << info.energy / MeV << "\n";
        file << "# pos_x_mm: " << x.x() / mm << "\n";
        file << "# pos_y_mm: " << x.y() / mm << "\n";
        file << "# pos_z_mm: " << x.z() / mm << "\n";
        file << "# dir_x: " << v.x() << "\n";
        file << "# dir_y: " << v.y() << "\n";
        file << "# dir_z: " << v.z() << "\n";
        file << "# TriggerThreshold_MeV: " << eTriggerThreshold / MeV << "\n";
        file << "# detectorType: " << std::string(Configuration::detectorType) << "\n";
        file << "# useOptics: " << (Configuration::useOptics ? "true" : "false") << "\n";
        file << "# ===== END METADATA =====\n";
        file.close();
    }

    particleGun->SetParticleDefinition(info.def);
    particleGun->SetParticleEnergy(info.energy);
    particleGun->SetParticlePosition(x);
    particleGun->SetParticleMomentumDirection(v);
    particleGun->SetParticleTime(0.0 * ns);
    particleGun->GeneratePrimaryVertex(evt);

    if (auto* ea = dynamic_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction())) {
        PrimaryRec rec;
        rec.index = static_cast<int>(ea->primBuf.size());
        rec.pdg = info.pdg;
        rec.name = info.name;
        rec.E_MeV = info.energy / MeV;
        rec.dir = v;
        rec.pos_mm = x / mm;
        rec.t0_ns = 0.0;
        ea->primBuf.emplace_back(std::move(rec));
    }
}