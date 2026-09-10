#ifndef FLUX_HH
#define FLUX_HH

#include <CLHEP/Units/SystemOfUnits.h>
#include "G4ParticleTable.hh"
#include <G4Types.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <unordered_map>

struct ParticleInfo {
    G4String name;
    G4int pdg;
    G4ParticleDefinition *def;
    G4double energy;
};


class Flux {
public:
    virtual ~Flux() = default;

    virtual ParticleInfo GenerateParticle();

protected:
    G4String particle;
    G4String configFile;
    G4double Emin{};
    G4double Emax{};

    virtual G4double SampleEnergy() = 0;

    static G4String Trim(const G4String &);

    G4String GetParam(const G4String &,
                      const G4String &,
                      const G4String &);
    G4double GetParam(const G4String &,
                      const G4String &,
                      G4double);

private:
    std::unordered_map<std::string, std::string> cache{};

    void LoadFileIfNeeded(const G4String &);
};


#endif //FLUX_HH
