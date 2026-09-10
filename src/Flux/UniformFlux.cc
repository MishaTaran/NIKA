#include "Flux/UniformFlux.hh"


UniformFlux::UniformFlux(G4double cThreshold) : eTriggerThreshold(cThreshold) {
    const G4String filepath = "../Flux_config/Uniform_params.txt";

    const G4String particleLine = GetParam(filepath, "particles", "");
    const G4String fracLine = GetParam(filepath, "fractions", "");
    const G4String EminLine = GetParam(filepath, "E_min", "");
    const G4String EmaxLine = GetParam(filepath, "E_max", "");

    particles = Split(particleLine);
    fractions = ParseDoubles(fracLine);
    EminVec = ParseDoubles(EminLine);
    EmaxVec = ParseDoubles(EmaxLine);

    const size_t n = particles.size();

    if (fractions.size() != n || EminVec.size() != n || EmaxVec.size() != n) {
        G4Exception("UniformFlux::UniformFlux", "BAD_CONFIG",
                    JustWarning, "Configuration lists have mismatched lengths. Using equal fractions.");
        fractions.assign(n, 1.0 / std::max<size_t>(1, n));
    }

    double sum = std::accumulate(fractions.begin(), fractions.end(), 0.0);
    if (sum <= 0) sum = 1.0;
    for (auto &f: fractions) f /= sum;
}


std::vector<G4String> UniformFlux::Split(const G4String &line) {
    std::vector<G4String> result;
    std::stringstream ss(line);
    G4String token;
    while (std::getline(ss, token, ',')) {
        token = Trim(token);
        if (!token.empty())
            result.push_back(token);
    }
    return result;
}


size_t UniformFlux::SampleIndex() const {
    const double r = G4UniformRand();
    double cumulative = 0.0;
    for (size_t i = 0; i < fractions.size(); ++i) {
        cumulative += fractions[i];
        if (r <= cumulative) return i;
    }
    return fractions.size() - 1;
}


std::vector<G4double> UniformFlux::ParseDoubles(const G4String &line) {
    std::vector<G4double> result;
    std::stringstream ss(line);
    G4String token;
    while (std::getline(ss, token, ',')) {
        token = Trim(token);
        if (!token.empty()) {
            try {
                result.push_back(std::stod(token));
            } catch (...) {
                result.push_back(0.0);
            }
        }
    }
    return result;
}

G4double UniformFlux::SampleEnergy() {
    const size_t idx = SampleIndex();
    Emin = std::max({EminVec[idx] * MeV, eTriggerThreshold});
    Emax = EmaxVec[idx] * MeV;
    particle = particles[idx];

    return Emin * std::pow(Emax / Emin, G4UniformRand());;
}