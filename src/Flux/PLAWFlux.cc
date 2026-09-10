#include "Flux/PLAWFlux.hh"


PLAWFlux::PLAWFlux(const G4double cThreshold) {
    particle = "gamma";

    configFile = "../Flux_config/PLAW_params.txt";
    alpha = GetParam(configFile, "alpha", 1.411103);

    Emin = std::max({GetParam(configFile, "E_min", 0.01) * MeV, cThreshold});
    Emax = GetParam(configFile, "E_max", 100.) * MeV;
}


double PLAWFlux::SampleEnergy() {
    double u = G4UniformRand();
    if (std::abs(alpha - 1.0) < 1e-12) {
        return Emin * std::pow(Emax / Emin, u);
    }
    double EminPow = std::pow(Emin, 1.0 - alpha);
    double EmaxPow = std::pow(Emax, 1.0 - alpha);
    double val = EminPow + u * (EmaxPow - EminPow);
    return std::pow(val, 1.0 / (1.0 - alpha));
}
