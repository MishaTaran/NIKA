#include "Flux/GalacticFlux.hh"

GalacticFlux::GalacticFlux(const G4double cThreshold) {

    configFile = "../Flux_config/Galactic_params.txt";
    particle = GetParam(configFile, "particle", "proton");
    phiMV = GetParam(configFile, "phiMV", 600);

    Emin = std::max({GetParam(configFile, "E_min", 1.) * MeV, cThreshold});
    Emax = GetParam(configFile, "E_max", 1000000.) * MeV;

    BuildCDF();
}


void GalacticFlux::BuildCDF() {
    constexpr G4int NBins = 1000;
    energyGrid.resize(NBins);
    cdfGrid.resize(NBins);

    const G4double logEmin = std::log(Emin);
    const G4double logEmax = std::log(Emax);

    for (int i = 0; i < NBins; i++) {
        const G4double lg = logEmin + (logEmax - logEmin) * i / (NBins - 1);
        energyGrid[i] = std::exp(lg);
    }

    G4double integral = 0.0;
    cdfGrid[0] = 0.0;
    for (int i = 1; i < NBins; i++) {
        const G4double E1 = energyGrid[i - 1];
        const G4double E2 = energyGrid[i];
        const G4double f1 = J_TOA_GeV(E1);
        const G4double f2 = J_TOA_GeV(E2);
        integral += 0.5 * (f1 + f2) * (E2 - E1);
        cdfGrid[i] = integral;
    }

    if (integral <= 0.0 || !std::isfinite(integral)) {
        G4Exception("GalacticFlux::BuildCDF", "BAD_INTEGRAL", FatalException,
                    "Integral of spectrum is non-positive or non-finite.");
    }
    for (int i = 0; i < NBins; i++) {
        cdfGrid[i] /= integral;
    }
    cdfGrid.back() = 1.0;
}


G4double GalacticFlux::J_Proton(const G4double E) const {
    constexpr G4double E0 = 1.0;
    constexpr G4double mp = 0.938272;

    const G4double ETot = E + mp;
    G4double beta = std::sqrt(1.0 - mp * mp / (ETot * ETot));
    if (beta <= 0) {
        beta = 1e-6;
    }
    G4double term1 = 2620.0 / (beta * beta) * std::pow(E / E0, 1.1);
    term1 *= std::pow((std::pow(E / E0, 0.98) + std::pow(0.7, 0.98)) / (1.0 + std::pow(0.7, 0.98)), -4.0);

    const G4double term2 = 30.0 * std::pow(E / E0, 2.0) * std::pow((E / E0 + 8.0) / 9.0, -12.0);

    return term1 + term2;
}


G4double GalacticFlux::J_Electron(const G4double E) const {
    constexpr G4double E0 = 1.0;
    constexpr G4double me = 0.000511;

    const G4double ETot = E + me;
    G4double beta_sq = 1.0 - me * me / (ETot * ETot);
    if (beta_sq <= 0) {
        beta_sq = 1e-6;
    }

    G4double term1 = 255.0 / beta_sq * std::pow(E / E0, -1.0);
    term1 *= std::pow((E / E0 + 0.63) / 1.63, -2.43);

    const G4double term2 = 6.4 * std::pow(E / E0, 2.0) * std::pow((E / E0 + 15.0) / 16.0, -26.0);

    return term1 + term2;
}


G4double GalacticFlux::J_Positron(const G4double E) const {
    constexpr G4double E0 = 1.0;
    constexpr G4double mp = 0.000511;

    const G4double ETot = E + mp;
    G4double beta = std::sqrt(1.0 - mp * mp / (ETot * ETot));
    if (beta <= 0) beta = 1e-6;

    G4double term1 = 25.0 / (beta * beta) * std::pow(E / E0, 0.1);
    term1 *= std::pow((std::pow(E / E0, 1.1) + std::pow(0.2, 1.1)) / (1.0 + std::pow(0.2, 1.1)), -3.31);

    const G4double term2 = 23.0 * std::pow(E / E0, 0.5) * std::pow((E / E0 + 2.2) / 3.2, -9.5);

    return term1 + term2;
}


G4double GalacticFlux::J_Alpha(const G4double E) const {
    constexpr G4double E0 = 1.0;
    constexpr G4double malpha = 3.727379;

    const G4double ETot = E + malpha;
    G4double beta_sq = 1.0 - malpha * malpha / (ETot * ETot);
    if (beta_sq <= 0) beta_sq = 1e-6;

    G4double term1 = 163.4 / beta_sq * std::pow(E / E0, 1.1);
    term1 *= std::pow((std::pow(E / E0, 0.97) + std::pow(0.58, 0.97)) / (1.0 + std::pow(0.58, 0.97)), -4.0);

    return term1;
}

G4double GalacticFlux::J_TOA_GeV(const G4double E) {
    constexpr G4double mp = 0.938272;
    constexpr G4double me = 0.000511;
    constexpr G4double malpha = 3.727379;
    G4double mass = 0;
    const G4double Z = particle == "alpha" ? 0.5 : 1;
    const G4double phiGeV = phiMV * 1e-3 * Z;

    const G4double ELis = E + phiGeV;
    if (ELis <= 0) return 0.0;
    if (particle == "alpha") {
        mass = malpha;
    } else if (particle == "e-" or particle == "e+") {
        mass = me;
    } else if (particle == "proton") {
        mass = mp;
    }
    const G4double num = E * (E + 2.0 * mass);
    const G4double den = ELis * (ELis + 2.0 * mass);
    if (den <= 0) {
        return 0.0;
    }

    G4double J_LIS = 0;
    if (particle == "alpha") {
        Emin = Emin * 4;
        Emax = Emax * 4;
        J_LIS = 4 * J_Alpha(ELis / 4);
    } else if (particle == "e-") {
        J_LIS = J_Electron(ELis);
    } else if (particle == "proton") {
        J_LIS = J_Proton(ELis);
    } else if (particle == "e+") {
        J_LIS = J_Positron(ELis);
    }
    return num / den * J_LIS;
}

G4double GalacticFlux::SampleEnergy() {
    const G4double u = G4UniformRand(); // равномерное [0,1)

    const auto it = std::lower_bound(cdfGrid.begin(), cdfGrid.end(), u);
    const int idx = std::max(1, static_cast<int>(it - cdfGrid.begin()));

    const G4double u1 = cdfGrid[idx - 1];
    const G4double u2 = cdfGrid[idx];
    const G4double E1 = energyGrid[idx - 1];
    const G4double E2 = energyGrid[idx];

    const G4double t = (u - u1) / (u2 - u1);
    return (E1 + t * (E2 - E1)) * GeV;
}
