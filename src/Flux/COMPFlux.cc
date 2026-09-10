#include "Flux/COMPFlux.hh"

COMPFlux::COMPFlux(const G4double cThreshold) {
    particle = "gamma";

    configFile = "../Flux_config/COMP_params.txt";
    alpha = GetParam(configFile, "alpha", 1.18511);
    E_Peak = GetParam(configFile, "E_Peak", 1.809619) * MeV;

    Emin = std::max({GetParam(configFile, "E_min", 0.01) * MeV, cThreshold});
    Emax = GetParam(configFile, "E_max", 50.) * MeV;

    BuildCDF();
}

void COMPFlux::BuildCDF() {
    const double k = (2.0 - alpha) / E_Peak;
    if (k <= 0.0) {
        G4Exception("COMPFlux::BuildCDF", "BAD_PARAM",
                    FatalException, "k <= 0 (alpha >= 2, распределение не определено)");
    }

    auto pdf = [&](double E) {
        return std::pow(E, -alpha) * std::exp(-k * E);
    };

    const int N = 2000;
    energyGrid.resize(N);
    cdfGrid.resize(N);

    double logEmin = std::log(Emin);
    double logEmax = std::log(Emax);

    std::vector<double> integrand(N);
    for (int i = 0; i < N; ++i) {
        double frac = double(i) / (N - 1);
        double Ei = std::exp(logEmin + frac * (logEmax - logEmin));
        energyGrid[i] = Ei;
        integrand[i] = pdf(Ei);
    }

    std::vector<double> cum(N, 0.0);
    for (int i = 1; i < N; ++i) {
        double dE = energyGrid[i] - energyGrid[i - 1];
        cum[i] = cum[i - 1] + 0.5 * (integrand[i] + integrand[i - 1]) * dE;
    }

    double Z = cum.back();
    if (!(Z > 0.0)) {
        G4Exception("COMPFlux::BuildCDF", "BAD_NORM",
                    FatalException, "Нормировка не посчиталась (Z <= 0)");
    }

    for (int i = 0; i < N; ++i) {
        cdfGrid[i] = cum[i] / Z;
    }

    cdfGrid.front() = 0.0;
    cdfGrid.back() = 1.0;
}


double COMPFlux::SampleEnergy() {
    double u = G4UniformRand();

    const auto it = std::lower_bound(cdfGrid.begin(), cdfGrid.end(), u);
    if (it == cdfGrid.begin()) return energyGrid.front();
    if (it == cdfGrid.end()) return energyGrid.back();

    const int idx = std::distance(cdfGrid.begin(), it);
    const double u1 = cdfGrid[idx - 1];
    const double u2 = cdfGrid[idx];
    const double E1 = energyGrid[idx - 1];
    double E2 = energyGrid[idx];

    const double t = (u - u1) / (u2 - u1);
    return E1 + t * (E2 - E1);
}
