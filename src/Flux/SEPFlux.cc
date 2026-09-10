#include "Flux/SEPFlux.hh"


SEPFlux::SEPFlux(const G4double cThreshold) {
    path = "../SEP_spectrum.CSV";
    particle = "proton";

    configFile = "../Flux_config/SEP_params.txt";
    year = static_cast<int>(GetParam(configFile, "year", 1998));
    order = static_cast<int>(GetParam(configFile, "order", 15));

    Emin = std::max({GetParam(configFile, "E_min", 0.1) * MeV, cThreshold});
    Emax = GetParam(configFile, "E_max", 1000.) * MeV;

    BuildCDF();
}


static std::vector<double> ExtractNumbers(const std::string &line) {
    static const std::regex re(R"(([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?))");
    std::vector<double> out;
    for (std::sregex_iterator it(line.begin(), line.end(), re), end; it != end; ++it) {
        out.push_back(std::stod((*it)[1].str()));
    }
    return out;
}


void SEPFlux::BuildCDF() {
    EList.clear();
    CDF.clear();

    std::ifstream in(path.c_str());
    if (!in) {
        G4Exception("SEPFlux::BuildCDF", "CSV_OPEN_FAIL",
                    JustWarning, ("Cannot open " + path).c_str());
        EList = {1. * MeV, 10. * MeV};
        CDF = {0.0, 1.0};
        return;
    }

    std::vector<Row> rows;
    rows.reserve(1884);

    std::string line;
    bool headerSkipped = false;
    while (std::getline(in, line)) {
        if (!headerSkipped) {
            headerSkipped = true;
            continue;
        }
        auto nums = ExtractNumbers(line);
        if (nums.size() < 4) continue;

        const int yr = static_cast<int>(std::llround(nums[0]));
        const double E_csv = nums[1];
        const double flx = nums[2];
        const int ord = static_cast<int>(std::llround(nums.back()));

        if (yr == year && ord == order) {
            const double E_G4 = E_csv * MeV;
            rows.push_back({E_G4, flx});
        }
    }
    in.close();

    if (rows.size() < 2) {
        G4Exception("PrimaryGeneratorAction::BuildCSVFluxCDF", "CSV_NO_ROWS",
                    JustWarning, "No matching rows for given year/order.");
        EList = {1. * MeV, 10. * MeV};
        CDF = {0.0, 1.0};
        return;
    }

    std::sort(rows.begin(), rows.end(),
              [](const Row &a, const Row &b) { return a.E_MeV < b.E_MeV; });

    {
        std::vector<Row> uniq;
        uniq.reserve(rows.size());
        for (const auto &r: rows) {
            if (!uniq.empty() && std::abs(r.E_MeV - uniq.back().E_MeV) < 1e-12 * MeV) continue;
            uniq.push_back(r);
        }
        rows.swap(uniq);
    }

    const size_t Ntot = rows.size();
    std::vector<double> E(Ntot), lE(Ntot), f_perE(Ntot, 0.0);
    for (size_t i = 0; i < Ntot; ++i) {
        E[i] = rows[i].E_MeV;
        lE[i] = std::log(E[i]);
    }

    auto clamp_pos = [](const double x) { return x > 0.0 && std::isfinite(x) ? x : 0.0; };

    for (size_t i = 0; i < Ntot; ++i) {
        f_perE[i] = clamp_pos(rows[i].flux);
    }

    bool useBounds = Emin > 0.0 && Emax > 0.0 && std::isfinite(Emin) && std::isfinite(Emax) && Emax > Emin;

    size_t iLo = 0, iHi = Ntot - 1;
    if (useBounds) {
        auto itLo = std::lower_bound(E.begin(), E.end(), Emin);
        if (itLo == E.begin()) {
            iLo = 0;
        } else if (itLo == E.end()) {
            iLo = Ntot - 1;
        } else {
            size_t j = static_cast<size_t>(itLo - E.begin());
            iLo = Emin - E[j - 1] <= E[j] - Emin ? j - 1 : j;
        }

        auto itHi = std::upper_bound(E.begin(), E.end(), Emax);
        if (itHi == E.begin()) {
            iHi = 0;
        } else if (itHi == E.end()) {
            iHi = Ntot - 1;
        } else {
            size_t j = static_cast<size_t>(itHi - E.begin());
            size_t jm1 = j - 1;
            if (j < Ntot) {
                iHi = Emax - E[jm1] <= E[j] - Emax ? jm1 : j;
            } else {
                iHi = jm1;
            }
        }

        if (iLo > iHi) std::swap(iLo, iHi);
        if (iHi == iLo && Ntot >= 2) {
            if (iHi + 1 < Ntot) ++iHi;
            else if (iLo > 0) --iLo;
        }
    }

    const size_t N = iHi - iLo + 1;
    if (N < 2) {
        G4Exception("PrimaryGeneratorAction::BuildCSVFluxCDF", "CSV_RANGE_TOO_NARROW",
                    JustWarning, "Energy range too narrow (fewer than 2 points).");
        EList = {1. * MeV, 10. * MeV};
        CDF = {0.0, 1.0};
        return;
    }

    std::vector<double> Es(N), lEs(N), f_sub(N);
    for (size_t k = 0; k < N; ++k) {
        Es[k] = E[iLo + k];
        lEs[k] = std::log(Es[k]);
        f_sub[k] = f_perE[iLo + k];
    }

    EList = Es;
    CDF.assign(N, 0.0);
    long double acc = 0.0L;
    for (size_t i = 0; i + 1 < N; ++i) {
        const long double gi = static_cast<long double>(f_sub[i]) * static_cast<long double>(Es[i]);
        const long double gip1 = static_cast<long double>(f_sub[i + 1]) * static_cast<long double>(Es[i + 1]);
        const long double dlnE = lEs[i + 1] - lEs[i];
        const long double wi = 0.5L * (gi + gip1) * dlnE;
        acc += wi;
        CDF[i + 1] = static_cast<double>(acc);
    }

    if (acc <= 0.0L || !std::isfinite(static_cast<double>(acc))) {
        for (size_t i = 0; i < N; ++i) CDF[i] = static_cast<double>(i) / static_cast<double>(N - 1);
        return;
    }

    for (auto &v: CDF) v /= static_cast<double>(acc);
    CDF.front() = 0.0;
    CDF.back() = 1.0;
}

G4double SEPFlux::SampleEnergy() {
    if (EList.size() < 2) return 1.0 * MeV;

    const G4double u = G4UniformRand();

    const auto it = std::upper_bound(CDF.begin(), CDF.end(), u);
    if (it == CDF.begin()) return EList.front();
    if (it == CDF.end()) return EList.back();
    const size_t j = std::distance(CDF.begin(), it);
    const size_t i = j - 1;

    const double C0 = CDF[i];
    const double C1 = CDF[j];
    const double t = (u - C0) / std::max(C1 - C0, 1e-12);

    const double lnE0 = std::log(EList[i]);
    const double lnE1 = std::log(EList[j]);
    const double lnEu = lnE0 + t * (lnE1 - lnE0);

    return std::exp(lnEu);
}
