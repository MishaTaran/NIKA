#include "CountRates.hh"


double fluxPLAW(const double E, const double A, double const alpha, const double E_piv) {
    return A * std::pow(E / E_piv, -alpha);
}

double fluxCOMP(const double E, const double A, const double alpha, const double E_piv, const double E_peak) {
    return A * std::pow(E / E_piv, -alpha) * std::exp((alpha - 2.0) * (E / E_peak));
}

// SEP: read CSV
static std::vector<double> readSepRow(int year, int order, const std::string& csvPath) {
    std::ifstream in(csvPath);
    if (!in.is_open()) {
        throw std::runtime_error("SEP: cannot open coefficients CSV: " + csvPath);
    }
    std::string line;
    std::getline(in, line);

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> cols;
        while (std::getline(ss, cell, ',')) {
            try {
                cols.push_back(std::stod(cell));
            }
            catch (...) {
                cols.push_back(std::numeric_limits<double>::quiet_NaN());
            }
        }
        if (cols.size() < 3) continue;

        int y = static_cast<int>(cols[0]);
        int k = static_cast<int>(cols[1]);
        if (y == year && k == order) {
            cols.pop_back();
            return cols;
        }
    }
    throw std::runtime_error("SEP: row not found for year=" + std::to_string(year) +
                             " order=" + std::to_string(order));
}

double fluxSEP(const double E, const int year, const int order, const std::string& csvPath) {
    static int cached_year = 0, cached_order = 0;
    static std::string cached_path;
    static std::vector<double> coeffs;

    if (year != cached_year || order != cached_order || csvPath != cached_path || coeffs.empty()) {
        coeffs = readSepRow(year, order, csvPath);
        cached_year = year;
        cached_order = order;
        cached_path = csvPath;
    }

    if (coeffs.size() < 4) return 0.0;
    const double x = std::log10(E);
    double poly = 0.0;
    for (size_t i = 2; i < cached_order + 3; i++) {
        poly += coeffs[i] * std::pow(x, static_cast<int>(i - 2));
    }
    return std::pow(10.0, poly);
}

// --- Table ---
double fluxTable(const double E, const std::string& csvPath) {
    static std::string cached_path;
    static std::vector<double> cached_energies;
    static std::vector<double> cached_fluxes;

    if (csvPath != cached_path || cached_energies.empty()) {
        std::ifstream in(csvPath);
        if (!in.is_open()) {
            throw std::runtime_error("Cannot open CSV file: " + csvPath);
        }

        std::string line;
        std::vector<double> energies;
        std::vector<double> fluxes;

        while (std::getline(in, line)) {
            if (line.empty()) continue;

            std::stringstream ss(line);
            std::string energy_str, flux_str;
            double energy, flux;

            std::getline(ss, energy_str, ',');
            std::getline(ss, flux_str, ',');

            try {
                energy = std::stod(energy_str);
                flux = std::stod(flux_str);
            }
            catch (...) {
                continue;
            }

            energies.push_back(energy);
            fluxes.push_back(flux);
        }

        if (energies.empty() || fluxes.empty()) {
            throw std::runtime_error("No valid data found in the CSV.");
        }

        cached_energies = std::move(energies);
        cached_fluxes = std::move(fluxes);
        cached_path = csvPath;
    }

    for (size_t i = 0; i < cached_energies.size(); ++i) {
        if (std::abs(cached_energies[i] - E) < 1e-6) {
            return cached_fluxes[i];
        }
    }

    for (size_t i = 1; i < cached_energies.size(); ++i) {
        if (cached_energies[i] > E) {
            double E1 = cached_energies[i - 1];
            double flux1 = cached_fluxes[i - 1];
            double E2 = cached_energies[i];
            double flux2 = cached_fluxes[i];

            double flux = flux1 + (flux2 - flux1) * (E - E1) / (E2 - E1);
            return flux;
        }
    }

    throw std::runtime_error("Energy is out of range in the CSV file.");
}

// --- Uniform ---
double fluxUniform(const double E, const double E_min, const double E_max) {
    return 1.0 / (E * std::log(E_max / E_min));
}


// --- Galactic ---
double J_Proton(const double E_GeV) {
    constexpr double E0 = 1.0;
    constexpr double mp = 0.938272;

    const double ETot = E_GeV + mp;
    double beta_sq = 1.0 - mp * mp / (ETot * ETot);
    if (beta_sq <= 0.0) beta_sq = 1e-12;

    double term1 = 2620.0 / beta_sq * std::pow(E_GeV / E0, 1.1);
    term1 *= std::pow((std::pow(E_GeV / E0, 0.98) + std::pow(0.7, 0.98)) / (1.0 + std::pow(0.7, 0.98)), -4.0);

    const double term2 = 30.0 * std::pow(E_GeV / E0, 2.0) * std::pow((E_GeV / E0 + 8.0) / 9.0, -12.0);

    return term1 + term2;
}


double J_Electron(const double E) {
    constexpr double E0 = 1.0;
    constexpr double me = 0.000511;

    const double ETot = E + me;
    double beta_sq = 1.0 - me * me / (ETot * ETot);
    if (beta_sq <= 0) {
        beta_sq = 1e-6;
    }

    double term1 = 255.0 / beta_sq * std::pow(E / E0, -1.0);
    term1 *= std::pow((E / E0 + 0.63) / 1.63, -2.43);

    const double term2 = 6.4 * std::pow(E / E0, 2.0) * std::pow((E / E0 + 15.0) / 16.0, -26.0);

    return term1 + term2;
}

double J_Positron(const double E) {
    constexpr double E0 = 1.0;
    constexpr double mp = 0.000511;

    const double ETot = E + mp;
    double beta = std::sqrt(1.0 - mp * mp / (ETot * ETot));
    if (beta <= 0) beta = 1e-6;

    double term1 = 25.0 / (beta * beta) * std::pow(E / E0, 0.1);
    term1 *= std::pow((std::pow(E / E0, 1.1) + std::pow(0.2, 1.1)) / (1.0 + std::pow(0.2, 1.1)), -3.31);

    const double term2 = 23.0 * std::pow(E / E0, 0.5) * std::pow((E / E0 + 2.2) / 3.2, -9.5);

    return term1 + term2;
}

double J_Alpha(const double E) {
    constexpr double E0 = 1.0;
    constexpr double malpha = 3.727379;

    const double ETot = E + malpha;
    double beta_sq = 1.0 - malpha * malpha / (ETot * ETot);
    if (beta_sq <= 0) beta_sq = 1e-6;

    double term1 = 163.4 / beta_sq * std::pow(E / E0, 1.1);
    term1 *= std::pow((std::pow(E / E0, 0.97) + std::pow(0.58, 0.97)) / (1.0 + std::pow(0.58, 0.97)), -4.0);

    return term1;
}

double fluxGalactic(const double E, const double phiMV, const std::string& name) {
    double mass = 0;
    const int Z = name == "alpha" ? 2 : 1;
    const double phiGV = phiMV * 1e-3 * Z;

    const double ELis = E + phiGV;
    if (ELis <= 0) return 0.0;
    if (name == "alpha") {
        mass = 3.727379;
    } else if (name == "e-" or name == "e+") {
        mass = 0.000511;
    } else if (name == "proton") {
        mass = 0.938272;
    }
    const double num = E * (E + 2.0 * mass);
    const double den = ELis * (ELis + 2.0 * mass);
    if (den <= 0) {
        return 0.0;
    }

    double J_LIS = 0;
    if (name == "alpha") {
        J_LIS = J_Alpha(ELis);
    } else if (name == "e-") {
        J_LIS = J_Electron(ELis);
    } else if (name == "proton") {
        J_LIS = J_Proton(ELis);
    } else if (name == "e+") {
        J_LIS = J_Positron(ELis);
    }
    return num / den * J_LIS;
}

// ---------------- Area ----------------

double AreaGen_cm2(const double halfX_mm, const double halfY_mm, const double halfZ_mm, const FluxDir dir) {
    const double to_cm2 = 0.01;  // mm^2 -> cm^2 (1 cm^2 = 100 mm^2)

    if (dir == FluxDir::Vertical_up) {
        return 4 * halfX_mm * halfY_mm * to_cm2;
    }
    if (dir == FluxDir::Vertical_down) {
        return 4 * halfX_mm * halfY_mm * to_cm2;
    }
    if (dir == FluxDir::Horizontal) {
        return 4 * halfY_mm * halfZ_mm * to_cm2;
    }
    double radius = sqrt(halfX_mm * halfX_mm + halfY_mm * halfY_mm + halfZ_mm * halfZ_mm) + 5;

    // Isotropic / Isotropic_up / Isotropic_down: sphere area 4*pi*r^2, projected = (1/4)*4*pi*r^2 = pi*r^2
    return M_PI * radius * radius * to_cm2;
}

double AreaRect_cm2(const double halfX_mm, const double halfY_mm, const double sizeZ_mm, const FluxDir dir) {
    const double to_cm2 = 0.01;  // mm^2 -> cm^2 (1 cm^2 = 100 mm^2)
    const double faceXY = 4.0 * halfX_mm * halfY_mm;
    const double faceXZ = 2.0 * halfX_mm * sizeZ_mm;
    const double faceYZ = 2.0 * halfY_mm * sizeZ_mm;

    if (dir == FluxDir::Vertical_up || dir == FluxDir::Vertical_down) {
        return faceXY * to_cm2;
    }
    if (dir == FluxDir::Horizontal) {
        return (faceXZ + faceYZ) * to_cm2;
    }
    // Isotropic: mean projected area = total surface / 4
    const double totalSurface = 2.0 * (faceXY + faceXZ + faceYZ);
    const double area = totalSurface / 4.0;
    return area * to_cm2;
}

// ---------------- Integral ----------------

static double simpson(const std::function<double(double)>& f, const double a, const double b) {
    const double c = 0.5 * (a + b);
    return (b - a) * (f(a) + 4.0 * f(c) + f(b)) / 6.0;
}

static double adaptiveSimpsonRec(const std::function<double(double)>& f,
                                 const double a, const double b, const double eps,
                                 const double whole, const int depth) {
    const double c = 0.5 * (a + b);
    const double left = simpson(f, a, c);
    const double right = simpson(f, c, b);
    const double delta = left + right - whole;
    if (depth <= 0 || std::fabs(delta) <= 15.0 * eps) {
        return left + right + delta / 15.0;
    }
    return adaptiveSimpsonRec(f, a, c, eps / 2.0, left, depth - 1) +
        adaptiveSimpsonRec(f, c, b, eps / 2.0, right, depth - 1);
}

double integrateAdaptiveSimpson(const std::function<double(double)>& f,
                                const double a, const double b,
                                const double rel_tol, const int max_depth) {
    if (a == b) return 0;
    const double initial = simpson(f, a, b);
    const double eps = rel_tol * std::max(1.0, std::fabs(initial));
    return adaptiveSimpsonRec(f, a, b, eps, initial, max_depth);
}

static inline double binEdgeLog(double Emin, double Emax, int nBins, int i) {
    // i in [0..nBins], logspace like np.logspace
    const double ratio = Emax / Emin;
    const double t = static_cast<double>(i) / static_cast<double>(nBins);
    return Emin * std::pow(ratio, t);
}


RateResult computeRate(const FluxType type,
                       const FluxParams& p,
                       EnergyRange eRange,
                       double A_eff_cm2,
                       const int N_histories,
                       const RateCounts& detCounts) {
    std::function<double(double)> f;

    switch (type) {
    case FluxType::PLAW:
        f = [=](const double E) {
            return fluxPLAW(E, p.A, p.alpha, p.E_piv);
        };
        break;
    case FluxType::COMP:
        f = [=](const double E) {
            return fluxCOMP(E, p.A, p.alpha, p.E_piv, p.E_peak);
        };
        break;
    case FluxType::SEP:
        f = [=](const double E) {
            return fluxSEP(E, p.sep_year, p.sep_order, p.sep_csv_path);
        };
        break;
    case FluxType::TABLE:
        f = [=](const double E) {
            return fluxTable(E, p.table_path);
        };
        break;
    case FluxType::UNIFORM:
        f = [=](const double E) {
            return fluxUniform(E, eRange.Emin, eRange.Emax);
        };
        break;
    case FluxType::GALACTIC: {
        eRange.Emin /= 1000.0;
        eRange.Emax /= 1000.0;
        A_eff_cm2 /= 10000.0;
        f = [=](const double E_GeV) {
            return fluxGalactic(E_GeV, p.phiMV, p.particle);
        };
        break;
    }
    default:
        throw std::runtime_error("Unknown flux type");
    }

    const double integral = integrateAdaptiveSimpson(f, eRange.Emin, eRange.Emax, 1e-6, 22);
    const double Ndot = A_eff_cm2 * integral;

    RateResult R;
    R.area = A_eff_cm2;
    R.integral = integral;
    R.Ndot = Ndot;
    R.rateTriggerNoVeto = N_histories > 0 ? static_cast<double>(detCounts.triggerNoVeto) * Ndot / N_histories : 0.0;
    const int bothDet = detCounts.triggerNoVeto + detCounts.triggerWithVeto;
    R.rateBoth = N_histories > 0 ? static_cast<double>(bothDet) * Ndot / N_histories : 0.0;
    return R;
}

RateResult computeRateReal(FluxType type,
                           const FluxParams& p,
                           EnergyRange eRange,
                           const std::vector<double>& Aeff,
                           int nBins) {
    if (nBins <= 0) throw std::runtime_error("computeRateReal: nBins <= 0");
    if (static_cast<int>(Aeff.size()) != nBins)
        throw std::runtime_error("computeRateReal: Aeff.size() != nBins");
    if (eRange.Emin <= 0.0 || eRange.Emax <= 0.0 || eRange.Emax <= eRange.Emin)
        throw std::runtime_error("computeRateReal: invalid energy range");

    std::function<double(double)> fluxF;


    double areaScale = 1.0;
    double energyScale = 1.0;

    switch (type) {
    case FluxType::PLAW:
        fluxF = [=](double E_MeV) {
            return fluxPLAW(E_MeV, p.A, p.alpha, p.E_piv);
        };
        break;
    case FluxType::COMP:
        fluxF = [=](double E_MeV) {
            return fluxCOMP(E_MeV, p.A, p.alpha, p.E_piv, p.E_peak);
        };
        break;
    case FluxType::SEP:
        fluxF = [=](double E_MeV) {
            return fluxSEP(E_MeV, p.sep_year, p.sep_order, p.sep_csv_path);
        };
        break;
    case FluxType::TABLE:
        fluxF = [=](double E_MeV) {
            return fluxTable(E_MeV, p.table_path);
        };
        break;
    case FluxType::UNIFORM:
        fluxF = [=](double E_MeV) {
            return fluxUniform(E_MeV, eRange.Emin, eRange.Emax);
        };
        break;
    case FluxType::GALACTIC:
        energyScale = 1.0 / 1000.0;
        areaScale = 1.0 / 10000.0;
        fluxF = [=](double E_GeV) {
            return fluxGalactic(E_GeV, p.phiMV, p.particle);
        };
        break;
    default:
        throw std::runtime_error("computeRateReal: unknown flux type");
    }

    double rateReal = 0.0;

    for (int i = 0; i < nBins; ++i) {
        const double e1 = binEdgeLog(eRange.Emin, eRange.Emax, nBins, i);
        const double e2 = binEdgeLog(eRange.Emin, eRange.Emax, nBins, i + 1);
        const double Ec = std::sqrt(e1 * e2);
        const double dE = e2 - e1;

        double A = Aeff[i];

        const double Earg = Ec * energyScale;
        const double dEarg = dE * energyScale;
        const double Aarg = A * areaScale;

        const double phi = fluxF(Earg);
        rateReal += phi * Aarg * dEarg;
    }

    RateResult R;
    R.rateRealTrigger = rateReal;
    return R;
}
