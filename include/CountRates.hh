#ifndef COUNTRATES_HH
#define COUNTRATES_HH

#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>

enum class FluxType { PLAW, COMP, SEP, UNIFORM, GALACTIC, TABLE };

struct EnergyRange {
    double Emin;
    double Emax;
};

struct FluxParams {
    // PLAW / COMP
    double A = 0.0;
    double alpha = 0.0;
    double E_piv = 1.0;
    double E_peak = 1.0; // only COMP

    // SEP
    int sep_year = 0;
    int sep_order = 0;
    std::string sep_csv_path; // path to CSV with coefficients

    // Galactic
    double phiMV = 600.0;
    std::string particle = "proton";

    // Table
    std::string table_path;
};

struct RateCounts {
    int triggerNoVeto = 0;    // N_det (Trigger && !Veto)
    int triggerWithVeto = 0;  // N_det (Trigger && Veto)
};

struct RateResult {
    double area = 0.0;
    double integral = 0.0;           // ∫ flux(E) dE
    double Ndot = 0.0;               // A_eff * integral
    double rateTriggerNoVeto = 0.0;   // triggerNoVeto / (N / Ndot)
    double rateBoth = 0.0;            // (triggerNoVeto+triggerWithVeto) / (N / Ndot)
    double rateRealTrigger = 0.0;     // ∫ flux(E) * Aeff(E) dE
};

double fluxPLAW(double E, double A, double alpha, double E_piv);

double fluxCOMP(double E, double A, double alpha, double E_piv, double E_peak);

double fluxSEP(double E, int year, int order, const std::string& csvPath);

double fluxTable(double E, const std::string& csvPath);

double fluxUniform(double E, double E_min, double E_max);

double fluxGalactic(double E, double phiMV, const std::string& particle);

double J_Proton(double E_GeV);
double J_Electron(double E_GeV);
double J_Positron(double E_GeV);
double J_Alpha(double E_GeV);

enum class FluxDir { Vertical_down, Vertical_up, Horizontal, Isotropic_up, Isotropic_down, Isotropic };

/** Area in cm^2 for a rectangular envelope: halfX_mm, halfY_mm (half extents), sizeZ_mm (full height). */
double AreaRect_cm2(double halfX_mm, double halfY_mm, double sizeZ_mm, FluxDir dir);

/** Area in cm^2 of the surface from which primaries are launched (same geometry as PrimaryGeneratorAction).
 *  halfY_mm = max(halfX, halfY) of envelope, sizeZ_mm = envelope full height,
 *  radiusVerticalDown_mm = circle radius for vertical_down, radiusSphere_mm = sphere radius for isotropic. */
double AreaGen_cm2(double halfX_mm, double halfY_mm, double halfZ_mm, FluxDir dir);


double integrateAdaptiveSimpson(const std::function<double(double)>& f,
                                double a, double b,
                                double rel_tol = 1e-6, int max_depth = 20);


RateResult computeRate(FluxType type,
                       const FluxParams& p,
                       EnergyRange eRange,
                       double A_eff_cm2,
                       int N_histories,
                       const RateCounts& detCounts);

RateResult computeRateReal(FluxType type,
                           const FluxParams& p,
                           EnergyRange eRange,
                           const std::vector<double>& Aeff,
                           int nBins);


#endif //COUNTRATES_HH
