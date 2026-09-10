#ifndef GALACTICFLUX_HH
#define GALACTICFLUX_HH

#include "Flux/Flux.hh"


class GalacticFlux : public Flux {
public:
    explicit GalacticFlux(G4double cThreshold);

private:
    G4double phiMV{};

    std::vector<G4double> energyGrid;
    std::vector<G4double> cdfGrid;

    void BuildCDF();

    [[nodiscard]] G4double J_Proton(G4double E) const;
    [[nodiscard]] G4double J_Electron(G4double E) const;
    [[nodiscard]] G4double J_Positron(G4double E) const;
    [[nodiscard]] G4double J_Alpha(G4double E) const;
    [[nodiscard]] G4double J_TOA_GeV(G4double E);

    G4double SampleEnergy() override;
};

#endif //GALACTICFLUX_HH
