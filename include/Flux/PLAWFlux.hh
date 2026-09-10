#ifndef PLAWFLUX_HH
#define PLAWFLUX_HH

#include "Flux.hh"

class PLAWFlux : public Flux {
public:
    explicit PLAWFlux(G4double cThreshold);

private:
    G4double alpha{};

    G4double SampleEnergy() override;
};


#endif //PLAWFLUX_HH
