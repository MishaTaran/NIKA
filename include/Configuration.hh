#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <G4String.hh>
#include <G4Types.hh>
#include <G4SystemOfUnits.hh>

namespace Configuration
{
    inline G4String detectorType{"CsI"};

    inline G4String fluxType{"Uniform"};
    inline G4String fluxDirection{"isotropic"};

    inline G4double eTriggerThreshold{0 * MeV};
    inline G4double eVetoThreshold{0 * MeV};

    constexpr G4double kPlaneThreshold  = 0.15;
    constexpr G4double kVetoThreshold   = 0.2;

    inline G4bool useOptics{false};
    inline G4int yieldScale{1};
    inline G4int oTriggerThreshold{0};
    inline G4int oSideVetoThreshold{0};
    inline G4int oBottomVetoThreshold{0};
    inline G4int oUpperVetoThreshold{0};

    inline G4bool polishedTyvek{false};
    inline G4double viewDeg{360 * deg};

    inline G4int nBins{1000};
    inline G4String outputFile{"GammaCube.root"};
    inline G4bool saveSecondaries{true};
    inline G4bool savePhotons{false};
}

#endif // CONFIGURATION_HH