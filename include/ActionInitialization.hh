#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include <G4VUserActionInitialization.hh>
#include <G4Types.hh>
#include <G4String.hh>

class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization(G4double a, G4double Emin, G4double Emax, bool saveSeeds = false);
    ActionInitialization(G4double a, G4double Emin, G4double Emax, 
                         const G4String& fluxDir, const G4String& fluxType, 
                         bool saveSeeds = false);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    G4double EminMeV;
    G4double EmaxMeV;
    G4double area;
    bool saveSeeds;
    
    G4String fluxDirection;
    G4String fluxType;
    bool useProvidedParams = false;
};

#endif