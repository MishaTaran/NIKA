#ifndef MAGNETICFIELD_HH
#define MAGNETICFIELD_HH

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

class MagneticField {
public:
    MagneticField();
    ~MagneticField() = default;

    [[nodiscard]] G4FieldManager* GetFieldManager() const;

    MagneticField(const MagneticField&) = delete;
    MagneticField& operator=(const MagneticField&) = delete;

private:
    G4UniformMagField* field = nullptr;
    G4Mag_UsualEqRhs* equation = nullptr;
    G4MagIntegratorStepper* stepper = nullptr;
    G4MagInt_Driver* driver = nullptr;
    G4ChordFinder* chordFinder = nullptr;
    G4FieldManager* fieldManager = nullptr;
};

#endif
