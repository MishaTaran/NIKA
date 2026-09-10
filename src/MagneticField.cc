#include "MagneticField.hh"

#include "G4DormandPrince745.hh"

MagneticField::MagneticField() {
    const G4ThreeVector B(0., -0.2 * tesla, 0.);
    const G4double minStep = 0.1 * mm;

    const G4double deltaIntersection = 1e-3 * mm;
    const G4double deltaOneStep = 1e-4 * mm;
    const G4double deltaChord = 0.1 * mm;

    field = new G4UniformMagField(B);
    fieldManager = new G4FieldManager(field);

    equation = new G4Mag_UsualEqRhs(field);
    stepper = new G4DormandPrince745(equation);

    driver = new G4MagInt_Driver(minStep, stepper, stepper->GetNumberOfVariables());

    chordFinder = new G4ChordFinder(driver);

    fieldManager->SetDeltaIntersection(deltaIntersection);
    fieldManager->SetDeltaOneStep(deltaOneStep);
    fieldManager->SetChordFinder(chordFinder);

    chordFinder->SetDeltaChord(deltaChord);
}

G4FieldManager* MagneticField::GetFieldManager() const {
    return fieldManager;
}
