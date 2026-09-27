#include "ActionInitialization.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"
#include "Configuration.hh"

using namespace Configuration;

ActionInitialization::ActionInitialization(const G4double a, const G4double Emin,
                                           const G4double Emax, bool saveSeeds)
    : EminMeV(Emin),
      EmaxMeV(Emax),
      area(a),
      saveSeeds(saveSeeds),
      fluxDirection(Configuration::fluxDirection),
      fluxType(Configuration::fluxType),
      useProvidedParams(false) {
    if (eTriggerThreshold > EmaxMeV) {
        G4Exception("ActionInitialization", "EnergyRange", FatalException,
                    "The energy threshold for a Trigger must be less than the maximum value in a given energy range");
    }
}

ActionInitialization::ActionInitialization(const G4double a, const G4double Emin,
                                           const G4double Emax,
                                           const G4String& fluxDir,
                                           const G4String& fluxType,
                                           bool saveSeeds)
    : EminMeV(Emin),
      EmaxMeV(Emax),
      area(a),
      saveSeeds(saveSeeds),
      fluxDirection(fluxDir),
      fluxType(fluxType),
      useProvidedParams(true) {
    if (eTriggerThreshold > EmaxMeV) {
        G4Exception("ActionInitialization", "EnergyRange", FatalException,
                    "The energy threshold for a Trigger must be less than the maximum value in a given energy range");
    }
}

void ActionInitialization::BuildForMaster() const {
    RunAction* runAct = new RunAction(area, EminMeV, EmaxMeV);
    SetUserAction(runAct);
}

void ActionInitialization::Build() const {
    RunAction* runAct = new RunAction(area, EminMeV, EmaxMeV);
    SetUserAction(runAct);

    EventAction* eventAct = new EventAction(runAct->analysisManager, runAct);
    SetUserAction(eventAct);

    PrimaryGeneratorAction* primaryGenerator;

    if (useProvidedParams) {
        primaryGenerator = new PrimaryGeneratorAction(fluxDirection, fluxType,
                                                      eTriggerThreshold, saveSeeds);
    } else {
        primaryGenerator = new PrimaryGeneratorAction(Configuration::fluxDirection,
                                                      Configuration::fluxType,
                                                      eTriggerThreshold,
                                                      saveSeeds);
    }
    SetUserAction(primaryGenerator);

    SteppingAction* stepAct = new SteppingAction();
    stepAct->EnableDebug(false);
    stepAct->SetDebugEventID(-1);
    stepAct->EnableGammaLog(true);   // <-- включить логирование гамма-взаимодействий
    stepAct->SetGammaLogEventID(-1); // -1 = все события; можно указать конкретный eventID
    SetUserAction(stepAct);
}