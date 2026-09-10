#ifndef RUNACTION_HH
#define RUNACTION_HH

#include <G4UserRunAction.hh>
#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <G4Accumulable.hh>
#include <G4AccumulableManager.hh>
#include <G4Run.hh>
#include <G4ios.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>
#include <G4AnalysisManager.hh>
#include <G4Threading.hh>
#include <iomanip>
#include <sstream>

#include "Sizes.hh"
#include "Configuration.hh"
#include "AnalysisManager.hh"

struct ParticleCounts {
    G4int triggerNoVeto = 0;
    G4int triggerWithVeto = 0;
};

class RunAction : public G4UserRunAction {
public:
    AnalysisManager *analysisManager;

    RunAction();
    RunAction(double Agen_cm2, double Emin_MeV, double Emax_MeV);
    ~RunAction() override;

    void BeginOfRunAction(const G4Run *) override;
    void EndOfRunAction(const G4Run *) override;

    // Существующие методы
    void AddTriggerNoVetoCount(const G4int v) { triggerNoVeto += v; }
    void AddTriggerWithVeto(const G4int v) { triggerWithVeto += v; }

    void AddTriggerNoVetoOpt(const G4int v) { triggerNoVetoOpt += v; }
    void AddTriggerWithVetoOpt(const G4int v) { triggerWithVetoOpt += v; }

    void AddGenerated(double E_MeV);
    void AddTriggerNoVeto(double E_MeV);
    void AddTriggerNoVetoOpt(double E_MeV);

    // НОВЫЕ МЕТОДЫ ДЛЯ КРИСТАЛЛА И ВЕТО
    void AddTriggerOnly(G4int v) { 
        // Увеличиваем счетчик событий "только кристалл"
        // Используем существующий triggerNoVeto для этой цели
        triggerNoVeto += v; 
    }
    
    void AddTriggerAndVeto(G4int v) { 
        // Увеличиваем счетчик событий "кристалл и вето"
        triggerWithVeto += v; 
    }
    
    void AddTriggeredTriggerOnly(double E_MeV) { 
        // Заполняем гистограмму для событий "только кристалл"
        AddTriggerNoVeto(E_MeV); 
    }

    // НОВЫЕ МЕТОДЫ ДЛЯ ОПТИКИ
    void AddTriggerOnlyOpt(G4int v) { 
        triggerNoVetoOpt += v; 
    }
    
    void AddTriggerAndVetoOpt(G4int v) { 
        triggerWithVetoOpt += v; 
    }
    
    void AddTriggeredTriggerOnlyOpt(double E_MeV) { 
        AddTriggerNoVetoOpt(E_MeV); 
    }

    [[nodiscard]] const ParticleCounts& GetCounts() const { return totals; }
    [[nodiscard]] const ParticleCounts& GetOptCounts() const { return totalsOpt; }

    [[nodiscard]] const std::vector<double>& GetEffArea() const { return effArea; }
    [[nodiscard]] const std::vector<double>& GetEffAreaOpt() const { return effAreaOpt; }

private:
    G4Accumulable<G4double> triggerNoVeto{0};
    G4Accumulable<G4double> triggerWithVeto{0};
    G4Accumulable<G4double> triggerNoVetoOpt{0};
    G4Accumulable<G4double> triggerWithVetoOpt{0};

    ParticleCounts totals{};
    ParticleCounts totalsOpt{};

    double EminMeV{0.0};
    double EmaxMeV{0.0};
    double area{0.0};

    double logEmin{0.0};
    double logEmax{0.0};
    double invDlogE{0.0};

    std::vector<G4Accumulable<G4double>> genCounts;
    std::vector<G4Accumulable<G4double>> trigCounts;
    std::vector<G4Accumulable<G4double>> trigOptCounts;
    std::vector<G4double> effArea;
    std::vector<G4double> effAreaOpt;

    [[nodiscard]] int FindBinLog(double E_MeV) const;
    [[nodiscard]] double BinCenterMeV(int i) const;
    [[nodiscard]] double BinWidthMeV(int i) const;

    void BookAccumulables();
    void FillDerivedHists();
};

#endif