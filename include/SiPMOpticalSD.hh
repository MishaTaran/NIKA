#ifndef SiPMOpticalSD_h
#define SiPMOpticalSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include <map>
#include <string>

// Перечисление для типов SiPM
enum class SiPMGroup {
    Unknown,
    Trigger,
    SideVeto,
    UpperVeto,
    BottomVeto
};

class SiPMOpticalSD : public G4VSensitiveDetector {
public:
    SiPMOpticalSD(const G4String& name);
    virtual ~SiPMOpticalSD() = default;

    // Методы Geant4
    virtual void Initialize(G4HCofThisEvent*) override;
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

    // Геттеры для общего количества фотоэлектронов (NPE)
    G4int GetNpeTrigger() const { return npeTrigger; }
    G4int GetNpeSideVeto() const { return npeSideVeto; }
    G4int GetNpeUpperVeto() const { return npeUpperVeto; }
    G4int GetNpeBottomVeto() const { return npeBottomVeto; }
    
    // Геттеры для per-channel счетчиков (оригинальные названия)
    const std::map<G4int, G4int>& GetPerChannelTrigger() const { return perChTrigger; }
    const std::map<G4int, G4int>& GetPerChannelSideVeto() const { return perChSideVeto; }
    const std::map<G4int, G4int>& GetPerChannelUpperVeto() const { return perChUpperVeto; }
    const std::map<G4int, G4int>& GetPerChannelBottomVeto() const { return perChBottomVeto; }

    // Сеттер для логического объема SiPM
    void SetSiPMLV(G4LogicalVolume* lv) { SiPMLV = lv; }

private:
    // Получение граничного процесса для оптических фотонов
    G4OpBoundaryProcess* GetBoundaryProcess();
    
    // Классификация SiPM по имени физического объема
    SiPMGroup ClassifyByPVName(const G4VPhysicalVolume* pv);
    
    // Получение номера канала (copy number) SiPM через Touchable
    G4int GetSiPMCopyNumber(const G4StepPoint* point);

    // Указатель на граничный процесс
    G4OpBoundaryProcess* boundary = nullptr;
    
    // Логический объем SiPM
    G4LogicalVolume* SiPMLV = nullptr;

    // Счетчики общего количества фотоэлектронов
    G4int npeTrigger = 0;
    G4int npeSideVeto = 0;
    G4int npeUpperVeto = 0;
    G4int npeBottomVeto = 0;

    // Per-channel счетчики
    std::map<G4int, G4int> perChTrigger;
    std::map<G4int, G4int> perChSideVeto;
    std::map<G4int, G4int> perChUpperVeto;
    std::map<G4int, G4int> perChBottomVeto;
};

#endif