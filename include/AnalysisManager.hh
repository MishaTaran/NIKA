#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include <G4AnalysisManager.hh>
#include <G4AccumulableManager.hh>
#include <G4Accumulable.hh>
#include <G4UnitsTable.hh>
#include <CLHEP/Units/SystemOfUnits.h>
#include <globals.hh>
#include <Sizes.hh>
#include <Configuration.hh>

class AnalysisManager {
public:
    explicit AnalysisManager(const std::string& fileName);
    AnalysisManager(const std::string&, int, double, double);
    ~AnalysisManager() = default;

    void Open();
    void Close();

    void FillEventRow(G4int eventID, G4int nPrimaries, G4int nInteractions, G4int nEdepHits);

    void FillPrimaryRow(G4int eventID, const G4String& primaryName,
                        G4double E_MeV, const G4ThreeVector& dir,
                        const G4ThreeVector& pos_mm);

    void FillEdepRow(G4int eventID, const G4String& detName, G4double edep_MeV);

    void FillTriggerHitRow(G4int eventID, G4int layer, G4double edep_MeV);

    void FillGenEnergyHist(G4double E_MeV, G4double weight = 1.0);
    void FillTrigEnergyHist(G4double E_MeV, G4double weight = 1.0);
    void FillTrigOptEnergyHist(G4double E_MeV, G4double weight = 1.0);

    void FillEffAreaHist(G4double E_MeV, G4double value);
    void FillEffAreaOptHist(G4double E_MeV, G4double value);

    void FillSensitivityHist(G4double E_MeV, G4double value);
    void FillSensitivityOptHist(G4double E_MeV, G4double value);

    void FillSiPMEventRow(int eventID, int npeC, int npeV, int npeB);
    void FillSiPMChannelRow(int eventID, const G4String& subdet, int ch, int npe);

    void FillInteractionRow(G4int eventID,
                            G4int trackID, G4int parentID,
                            const G4String& process,
                            const G4String& volumeName,
                            const G4ThreeVector& x_mm,
                            G4int secIndex, const G4String& secName,
                            G4double secE_MeV, const G4ThreeVector& secDir);

    void FillPhotonCountRow(G4int eventID,
                            G4int npeTrigger, G4int npeVeto,
                            G4int npeBottomVeto);

    void FillPhotonRow(G4int eventID, G4int photonID, const G4String& det_name, G4int det_ch,
                       G4double energy_eV, G4double x_mm, G4double y_mm, G4double z_mm);

private:
    G4String fileName;

    G4int eventNT{-1};
    G4int primaryNT{-1};
    G4int interactionsNT{-1};
    G4int photonsCountNT{-1};
    G4int photonsNT{-1};
    G4int edepNT{-1};
    G4int triggerHitsNT{-1};

    G4int SiPMEventNT{-1};
    G4int SiPMChannelNT{-1};

    G4int genEnergyHist{-1};
    G4int trigEnergyHist{-1};
    G4int trigOptEnergyHist{-1};

    G4int effAreaHist{-1};
    G4int effAreaOptHist{-1};

    G4int sensitivityHist{-1};
    G4int sensitivityOptHist{-1};

    G4int nBins{1000};
    G4double xMin{0};
    G4double xMax{1000 * MeV};

    void Book();
};

#endif // ANALYSISMANAGER_HH