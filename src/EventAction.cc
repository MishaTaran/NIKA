#include "EventAction.hh"

using namespace Sizes;
using namespace Configuration;

EventAction::EventAction(AnalysisManager* an, RunAction* r) : analysisManager(an), run(r) {
    detMap = {
        {"TriggerSD/EdepHits", 0, "Trigger"},
        {"SideVetoSD/EdepHits", 1, "SideVeto"},
        {"UpperVetoSD/EdepHits", 2, "UpperVeto"},
        {"BottomVetoSD/EdepHits", 3, "BottomVeto"},
    };
    HCIDs.assign(detMap.size(), -1);
}

void EventAction::BeginOfEventAction(const G4Event*) {
    nPrimaries = 0;
    nInteractions = 0;
    nEdepHits = 0;

    // Инициализация флагов для энергетических детекторов
    hasTrigger = false;
    hasSideVeto = false;
    hasBottomVeto = false;
    hasUpperVeto = false;

    // Инициализация флагов для оптических детекторов
    hasTriggerOpt = false;
    hasSideVetoOpt = false;
    hasBottomVetoOpt = false;
    hasUpperVetoOpt = false;

    // Инициализация слоев триггера
    for (int i = 0; i < 3; ++i) {
        hasTriggerLayer[i] = false;
        stripE[i] = 0.0;
    }
    detE.clear();
}

void EventAction::EndOfEventAction(const G4Event* evt) {
    const int eventID = evt->GetEventID();

    WritePrimaries_(eventID);
    nPrimaries = static_cast<int>(primBuf.size());

    double primaryE_MeV = -1.0;
    if (!primBuf.empty()) {
        primaryE_MeV = primBuf.front().E_MeV;
        if (run) {
            run->AddGenerated(primaryE_MeV);
        }
    }
    primBuf.clear();

    nInteractions = WriteInteractions_(eventID);
    interBuf.clear();

    if (savePhotons) {
        nPhotons = WritePhotons_(eventID);
        photonBuf.clear();
        WritePhotonsCount_(eventID);
        photonCountBuf = {0, 0, 0};
    }

    nEdepHits = WriteEdepFromSD_(evt, eventID);

    if (saveSecondaries) {
        analysisManager->FillEventRow(eventID, nPrimaries, nInteractions, nEdepHits);
    }

    // Счетчики для энергетических детекторов
    if (run && hasTrigger && !hasSideVeto && !hasBottomVeto && !hasUpperVeto) run->AddTriggerOnly(1);
    if (run && hasTrigger && hasSideVeto && hasBottomVeto && hasUpperVeto) run->AddTriggerAndVeto(1);

    if (primaryE_MeV > 0.0) {
        if (run && hasTrigger && !hasSideVeto && !hasBottomVeto && !hasUpperVeto) {
            run->AddTriggeredTriggerOnly(primaryE_MeV);
        }
    }

    // Счетчики для оптических детекторов
    if (useOptics) {
        WriteSiPMFromSD_(eventID);

        if (run && hasTriggerOpt && !hasSideVetoOpt && !hasBottomVetoOpt && !hasUpperVetoOpt) run->AddTriggerOnlyOpt(1);
        if (run && hasTriggerOpt && hasSideVetoOpt && hasBottomVetoOpt && hasUpperVetoOpt) run->AddTriggerAndVetoOpt(1);

        if (primaryE_MeV > 0.0) {
            if (run && hasTriggerOpt && !hasSideVetoOpt && !hasBottomVetoOpt && !hasUpperVetoOpt)
                run->AddTriggeredTriggerOnlyOpt(primaryE_MeV);
        }

        // === NEW: вывод счётчиков SiPM в терминал ===
        //PrintSiPMCounts_(eventID, primaryE_MeV);
    }
}

void EventAction::WritePrimaries_(int eventID) {
    for (const auto& p : primBuf) {
        analysisManager->FillPrimaryRow(eventID, p.name, p.E_MeV, p.dir, p.pos_mm);
    }
}

int EventAction::WriteInteractions_(int eventID) {
    if (saveSecondaries) {
        for (const auto& r : interBuf) {
            analysisManager->FillInteractionRow(eventID,
                                                r.trackID, r.parentID,
                                                r.process, r.volumeName, r.pos_mm,
                                                r.secIndex, r.secName,
                                                r.secE_MeV, r.secDir);
        }
    }
    return static_cast<int>(interBuf.size());
}

int EventAction::WritePhotonsCount_(int eventID) {
    if (savePhotons) {
        analysisManager->FillPhotonCountRow(eventID, photonCountBuf[0], photonCountBuf[1], photonCountBuf[2]);
    }
    return static_cast<int>(photonCountBuf.size());
}

int EventAction::WritePhotons_(int eventID) {
    for (const auto& photon : photonBuf) {
        analysisManager->FillPhotonRow(eventID, photon.photonID, photon.detName, photon.detCh, photon.energy,
                                       photon.pos_mm.x(), photon.pos_mm.y(), photon.pos_mm.z());
    }
    return static_cast<int>(photonBuf.size());
}

int EventAction::WriteEdepFromSD_(const G4Event* evt, int eventID) {
    auto* hce = evt->GetHCofThisEvent();
    if (!hce) return 0;

    auto* sdm = G4SDManager::GetSDMpointer();

    for (size_t i = 0; i < detMap.size(); ++i) {
        if (HCIDs[i] < 0) {
            const auto& hcName = std::get<0>(detMap[i]);
            HCIDs[i] = sdm->GetCollectionID(hcName);
        }
    }

    int nHitsTotal = 0;

    for (size_t i = 0; i < detMap.size(); ++i) {
        const int hcID = HCIDs[i];
        if (hcID < 0) continue;

        auto* hc = dynamic_cast<SDHitCollection*>(hce->GetHC(hcID));
        if (!hc) continue;

        const auto& det_name = std::get<2>(detMap[i]);

        const auto N = hc->GetSize();
        for (unsigned j = 0; j < N; ++j) {
            auto* h = (*hc)[j];
            double edep_MeV = h->edep / MeV;

            if ((det_name == "Veto" or det_name == "BottomVeto") and edep_MeV <= eVetoThreshold) {
                edep_MeV = 0;
            }
            if (det_name == "Trigger" and edep_MeV <= eTriggerThreshold) {
                edep_MeV = 0;
            }
            if (edep_MeV > 0.0) {
                if (det_name == "Trigger") {
                    MarkTrigger();
                    // === пишем хит в trigger_hits с номером слоя ===
                    const int layer = (h->layer >= 0) ? h->layer : h->volumeID;
                    analysisManager->FillTriggerHitRow(eventID, layer, edep_MeV);
                }
                else if (det_name == "SideVeto" or det_name == "BottomVeto" or det_name == "UpperVeto") {
                    MarkVeto();
                }
                analysisManager->FillEdepRow(eventID, det_name, edep_MeV);
            }
        }
        nHitsTotal += static_cast<int>(N);
    }

    return nHitsTotal;
}

void EventAction::WriteSiPMFromSD_(int eventID) {
    auto* sdm = G4SDManager::GetSDMpointer();
    if (!sdm) return;

    auto* sdBase = sdm->FindSensitiveDetector("SiPMOpticalSD", false);
    auto* sipmSD = dynamic_cast<SiPMOpticalSD*>(sdBase);
    if (!sipmSD) return;

    int npeC = sipmSD->GetNpeTrigger();
    int npeS = sipmSD->GetNpeSideVeto();
    int npeU = sipmSD->GetNpeUpperVeto();
    int npeB = sipmSD->GetNpeBottomVeto();

    npeC = npeC > oTriggerThreshold ? npeC : 0;
    npeS = npeS > oSideVetoThreshold ? npeS : 0;
    npeU = npeU > oUpperVetoThreshold ? npeU : 0;
    npeB = npeB > oBottomVetoThreshold ? npeB : 0;

    if (npeC > 0) MarkTriggerOpt();

    // Любое veto (SideVeto, UpperVeto, BottomVeto) считается veto
    if (npeS > 0 || npeU > 0 || npeB > 0) MarkVetoOpt();

    int npeVetoTotal = npeS + npeU;  // SideVeto + UpperVeto
    analysisManager->FillSiPMEventRow(eventID, npeC, npeVetoTotal, npeB);

    // ---- Trigger каналы ----
    for (const auto& kv : sipmSD->GetPerChannelTrigger()) {
        analysisManager->FillSiPMChannelRow(eventID, "Trigger", kv.first, kv.second);
    }

    // ---- SideVeto каналы ----
    for (const auto& kv : sipmSD->GetPerChannelSideVeto()) {
        analysisManager->FillSiPMChannelRow(eventID, "SideVeto", kv.first, kv.second);
    }

    // ---- UpperVeto каналы ----
    for (const auto& kv : sipmSD->GetPerChannelUpperVeto()) {
        analysisManager->FillSiPMChannelRow(eventID, "UpperVeto", kv.first, kv.second);
    }

    // ---- BottomVeto каналы ----
    for (const auto& kv : sipmSD->GetPerChannelBottomVeto()) {
        analysisManager->FillSiPMChannelRow(eventID, "BottomVeto", kv.first, kv.second);
    }
}

// ============================================================
// NEW: вывод счётчиков SiPM в терминал
// ============================================================
void EventAction::PrintSiPMCounts_(int eventID, double primaryE_MeV) {
    auto* sdm = G4SDManager::GetSDMpointer();
    if (!sdm) return;

    auto* sdBase = sdm->FindSensitiveDetector("SiPMOpticalSD", false);
    auto* sipmSD = dynamic_cast<SiPMOpticalSD*>(sdBase);
    if (!sipmSD) return;

    const int npeTrigger    = sipmSD->GetNpeTrigger();
    const int npeSideVeto   = sipmSD->GetNpeSideVeto();
    const int npeUpperVeto  = sipmSD->GetNpeUpperVeto();
    const int npeBottomVeto = sipmSD->GetNpeBottomVeto();

    const auto& layers = sipmSD->GetNpeTriggerLayers();
    const int npeL0 = layers[0];
    const int npeL1 = layers[1];
    const int npeL2 = layers[2];

    G4cout << "========================================\n"
           << " Event " << eventID
           << "  E_primary = " << primaryE_MeV << " MeV\n"
           << "----------------------------------------\n"
           << "  Trigger     : " << npeTrigger
           << "   (L0=" << npeL0
           << ", L1=" << npeL1
           << ", L2=" << npeL2 << ")\n"
           << "  SideVeto    : " << npeSideVeto << "\n"
           << "  UpperVeto   : " << npeUpperVeto << "\n"
           << "  BottomVeto  : " << npeBottomVeto << "\n"
           << "========================================"
           << G4endl;
}