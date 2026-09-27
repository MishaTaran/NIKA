#include "SiPMOpticalSD.hh"
#include "EventAction.hh"
#include <G4EventManager.hh>

SiPMOpticalSD::SiPMOpticalSD(const G4String& name)
    : G4VSensitiveDetector(name) {}

void SiPMOpticalSD::Initialize(G4HCofThisEvent*) {
    npeTrigger = npeSideVeto = npeUpperVeto = npeBottomVeto = 0;
    npeTriggerLayer = {0, 0, 0};   // === NEW ===
    perChTrigger.clear();
    perChSideVeto.clear();
    perChUpperVeto.clear();
    perChBottomVeto.clear();
}

G4OpBoundaryProcess* SiPMOpticalSD::GetBoundaryProcess() {
    if (boundary) return boundary;
    auto* pm = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
    if (!pm) return nullptr;

    auto* plist = pm->GetProcessList();
    if (!plist) return nullptr;

    for (int i = 0; i < plist->size(); ++i) {
        auto* p = (*plist)[i];
        if (p && p->GetProcessName() == "OpBoundary") {
            boundary = static_cast<G4OpBoundaryProcess*>(p);
            break;
        }
    }
    return boundary;
}

SiPMGroup SiPMOpticalSD::ClassifyByPVName(const G4VPhysicalVolume* pv) {
    if (!pv) return SiPMGroup::Unknown;
    const std::string name = pv->GetName();
    if (name.find("TriggerSiPM") != std::string::npos) return SiPMGroup::Trigger;
    if (name.find("SideVetoSiPM") != std::string::npos) return SiPMGroup::SideVeto;
    if (name.find("UpperVetoSiPM") != std::string::npos) return SiPMGroup::UpperVeto;
    if (name.find("BottomVetoSiPM") != std::string::npos) return SiPMGroup::BottomVeto;

    return SiPMGroup::Unknown;
}

// ============================================================
// ФУНКЦИЯ ДЛЯ ПОЛУЧЕНИЯ НОМЕРА КОПИИ SiPM ЧЕРЕЗ TOUCHABLE
// ============================================================
int SiPMOpticalSD::GetSiPMCopyNumber(const G4StepPoint* point) {
    if (!point) return -1;

    auto touchable = point->GetTouchableHandle();
    if (!touchable) return -1;

    // Проходим по всей иерархии от корня к текущему объему
    for (int i = touchable->GetHistoryDepth() - 1; i >= 0; --i) {
        auto* pv = touchable->GetVolume(i);
        if (!pv) continue;

        G4String name = pv->GetName();
        // Ищем наш SiPM по имени (формат: "TriggerSiPM_PVP_N", "SideVetoSiPM_PVP_N" и т.д.)
        if (name.find("SiPM_PVP_") != std::string::npos) {
            return pv->GetCopyNo();
        }
    }

    return -1;
}

G4bool SiPMOpticalSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
    if (!step) return false;

    auto* track = step->GetTrack();
    if (!track) return false;

    auto* def = track->GetDefinition();
    if (!def) return false;

    if (def != G4OpticalPhoton::Definition()) {
        return false;
    }

    auto* post = step->GetPostStepPoint();
    auto* pre = step->GetPreStepPoint();
    if (!pre || !post) return false;

    if (post->GetStepStatus() != fGeomBoundary) {
        return false;
    }

    auto* b = GetBoundaryProcess();
    if (!b) return false;

    G4OpBoundaryProcessStatus status = b->GetStatus();
    if (status != Detection && status != Absorption) {
        return false;
    }
    auto* prePV = pre->GetPhysicalVolume();
    auto* postPV = post->GetPhysicalVolume();

    // ============================================================
    // ПОЛУЧАЕМ НОМЕР КОПИИ ЧЕРЕЗ TOUCHABLE
    // ============================================================
    int ch = GetSiPMCopyNumber(pre);
    if (ch < 0) ch = GetSiPMCopyNumber(post);
    // Если не получилось через Touchable, пробуем прямой доступ к copy number
    if (ch < 0) {
        if (prePV) ch = prePV->GetCopyNo();
        if (ch < 0 && postPV) ch = postPV->GetCopyNo();
    }

    // ============================================================
    // ОПРЕДЕЛЯЕМ ГРУППУ SiPM
    // ============================================================
    SiPMGroup grp = SiPMGroup::Unknown;

    // Способ 1: По имени физического объема
    if (prePV) {
        grp = ClassifyByPVName(prePV);
    }
    if (grp == SiPMGroup::Unknown && postPV) {
        grp = ClassifyByPVName(postPV);
    }

    // Способ 2: По логическому объему
    if (grp == SiPMGroup::Unknown && SiPMLV) {
        auto* preLV = prePV ? prePV->GetLogicalVolume() : nullptr;
        auto* postLV = postPV ? postPV->GetLogicalVolume() : nullptr;
        if (preLV == SiPMLV || postLV == SiPMLV) {
            grp = ClassifyByPVName(prePV);
            if (grp == SiPMGroup::Unknown) grp = ClassifyByPVName(postPV);
        }
    }

    // Способ 3: По номеру копии (если известны диапазоны)
    if (grp == SiPMGroup::Unknown && ch >= 0) {
        // 24 Trigger (0-23), 8 SideVeto (24-31),
        // 4 UpperVeto (32-35), 4 BottomVeto (36-39)
        if (ch < 24) grp = SiPMGroup::Trigger;
        else if (ch < 32) grp = SiPMGroup::SideVeto;
        else if (ch < 36) grp = SiPMGroup::UpperVeto;
        else if (ch < 40) grp = SiPMGroup::BottomVeto;
    }

    // ============================================================
    // СОХРАНЯЕМ РЕЗУЛЬТАТЫ
    // ============================================================
    G4String detName = "";
    if (grp == SiPMGroup::Trigger) {
        detName = "Trigger";
        ++npeTrigger;

        // === NEW: определяем слой триггера по номеру канала ===
        // Порядок размещения в Detector::PlaceAllSiPMs():
        //   layer 0 -> copyNo 0..7
        //   layer 1 -> copyNo 8..15
        //   layer 2 -> copyNo 16..23
        int layer = -1;
        if (ch >= 0) {
            if (ch < 8)       layer = 0;
            else if (ch < 16) layer = 1;
            else if (ch < 24) layer = 2;
        }
        if (layer >= 0 && layer < 3) {
            ++npeTriggerLayer[layer];
        }

        if (ch >= 0) ++perChTrigger[ch];
    } else if (grp == SiPMGroup::SideVeto) {
        detName = "SideVeto";
        ++npeSideVeto;
        if (ch >= 0) ++perChSideVeto[ch];
    } else if (grp == SiPMGroup::UpperVeto) {
        detName = "UpperVeto";
        ++npeUpperVeto;
        if (ch >= 0) ++perChUpperVeto[ch];
    } else if (grp == SiPMGroup::BottomVeto) {
        detName = "BottomVeto";
        ++npeBottomVeto;
        if (ch >= 0) ++perChBottomVeto[ch];
    } else {
        return false;
    }

    if (Configuration::savePhotons) {
        if (auto* ea = dynamic_cast<EventAction*>(
            G4EventManager::GetEventManager()->GetUserEventAction())) {
            PhotonRec rec;
            rec.photonID = track->GetTrackID();
            rec.detName = detName;
            rec.detCh = ch;
            rec.energy = track->GetTotalEnergy() / eV;
            rec.pos_mm = post->GetPosition();
            ea->photonBuf.emplace_back(std::move(rec));
        }
    }
    track->SetTrackStatus(fStopAndKill);
    return true;
}