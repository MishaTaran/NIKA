#include "Flux/Flux.hh"

ParticleInfo Flux::GenerateParticle() {
    auto *pt = G4ParticleTable::GetParticleTable();

    ParticleInfo info;
    info.energy = SampleEnergy();
    info.def = pt->FindParticle(particle);
    info.name = particle;
    info.pdg = info.def->GetPDGEncoding();
    return info;
}

G4String Flux::Trim(const G4String &_s) {
    const size_t start = _s.find_first_not_of(" \t\r\n");
    if (start == G4String::npos) return "";
    const size_t end = _s.find_last_not_of(" \t\r\n");
    return _s.substr(start, end - start + 1);
}

void Flux::LoadFileIfNeeded(const G4String &filepath) {
    std::ifstream fin(filepath);
    if (!fin.is_open()) {
        G4Exception("FluxBase::loadFileIfNeeded", "FILE_OPEN_FAIL",
                    JustWarning, ("Cannot open " + filepath).c_str());
        return;
    }

    G4String line;
    while (std::getline(fin, line)) {
        size_t pos = line.find(':');
        if (pos == G4String::npos)
            pos = line.find('=');
        if (pos == G4String::npos) continue;

        G4String key = Trim(line.substr(0, pos));
        G4String val = Trim(line.substr(pos + 1));
        if (!key.empty() && !val.empty())
            cache[key] = val;
    }
    fin.close();
}

G4double Flux::GetParam(const G4String &filepath,
                        const G4String &key,
                        G4double defaultValue) {
    LoadFileIfNeeded(filepath);
    try {
        if (cache.count(key)) {
            return std::stod(cache[key]);
        }
    } catch (...) {
    }
    return defaultValue;
}

G4String Flux::GetParam(const G4String &filepath,
                        const G4String &key,
                        const G4String &defaultValue) {
    LoadFileIfNeeded(filepath);
    try {
        if (cache.count(key)) {
            return cache[key];
        }
    } catch (...) {
    }
    return defaultValue;
}
