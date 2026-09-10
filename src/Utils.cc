#include "Utils.hh"

std::string Utils::Trim(std::string s_) {
    auto not_space = [](unsigned char c) {
        return !std::isspace(c);
    };
    s_.erase(s_.begin(), std::find_if(s_.begin(), s_.end(), not_space));
    s_.erase(std::find_if(s_.rbegin(), s_.rend(), not_space).base(), s_.end());
    return s_;
}

bool Utils::StartsWith(const std::string& s_, const char* prefix) {
    const size_t n = std::strlen(prefix);
    return s_.size() >= n && s_.compare(0, n, prefix) == 0;
}

Utils::ConstMap Utils::ReadConstFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Can't open const file: " + filename);

    ConstMap out;
    std::string line;
    size_t lineno = 0;

    while (std::getline(file, line)) {
        ++lineno;

        // Strip comments starting with '#'
        const auto hashPos = line.find('#');
        if (hashPos != std::string::npos) line = line.substr(0, hashPos);

        line = Trim(line);
        if (line.empty()) continue;

        // Allow KEY=VALUE form by replacing '=' with space
        for (char& ch : line) if (ch == '=') ch = ' ';

        std::istringstream iss(line);
        std::string key, valTok, unitTok;

        if (!(iss >> key >> valTok)) {
            continue;
        }
        if (iss >> unitTok) {
            // ok
        } else {
            unitTok.clear(); // dimensionless
        }

        key = Trim(key);
        valTok = Trim(valTok);
        unitTok = Trim(unitTok);

        double v = 0.0;
        try {
            v = std::stod(valTok);
        }
        catch (...) {
            throw std::runtime_error("Bad numeric value in " + filename + ":" + std::to_string(lineno));
        }

        const G4double factor = UnitFactor(unitTok);
        out[key] = v * factor;
    }

    file.close();
    return out;
}

G4double Utils::GetRequired(const ConstMap& m, const std::string& key) {
    auto it = m.find(key);
    if (it == m.end()) throw std::runtime_error("Missing required optical const: " + key);
    return it->second;
}

G4double Utils::GetOr(const ConstMap& m, const std::string& key, G4double def) {
    auto it = m.find(key);
    return it == m.end() ? def : it->second;
}

G4double Utils::UnitFactor(const std::string& unitToken) {
    if (unitToken.empty() || unitToken == "-") return 1.0;

    // length
    if (unitToken == "mm") return mm;
    if (unitToken == "cm") return cm;
    if (unitToken == "m") return m;
    if (unitToken == "um") return um;

    // time
    if (unitToken == "ps") return ps;
    if (unitToken == "ns") return ns;
    if (unitToken == "us") return us;
    if (unitToken == "ms") return ms;
    if (unitToken == "s") return s;

    // energy
    if (unitToken == "eV") return eV;
    if (unitToken == "keV") return keV;
    if (unitToken == "MeV") return MeV;

    // inverse energy (for yields)
    if (unitToken == "1/eV") return 1.0 / eV;
    if (unitToken == "1/keV") return 1.0 / keV;
    if (unitToken == "1/MeV") return 1.0 / MeV;

    // inverse energy (for Briks constant)
    if (unitToken == "mm/eV") return mm / eV;
    if (unitToken == "mm/keV") return mm / keV;
    if (unitToken == "mm/MeV") return mm / MeV;

    // common typo support
    if (unitToken == "1/MEV") return 1.0 / MeV;

    throw std::runtime_error("Unknown unit token: '" + unitToken + "'");
}

Utils::Table Utils::ReadCSV(const std::string& filename,
                            G4double valueScale,
                            bool clampNonNegative) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Can't open CSV file: " + filename);

    struct Row {
        G4double E;
        G4double V;
    };
    std::vector<Row> rows;
    rows.reserve(512);

    std::string line;
    while (std::getline(file, line)) {
        line = Trim(line);
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        if (line.size() >= 2 && line[0] == '/' && line[1] == '/') continue;

        std::stringstream ss(line);
        std::string col1, col2;
        if (!std::getline(ss, col1, ',')) continue;
        if (!std::getline(ss, col2)) continue;

        col1 = Trim(col1);
        col2 = Trim(col2);
        if (col1.empty() || col2.empty()) continue;

        const double e_eV = std::stod(col1);
        double v = std::stod(col2);
        if (clampNonNegative && v < 0.0) v = 0.0;

        rows.push_back(Row{
                           static_cast<G4double>(e_eV) * eV,
                           v * valueScale
                       });
    }
    file.close();

    if (rows.size() < 2) throw std::runtime_error("CSV has <2 valid rows: " + filename);

    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
        return a.E < b.E;
    });

    // de-dup energies
    const G4double eps = 1e-12 * eV;
    std::vector<Row> uniq;
    uniq.reserve(rows.size());
    for (const auto& r : rows) {
        if (uniq.empty() || std::abs(r.E - uniq.back().E) > eps) uniq.push_back(r);
    }
    if (uniq.size() < 2) throw std::runtime_error("CSV energies not increasing after de-dup: " + filename);

    Table t;
    t.E.reserve(uniq.size());
    t.V.reserve(uniq.size());
    for (const auto& r : uniq) {
        t.E.push_back(r.E);
        t.V.push_back(r.V);
    }

    ValidateSameSize(t, filename);
    return t;
}

Utils::EmissionTables Utils::ReadEmissionCSV(const std::string& filename,
                                             G4double valueScale,
                                             bool clampNonNegative) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Can't open emission CSV file: " + filename);

    struct Row {
        G4double E;
        G4double v1;
        G4double v2;
    };

    std::vector<Row> rows;
    rows.reserve(512);

    std::string line;
    while (std::getline(file, line)) {
        line = Trim(line);
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        if (line.size() >= 2 && line[0] == '/' && line[1] == '/') continue;

        std::stringstream ss(line);
        std::string colE, col1, col2;

        if (!std::getline(ss, colE, ',')) continue;
        if (!std::getline(ss, col1, ',')) continue;
        std::getline(ss, col2);

        colE = Trim(colE);
        col1 = Trim(col1);
        col2 = Trim(col2);

        if (colE.empty() || col1.empty()) continue;

        double e_eV = 0.0;
        double v1 = 0.0;

        try {
            e_eV = std::stod(colE);
            v1 = std::stod(col1);
        }
        catch (...) {
            continue;
        }

        double v2 = v1;
        if (!col2.empty()) {
            try {
                v2 = std::stod(col2);
            }
            catch (...) {
                v2 = v1;
            }
        }

        if (clampNonNegative) {
            if (v1 < 0.0) v1 = 0.0;
            if (v2 < 0.0) v2 = 0.0;
        }

        rows.push_back(Row{
                           e_eV * eV,
                           v1 * valueScale,
                           v2 * valueScale
                       });
    }

    file.close();

    if (rows.size() < 2) throw std::runtime_error("Emission CSV has <2 valid rows: " + filename);

    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
        return a.E < b.E;
    });

    const G4double eps = 1e-12 * eV;
    std::vector<Row> uniq;
    uniq.reserve(rows.size());

    for (const auto& r : rows) {
        if (uniq.empty() || std::abs(r.E - uniq.back().E) > eps) {
            uniq.push_back(r);
        }
    }

    if (uniq.size() < 2) throw std::runtime_error("Emission CSV energies not increasing after de-dup: " + filename);

    EmissionTables out;
    out.c1.E.reserve(uniq.size());
    out.c1.V.reserve(uniq.size());
    out.c2.E.reserve(uniq.size());
    out.c2.V.reserve(uniq.size());

    for (const auto& r : uniq) {
        out.c1.E.push_back(r.E);
        out.c1.V.push_back(r.v1);
        out.c2.E.push_back(r.E);
        out.c2.V.push_back(r.v2);
    }

    ValidateSameSize(out.c1, filename + " [component1]");
    ValidateSameSize(out.c2, filename + " [component2]");

    return out;
}

void Utils::ValidateSameSize(const Table& t, const std::string& filename) {
    if (t.E.size() != t.V.size() || t.E.size() < 2) {
        throw std::runtime_error("Bad table sizes in: " + filename);
    }
    for (size_t i = 1; i < t.E.size(); ++i) {
        if (!(t.E[i] > t.E[i - 1])) {
            throw std::runtime_error("Energies are not strictly increasing in: " + filename);
        }
    }
}

void Utils::NormalizeMaxToOne(Table& t) {
    G4double mx = 0.0;
    for (auto v : t.V) if (v > mx) mx = v;
    if (mx <= 0.0) return;
    for (auto& v : t.V) v /= mx;
}

void Utils::ClampNonNegative(Table& t) {
    for (auto& v : t.V) if (v < 0.0) v = 0.0;
}

Utils::Table Utils::MakeConstantTable(G4double Emin_eV, G4double Emax_eV, G4double value) {
    Table t;
    t.E = {Emin_eV * eV, Emax_eV * eV};
    t.V = {value, value};
    return t;
}

void Utils::ApplyMaterialTable(G4Material* mat,
                               const Table& rindex,
                               const Table* abslength) {
    if (!mat) return;
    auto* mpt = mat->GetMaterialPropertiesTable();
    if (!mpt) mpt = new G4MaterialPropertiesTable();

    mpt->AddProperty("RINDEX", rindex.E, rindex.V, rindex.E.size());
    if (abslength) {
        mpt->AddProperty("ABSLENGTH", abslength->E, abslength->V, abslength->E.size());
    }
    mat->SetMaterialPropertiesTable(mpt);
}

void Utils::AddConstIfPresent(G4MaterialPropertiesTable* mpt,
                              const ConstMap& c,
                              const std::string& key) {
    auto it = c.find(key);
    if (it != c.end()) {
        mpt->AddConstProperty(key.c_str(), it->second);
    }
}

void Utils::ApplyScintillation(G4Material* mat,
                               G4MaterialPropertiesTable* mpt,
                               const ConstMap& c,
                               const Table& scintComponent1,
                               const Table& scintComponent2,
                               bool requireYield) {
    if (!mat || !mpt) return;

    mpt->AddProperty("SCINTILLATIONCOMPONENT1",
                     scintComponent1.E, scintComponent1.V, scintComponent1.E.size());

    const bool hasYield2Key = (c.find("SCINTILLATIONYIELD2") != c.end());
    const bool hasTau2Key = (c.find("SCINTILLATIONTIMECONSTANT2") != c.end());

    const G4double y1 = GetOr(c, "SCINTILLATIONYIELD1", 1.0);
    const G4double y2 = GetOr(c, "SCINTILLATIONYIELD2", 0.0);

    G4bool has2 = false;
    if ((hasYield2Key || hasTau2Key) && (y2 > 0.0 || y1 < 1.0)) {
        has2 = true;
    }

    if (has2) {
        mpt->AddProperty("SCINTILLATIONCOMPONENT2",
                         scintComponent2.E, scintComponent2.V, scintComponent2.E.size());
    }

    if (requireYield) {
        const auto it = c.find("SCINTILLATIONYIELD");
        if (it == c.end()) {
            G4cerr << "[Utils::ApplyScintillation] ⚠ WARNING: SCINTILLATIONYIELD not found!" << G4endl;
            mpt->AddConstProperty("SCINTILLATIONYIELD", 10000.0);
        } else {
            mpt->AddConstProperty("SCINTILLATIONYIELD", it->second);
            G4cout << "[Utils::ApplyScintillation] SCINTILLATIONYIELD = " << it->second << G4endl;
        }
    }

    AddConstIfPresent(mpt, c, "RESOLUTIONSCALE");
    AddConstIfPresent(mpt, c, "SCINTILLATIONTIMECONSTANT1");
    AddConstIfPresent(mpt, c, "SCINTILLATIONYIELD1");

    if (has2) {
        AddConstIfPresent(mpt, c, "SCINTILLATIONTIMECONSTANT2");

        if (hasYield2Key) {
            AddConstIfPresent(mpt, c, "SCINTILLATIONYIELD2");
        } else {
            G4double y2_derived = 1.0 - y1;
            if (y2_derived < 0.0) y2_derived = 0.0;
            if (y2_derived > 1.0) y2_derived = 1.0;
            mpt->AddConstProperty("SCINTILLATIONYIELD2", y2_derived);
        }
    }

    mat->SetMaterialPropertiesTable(mpt);
}

void Utils::AddSurfaceConstIfPresent(G4MaterialPropertiesTable* mpt,
                                     const ConstMap& c,
                                     const std::string& key) {
    auto it = c.find(key);
    if (it != c.end()) {
        mpt->AddConstProperty(key.c_str(), it->second);
    }
}

void Utils::ApplySurface(G4OpticalSurface* surf,
                         const ConstMap& c,
                         const Table* reflectivity,
                         const Table* efficiency) {
    if (!surf) return;

    if (auto it = c.find("SIGMA_ALPHA"); it != c.end()) surf->SetSigmaAlpha(it->second);

    auto* mpt = surf->GetMaterialPropertiesTable();
    if (!mpt) mpt = new G4MaterialPropertiesTable();

    const G4double Emin = GetOr(c, "E_MIN_eV", 2.0) * eV / eV;

    const G4double Emax = GetOr(c, "E_MAX_eV", 3.4) * eV / eV;

    if (reflectivity) {
        mpt->AddProperty("REFLECTIVITY", reflectivity->E, reflectivity->V, reflectivity->E.size());
    } else if (auto it = c.find("REFLECTIVITY"); it != c.end()) {
        auto t = MakeConstantTable(Emin, Emax, it->second);
        mpt->AddProperty("REFLECTIVITY", t.E.data(), t.V.data(), (G4int)t.E.size());
    }

    if (efficiency) {
        mpt->AddProperty("EFFICIENCY", efficiency->E, efficiency->V, efficiency->E.size());
    } else if (auto it = c.find("EFFICIENCY"); it != c.end()) {
        auto t = MakeConstantTable(Emin, Emax, it->second);
        mpt->AddProperty("EFFICIENCY", t.E.data(), t.V.data(), (G4int)t.E.size());
    }

    // Unified model lobe constants (optional)
    AddSurfaceConstIfPresent(mpt, c, "SPECULARLOBECONSTANT");
    AddSurfaceConstIfPresent(mpt, c, "SPECULARSPIKECONSTANT");
    AddSurfaceConstIfPresent(mpt, c, "BACKSCATTERCONSTANT");
    AddSurfaceConstIfPresent(mpt, c, "DIFFUSELOBECONSTANT");

    surf->SetMaterialPropertiesTable(mpt);
}

void Utils::ApplyBirksIfPresent(G4Material* mat, const ConstMap& c) {
    if (!mat) return;
    auto it = c.find("BIRKSCONSTANT");
    if (it == c.end()) return;

    mat->GetIonisation()->SetBirksConstant(it->second * (mm / MeV));
}
