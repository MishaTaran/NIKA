#ifndef UTILS_HH
#define UTILS_HH

#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4OpticalSurface.hh>
#include <G4SystemOfUnits.hh>
#include <G4Types.hh>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

class Utils {
public:
    using ConstMap = std::unordered_map<std::string, G4double>;

    static ConstMap ReadConstFile(const std::string& filename);

    static G4double GetRequired(const ConstMap& m, const std::string& key);

    static G4double GetOr(const ConstMap& m, const std::string& key, G4double def);

    struct Table {
        std::vector<G4double> E;
        std::vector<G4double> V;
    };

    struct EmissionTables {
        Table c1;
        Table c2;
    };

    static Table ReadCSV(const std::string& filename,
                         G4double valueScale = 1.0,
                         bool clampNonNegative = false);

    static EmissionTables ReadEmissionCSV(const std::string& filename,
                                          G4double valueScale = 1.0,
                                          bool clampNonNegative = true);

    static void NormalizeMaxToOne(Table& t);
    static void ClampNonNegative(Table& t);

    static void ApplyMaterialTable(G4Material* mat,
                                   const Table& rindex,
                                   const Table* abslength = nullptr);

    static void ApplyScintillation(G4Material* mat,
                                   G4MaterialPropertiesTable* mpt,
                                   const ConstMap& c,
                                   const Table& scintComponent1,
                                   const Table& scintComponent2,
                                   bool requireYield = true);

    static void ApplySurface(G4OpticalSurface* surf,
                             const ConstMap& c,
                             const Table* reflectivity = nullptr,
                             const Table* efficiency = nullptr);

    static Table MakeConstantTable(G4double Emin_eV, G4double Emax_eV, G4double value);

    static void ApplyBirksIfPresent(G4Material* mat, const ConstMap& c);

private:
    static std::string Trim(std::string s_);
    static bool StartsWith(const std::string& s_, const char* prefix);

    static G4double UnitFactor(const std::string& unitToken);

    static void AddConstIfPresent(G4MaterialPropertiesTable* mpt,
                                  const ConstMap& c,
                                  const std::string& key);

    static void AddSurfaceConstIfPresent(G4MaterialPropertiesTable* mpt,
                                         const ConstMap& c,
                                         const std::string& key);

    static void ValidateSameSize(const Table& t, const std::string& filename);
};

#endif //UTILS_HH