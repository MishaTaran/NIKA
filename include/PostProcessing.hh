#ifndef POSTPROCESSING_HH
#define POSTPROCESSING_HH

#include <string>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>
#include <sstream>
#include <filesystem>

#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TLeafC.h>
#include <TH1.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TError.h>

#include "Configuration.hh"

class TFile;
class TH1;

class PostProcessing {
public:
    PostProcessing(std::string outputFolderName,
                   double eMinMeV,
                   double eMaxMeV,
                   std::string particle);

    ~PostProcessing();

    void ExtractNtData();
    void SaveEffArea();
    void SaveSensitivity();

    void SaveTrigEdepCsv();
    void SaveEdepCsv();

    void SaveOpticsCsv();

private:
    std::string outputFolderName;
    double eMinMeV;
    double eMaxMeV;
    std::string particleName;

    std::unique_ptr<TFile> rootFile;

    std::string postProcessingDir;
    std::string runDir;
    std::string effectiveAreaDir;
    std::string sensitivityDir;
    std::string csvDir;
    std::string histogramsDir;

    void OpenRootFile();
    void PrepareOutputDirs();

    void ExportTreeToCsv(const std::string& treeName,
                         const std::string& csvPath);

    TH1* GetHistOrThrow(const std::string& histName);

    void SaveHistPng(const std::string& histName,
                     const std::string& outPngPath,
                     const std::string& plotTitle,
                     const std::string& yTitle,
                     bool filled,
                     bool useTrig = false);

    static double GeomCenter(double eLow, double eHigh);
    static double EffAreaErrFromCounts(double n0, double n, double effArea);
};


#endif //POSTPROCESSING_HH
