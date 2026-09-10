#include "PostProcessing.hh"

namespace fs = std::filesystem;
using namespace Configuration;

PostProcessing::PostProcessing(std::string outputFolderName,
                               double eMinMeV,
                               double eMaxMeV,
                               std::string particle) : outputFolderName(std::move(outputFolderName)),
                                                       eMinMeV(eMinMeV),
                                                       eMaxMeV(eMaxMeV),
                                                       particleName(std::move(particle)) {
    ROOT::EnableImplicitMT();
    gROOT->SetBatch(kTRUE);

    gErrorIgnoreLevel = kWarning;

    OpenRootFile();
    PrepareOutputDirs();
}

PostProcessing::~PostProcessing() = default;

void PostProcessing::OpenRootFile() {
    rootFile.reset(TFile::Open(outputFile.c_str(), "READ"));
    if (!rootFile || rootFile->IsZombie()) {
        throw std::runtime_error("Failed to open ROOT file: " + outputFile);
    }
}

void PostProcessing::PrepareOutputDirs() {
    fs::path rootPath(outputFile.data());
    fs::path rootDir = rootPath.parent_path();

    postProcessingDir = (rootDir / "post_processing").string();
    runDir = (fs::path(postProcessingDir) / outputFolderName).string();

    effectiveAreaDir = (fs::path(runDir) / "effective_area").string();
    sensitivityDir = (fs::path(runDir) / "sensitivity").string();
    csvDir = (fs::path(runDir) / "CSV").string();
    histogramsDir = (fs::path(runDir) / "Histograms").string();

    fs::create_directories(postProcessingDir);

    if (fs::exists(runDir)) {
        fs::remove_all(runDir);
    }

    fs::create_directories(histogramsDir);
}

TH1* PostProcessing::GetHistOrThrow(const std::string& histName) {
    TH1* h = nullptr;
    rootFile->GetObject(histName.c_str(), h);
    if (!h) {
        throw std::runtime_error("Histogram not found: " + histName);
    }
    return h;
}

double PostProcessing::GeomCenter(double eLow, double eHigh) {
    if (eLow <= 0.0 || eHigh <= 0.0) return 0.0;
    return std::sqrt(eLow * eHigh);
}

double PostProcessing::EffAreaErrFromCounts(double n0, double n, double effArea) {
    if (n0 <= 0.0 || n <= 0.0) return 0.0;
    double val = (n0 - n) / (n * n0);
    if (val <= 0.0) return 0.0;
    return effArea * std::sqrt(val);
}

static short GetColorForParticle(const std::string& particleName) {
    if (particleName == "gamma") return kGreen + 2;
    if (particleName == "e-") return kOrange + 7;
    if (particleName == "proton") return kBlue + 1;
    if (particleName == "alpha") return kRed + 1;
    return kBlack;
}

void PostProcessing::SaveHistPng(const std::string& histName,
                                 const std::string& outPngPath,
                                 const std::string& plotTitle,
                                 const std::string& yTitle,
                                 const bool filled,
                                 const bool useTrig) {
    TH1* h = GetHistOrThrow(histName);

    TCanvas canvas("canvas", "canvas", 1920, 1080);
    bool isEffArea = (histName.find("effArea") != std::string::npos);
    if (!isEffArea) {
        canvas.SetLogx(true);
    }

    const short color = GetColorForParticle(particleName);

    h->SetTitle(plotTitle.c_str());
    h->GetXaxis()->SetTitle("Energy");
    h->GetYaxis()->SetTitle(yTitle.c_str());

    if (eMinMeV > eTriggerThreshold && eMaxMeV > eMinMeV) {
        h->GetXaxis()->SetRangeUser(eMinMeV, eMaxMeV);
    } else {
        h->GetXaxis()->SetRangeUser(eTriggerThreshold, eMaxMeV);
    }

    h->SetLineColor(color);
    h->SetLineWidth(filled ? 1 : 2);

    if (filled) {
        h->SetFillColorAlpha(color, 0.35);
        h->SetFillStyle(1001);
        h->SetMarkerColor(color);
    } else {
        h->SetFillStyle(0);
    }

    h->Draw("HIST");
    canvas.SaveAs(outPngPath.c_str());
}


void PostProcessing::ExportTreeToCsv(const std::string& treeName,
                                     const std::string& csvPath) {
    TTree* tree = nullptr;
    rootFile->GetObject(treeName.c_str(), tree);
    if (!tree) {
        throw std::runtime_error("TTree/NTuple not found: " + treeName);
    }

    auto* leaves = tree->GetListOfLeaves();
    if (!leaves || leaves->GetEntries() == 0) {
        throw std::runtime_error("No leaves found in tree: " + treeName);
    }

    std::ofstream out(csvPath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: " + csvPath);
    }

    const int nLeaves = leaves->GetEntries();
    for (int i = 0; i < nLeaves; ++i) {
        auto* leaf = dynamic_cast<TLeaf*>(leaves->At(i));
        out << (leaf ? leaf->GetName() : "unknown");
        if (i + 1 != nLeaves) out << ",";
    }
    out << "\n";

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        for (int i = 0; i < nLeaves; ++i) {
            auto* leaf = dynamic_cast<TLeaf*>(leaves->At(i));
            if (!leaf) {
                out << "";
                if (i + 1 != nLeaves) out << ",";
                continue;
            }

            if (leaf->InheritsFrom(TLeafC::Class())) {
                auto* lc = dynamic_cast<TLeafC*>(leaf);
                std::string s = lc->GetValueString();

                bool needQuotes = s.find(',') != std::string::npos ||
                    s.find('"') != std::string::npos ||
                    s.find('\n') != std::string::npos;
                if (needQuotes) {
                    std::string esc;
                    esc.reserve(s.size() + 8);
                    for (char c : s) {
                        if (c == '"') esc += "\"\"";
                        else esc += c;
                    }
                    out << "\"" << esc << "\"";
                } else {
                    out << s;
                }
            } else {
                int nData = leaf->GetNdata();
                if (nData <= 1) {
                    double v = leaf->GetValue(0);
                    out << std::setprecision(17) << v;
                } else {
                    std::ostringstream cell;
                    cell << std::setprecision(17);
                    for (int k = 0; k < nData; ++k) {
                        if (k) cell << ";";
                        cell << leaf->GetValue(k);
                    }
                    out << "\"" << cell.str() << "\"";
                }
            }

            if (i + 1 != nLeaves) out << ",";
        }
        out << "\n";
    }

    out.close();
}

void PostProcessing::ExtractNtData() {
    fs::create_directories(csvDir);

    ExportTreeToCsv("edep", (fs::path(csvDir) / "edep.csv").string());
    ExportTreeToCsv("primary", (fs::path(csvDir) / "primary.csv").string());

    if (saveSecondaries) {
        ExportTreeToCsv("event", (fs::path(csvDir) / "event.csv").string());
        ExportTreeToCsv("interactions", (fs::path(csvDir) / "interactions.csv").string());
    }

    if (useOptics) {
        ExportTreeToCsv("sipm_event", (fs::path(csvDir) / "sipm_event.csv").string());
        ExportTreeToCsv("sipm_ch", (fs::path(csvDir) / "sipm_ch.csv").string());
    }

    if (savePhotons) {
        ExportTreeToCsv("photons_count", (fs::path(csvDir) / "photons_count.csv").string());
        ExportTreeToCsv("photons", (fs::path(csvDir) / "photons.csv").string());
    }
}

void PostProcessing::SaveEffArea() {
    fs::create_directories(effectiveAreaDir);

    TH1* gen = GetHistOrThrow("genEnergyHist");
    TH1* trig = GetHistOrThrow("trigEnergyHist");
    TH1* trigOpt = trig;
    TH1* effArea = GetHistOrThrow("effAreaHist");
    TH1* effAreaOpt = effArea;
    if (useOptics) {
        trigOpt = GetHistOrThrow("trigOptEnergyHist");
        effAreaOpt = GetHistOrThrow("effAreaOptHist");
    }

    if (trig->GetNbinsX() != nBins || effArea->GetNbinsX() != nBins) {
        throw std::runtime_error("Histogram binning mismatch among genEnergyHist/trigEnergyHist/effAreaHist");
    }

    SaveHistPng("effAreaHist", (fs::path(effectiveAreaDir) / "effective_area.png").string(),
                "Effective Area vs Energy", "Effective Area [cm^{2}]", false, true);

    SaveHistPng("genEnergyHist", (fs::path(histogramsDir) / "genEnergyHist.png").string(),
                "Initial Energy", "Counts", true, false);

    SaveHistPng("trigEnergyHist", (fs::path(histogramsDir) / "trigEnergyHist.png").string(),
                "N_{trig} vs Energy", "Counts", true, false);

    if (useOptics) {
        SaveHistPng("effAreaOptHist", (fs::path(effectiveAreaDir) / "effective_area_opt.png").string(),
                    "Effective Area vs Energy", "Effective Area [cm^{2}]", false, true);
        SaveHistPng("trigOptEnergyHist", (fs::path(histogramsDir) / "trigOptEnergyHist.png").string(),
                    "N_{trig, opt} vs Energy", "Counts", true, false);
    }

    std::ofstream out((fs::path(effectiveAreaDir) / "effective_area_by_energy.csv").string());
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: effective_area_by_energy.csv");
    }

    out <<
        "E_low,E_high,E_width,E_center_geom,N0_i,N_i,N_i_opt,effective_area,effective_area_opt,effective_area_err,effective_area_opt_err\n";
    out << std::setprecision(17);

    auto* ax = gen->GetXaxis();

    for (int i = 1; i <= nBins; ++i) {
        double eLow = ax->GetBinLowEdge(i);
        double eHigh = ax->GetBinUpEdge(i);
        double eWidth = eHigh - eLow;
        double eCenterGeom = GeomCenter(eLow, eHigh);

        double n0 = gen->GetBinContent(i);
        double n = trig->GetBinContent(i);
        double nOpt = trigOpt->GetBinContent(i);

        double aeff = effArea->GetBinContent(i);
        double aeffOpt = effAreaOpt->GetBinContent(i);
        double aeffErr = EffAreaErrFromCounts(n0, n, aeff);
        double aeffErrOpt = EffAreaErrFromCounts(n0, nOpt, aeffOpt);

        if (!useOptics) {
            nOpt = aeffOpt = aeffErrOpt = 0;
        }

        out << eLow << ","
            << eHigh << ","
            << eWidth << ","
            << eCenterGeom << ","
            << n0 << ","
            << n << ","
            << nOpt << ","
            << aeff << ","
            << aeffOpt << ","
            << aeffErr << ","
            << aeffErrOpt << "\n";
    }

    out.close();
}

void PostProcessing::SaveSensitivity() {
    fs::create_directories(sensitivityDir);

    TH1* gen = GetHistOrThrow("genEnergyHist");
    TH1* trig = GetHistOrThrow("trigEnergyHist");
    TH1* sens = GetHistOrThrow("sensitivityHist");
    TH1* trigOpt = trig;
    TH1* sensOpt = sens;
    if (useOptics) {
        trigOpt = GetHistOrThrow("trigOptEnergyHist");
        sensOpt = GetHistOrThrow("sensitivityOptHist");
    }

    if (trig->GetNbinsX() != nBins || sens->GetNbinsX() != nBins) {
        throw std::runtime_error("Histogram binning mismatch among genEnergyHist/trigEnergyHist/sensitivityHist");
    }

    SaveHistPng("sensitivityHist", (fs::path(sensitivityDir) / "sensitivity.png").string(),
                "Sensitivity vs Energy", "Sensitivity [cm^{2} \\cdot sr]", false, true);

    SaveHistPng("genEnergyHist", (fs::path(histogramsDir) / "genEnergyHist.png").string(),
                "Initial Energy", "Counts", true, false);

    SaveHistPng("trigEnergyHist", (fs::path(histogramsDir) / "trigEnergyHist.png").string(),
                "N_{trig} vs Energy", "Counts", true, false);

    if (useOptics) {
        SaveHistPng("sensitivityOptHist", (fs::path(sensitivityDir) / "sensitivity_opt.png").string(),
                    "Sensitivity vs Energy", "Sensitivity [cm^{2} \\cdot sr]", false, true);
        SaveHistPng("trigOptEnergyHist", (fs::path(histogramsDir) / "trigOptEnergyHist.png").string(),
                    "N_{trig, opt} vs Energy", "Counts", true, false);
    }

    std::ofstream out((fs::path(sensitivityDir) / "sensitivity_by_energy.csv").string());
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: sensitivity_by_energy.csv");
    }

    out << "E_low,E_high,E_width,E_center_geom,N0_i,N_i,N_i_opt,sensitivity,sensitivity_opt\n";
    out << std::setprecision(17);

    auto* ax = gen->GetXaxis();

    for (int i = 1; i <= nBins; ++i) {
        double eLow = ax->GetBinLowEdge(i);
        double eHigh = ax->GetBinUpEdge(i);
        double eWidth = eHigh - eLow;
        double eCenterGeom = GeomCenter(eLow, eHigh);

        double n0 = gen->GetBinContent(i);
        double n = trig->GetBinContent(i);
        double nOpt = trigOpt->GetBinContent(i);

        double s_ = sens->GetBinContent(i);
        double sOpt = sensOpt->GetBinContent(i);
        if (!useOptics) {
            nOpt = sOpt = 0;
        }

        out << eLow << ","
            << eHigh << ","
            << eWidth << ","
            << eCenterGeom << ","
            << n0 << ","
            << n << ","
            << nOpt << ","
            << s_ << ","
            << sOpt << "\n";
    }

    out.close();
}

void PostProcessing::SaveTrigEdepCsv() {
    TTree* primary = nullptr;
    rootFile->GetObject("primary", primary);
    if (!primary) {
        throw std::runtime_error("TTree not found: primary");
    }

    Int_t eventID_p = 0;
    double E0 = 0.0;

    primary->SetBranchStatus("*", false);
    primary->SetBranchStatus("eventID", true);
    primary->SetBranchStatus("E_MeV", true);

    primary->SetBranchAddress("eventID", &eventID_p);
    primary->SetBranchAddress("E_MeV", &E0);

    std::unordered_map<int, double> e0ByEvent;
    e0ByEvent.reserve(std::max<Long64_t>(1, primary->GetEntries()));

    const Long64_t nP = primary->GetEntries();
    for (Long64_t i = 0; i < nP; ++i) {
        primary->GetEntry(i);
        e0ByEvent[eventID_p] = E0;
    }

    TTree* edep = nullptr;
    rootFile->GetObject("edep", edep);
    if (!edep) {
        throw std::runtime_error("TTree not found: edep");
    }

    Int_t eventID_e = 0;
    Char_t det_name[64] = {};
    double edep_MeV = 0.0;

    edep->SetBranchStatus("*", false);
    edep->SetBranchStatus("eventID", true);
    edep->SetBranchStatus("det_name", true);
    edep->SetBranchStatus("edep_MeV", true);

    edep->SetBranchAddress("eventID", &eventID_e);
    edep->SetBranchAddress("det_name", det_name);
    edep->SetBranchAddress("edep_MeV", &edep_MeV);

    struct Agg {
        double trigger = 0.0;
        double sideVeto = 0.0;
        double bottomVeto = 0.0;
    };

    std::unordered_map<int, Agg> agg;
    agg.reserve(std::max<Long64_t>(1, edep->GetEntries()));

    const Long64_t nE = edep->GetEntries();
    for (Long64_t i = 0; i < nE; ++i) {
        edep->GetEntry(i);

        auto& a = agg[eventID_e];

        if (std::strcmp(det_name, "Trigger") == 0) {
            a.trigger += edep_MeV;
        } else if (std::strcmp(det_name, "SideVeto") == 0) {
            a.sideVeto += edep_MeV;
        } else if (std::strcmp(det_name, "BottomVeto") == 0) {
            a.bottomVeto += edep_MeV;
        }
    }

    const std::string outPath = (fs::path(histogramsDir) / "trig_edep.csv").string();
    std::ofstream out(outPath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: " + outPath);
    }

    out << "eventID,E0,Trigger_only_edep\n";
    out << std::setprecision(17);

    std::vector<int> events;
    events.reserve(e0ByEvent.size());
    for (const auto& kv : e0ByEvent) events.push_back(kv.first);
    std::sort(events.begin(), events.end());

    for (int evt : events) {
        const double e0 = e0ByEvent.at(evt);

        double trigger_only = 0.0;
        auto it = agg.find(evt);
        if (it != agg.end()) {
            const double t = it->second.trigger;
            const double sv = it->second.sideVeto;
            const double bv = it->second.bottomVeto;

            if (t > 0.0 && sv + bv == 0.0) {
                trigger_only = t;
            }
        }

        out << evt << "," << e0 << "," << trigger_only << "\n";
    }

    out.close();
}


void PostProcessing::SaveEdepCsv() {
    TTree* edep = nullptr;
    rootFile->GetObject("edep", edep);
    if (!edep) {
        throw std::runtime_error("TTree not found: edep");
    }

    Int_t eventID_e = 0;
    Char_t det_name[64] = {};
    double edep_MeV = 0.0;

    edep->SetBranchStatus("*", false);
    edep->SetBranchStatus("eventID", true);
    edep->SetBranchStatus("det_name", true);
    edep->SetBranchStatus("edep_MeV", true);

    edep->SetBranchAddress("eventID", &eventID_e);
    edep->SetBranchAddress("det_name", det_name);
    edep->SetBranchAddress("edep_MeV", &edep_MeV);

    struct DetectorEdep {
        double trigger = 0.0;
        double sideVeto = 0.0;
        double bottomVeto = 0.0;
    };

    std::unordered_map<int, DetectorEdep> edepMap;
    edepMap.reserve(std::max<Long64_t>(1, edep->GetEntries()));

    const Long64_t nEntries = edep->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        edep->GetEntry(i);

        auto& deps = edepMap[eventID_e];

        if (std::strcmp(det_name, "Trigger") == 0) {
            deps.trigger += edep_MeV;
        } else if (std::strcmp(det_name, "SideVeto") == 0) {
            deps.sideVeto += edep_MeV;
        } else if (std::strcmp(det_name, "BottomVeto") == 0) {
            deps.bottomVeto += edep_MeV;
        }
    }

    const std::string outPath = (fs::path(histogramsDir) / "edep.csv").string();
    std::ofstream out(outPath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: " + outPath);
    }

    out << "eventID,Trigger,Trigger_edep_MeV,SideVeto_edep_MeV,BottomVeto_edep_MeV\n";
    out << std::setprecision(17);

    std::vector<int> eventIDs;
    eventIDs.reserve(edepMap.size());
    for (const auto& kv : edepMap) {
        eventIDs.push_back(kv.first);
    }
    std::sort(eventIDs.begin(), eventIDs.end());

    for (int evtID : eventIDs) {
        const auto& deps = edepMap.at(evtID);

        int trig = deps.trigger > 0.0 &&
                   deps.sideVeto == 0.0 &&
                   deps.bottomVeto == 0.0
                       ? 1
                       : 0;

        out << evtID << ","
            << trig << ","
            << deps.trigger << ","
            << deps.sideVeto << ","
            << deps.bottomVeto << "\n";
    }

    out.close();
}

void PostProcessing::SaveOpticsCsv() {
    std::string opticDir = (fs::path(runDir) / "optic").string();
    fs::create_directories(opticDir);

    std::unordered_map<int, int> edepTriggerMap;

    TTree* edep = nullptr;
    rootFile->GetObject("edep", edep);
    if (edep) {
        Int_t eventID_e = 0;
        Char_t det_name[64] = {};
        double edep_MeV = 0.0;

        edep->SetBranchStatus("*", true);

        edep->SetBranchAddress("eventID", &eventID_e);
        edep->SetBranchAddress("det_name", det_name);
        edep->SetBranchAddress("edep_MeV", &edep_MeV);

        struct DetectorEdep {
            double trigger = 0.0;
            double sideVeto = 0.0;
            double bottomVeto = 0.0;
        };

        std::unordered_map<int, DetectorEdep> edepMap;
        const Long64_t nEntries = edep->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            edep->GetEntry(i);

            auto& deps = edepMap[eventID_e];

            if (std::strcmp(det_name, "Trigger") ==0 ) {
                deps.trigger += edep_MeV;
            } else if (std::strcmp(det_name, "SideVeto") == 0) {
                deps.sideVeto += edep_MeV;
            } else if (std::strcmp(det_name, "BottomVeto") == 0) {
                deps.bottomVeto += edep_MeV;
            }
        }

        for (const auto& [evtID, deps] : edepMap) {
            edepTriggerMap[evtID] = deps.trigger > 0.0 &&
                                    deps.sideVeto == 0.0 &&
                                    deps.bottomVeto == 0.0 ? 1 : 0;
        }
    }

    TTree* sipmEvent = nullptr;
    rootFile->GetObject("sipm_event", sipmEvent);
    if (!sipmEvent) {
        throw std::runtime_error("TTree not found: sipm_event");
    }

    Int_t eventID = 0;
    Int_t npe_Trigger = 0;
    Int_t npe_veto = 0;
    Int_t npe_bottom_veto = 0;

    sipmEvent->SetBranchStatus("*", true);

    sipmEvent->SetBranchAddress("eventID", &eventID);
    sipmEvent->SetBranchAddress("npe_Trigger", &npe_Trigger);
    sipmEvent->SetBranchAddress("npe_veto", &npe_veto);
    sipmEvent->SetBranchAddress("npe_bottom_veto", &npe_bottom_veto);

    struct EventInfo {
        Int_t trigger_npe = 0;
        Int_t veto_npe = 0;
        Int_t bottom_veto_npe = 0;
        Int_t trigger = 0;
    };

    std::unordered_map<Int_t, EventInfo> eventMap;
    eventMap.reserve(std::max<Long64_t>(1, sipmEvent->GetEntries()));

    const Long64_t nEvents = sipmEvent->GetEntries();
    for (Long64_t i = 0; i < nEvents; ++i) {
        sipmEvent->GetEntry(i);

        EventInfo info;
        info.trigger_npe = npe_Trigger;
        info.veto_npe = npe_veto;
        info.bottom_veto_npe = npe_bottom_veto;
        info.trigger = npe_Trigger > 0 && npe_veto + npe_bottom_veto == 0 ? 1 : 0;

        eventMap[eventID] = info;
    }

    TTree* sipmCh = nullptr;
    rootFile->GetObject("sipm_ch", sipmCh);
    if (!sipmCh) {
        throw std::runtime_error("TTree not found: sipm_ch");
    }

    Int_t ch_eventID = 0;
    char subdet[64] = {};
    Int_t ch = 0;
    Int_t npe = 0;

    sipmCh->SetBranchStatus("*", true);
    sipmCh->SetBranchAddress("eventID", &ch_eventID);
    sipmCh->SetBranchAddress("subdet", subdet);
    sipmCh->SetBranchAddress("ch", &ch);
    sipmCh->SetBranchAddress("npe", &npe);

    using ChannelMap = std::unordered_map<Int_t, std::unordered_map<Int_t, Int_t>>;
    ChannelMap triggerChannels;
    ChannelMap vetoChannels;
    ChannelMap bottomVetoChannels;

    std::set<Int_t> allTriggerChannels;
    std::set<Int_t> allVetoChannels;
    std::set<Int_t> allBottomVetoChannels;

    const Long64_t nChEntries = sipmCh->GetEntries();
    for (Long64_t i = 0; i < nChEntries; ++i) {
        sipmCh->GetEntry(i);

        std::string subdetStr(subdet);


        if (subdetStr == "Trigger" || subdetStr == "Trigger") {
            triggerChannels[ch_eventID][ch] = npe;
            allTriggerChannels.insert(ch);
        } else if (subdetStr == "SideVeto" || subdetStr == "Veto") {
            vetoChannels[ch_eventID][ch] = npe;
            allVetoChannels.insert(ch);
        } else if (subdetStr == "BottomVeto") {
            bottomVetoChannels[ch_eventID][ch] = npe;
            allBottomVetoChannels.insert(ch);
        }
    }

    std::vector sortedTriggerChannels(allTriggerChannels.begin(), allTriggerChannels.end());
    std::vector sortedVetoChannels(allVetoChannels.begin(), allVetoChannels.end());
    std::vector sortedBottomVetoChannels(allBottomVetoChannels.begin(), allBottomVetoChannels.end());

    std::sort(sortedTriggerChannels.begin(), sortedTriggerChannels.end());
    std::sort(sortedVetoChannels.begin(), sortedVetoChannels.end());
    std::sort(sortedBottomVetoChannels.begin(), sortedBottomVetoChannels.end());

    std::ofstream trigOptFile(fs::path(opticDir) / "trig_opt.csv");
    if (!trigOptFile.is_open()) {
        throw std::runtime_error("Cannot open output CSV: trig_opt.csv");
    }

    trigOptFile << "eventID,trigger_opt,trigger_edep,Trigger_npe,Veto_npe,BottomVeto_npe\n";

    std::vector<Int_t> allEventIDs;
    allEventIDs.reserve(eventMap.size());
    for (const auto& kv : eventMap) {
        allEventIDs.push_back(kv.first);
    }
    std::sort(allEventIDs.begin(), allEventIDs.end());

    for (Int_t evtID : allEventIDs) {
        const auto& info = eventMap[evtID];

        int trigger_edep = 0;
        auto it = edepTriggerMap.find(evtID);
        if (it != edepTriggerMap.end()) {
            trigger_edep = it->second;
        }

        trigOptFile << evtID << ","
            << info.trigger << ","
            << trigger_edep << ","
            << info.trigger_npe << ","
            << info.veto_npe << ","
            << info.bottom_veto_npe << "\n";
    }
    trigOptFile.close();

    if (!sortedTriggerChannels.empty()) {
        std::ofstream triggerFile(fs::path(opticDir) / "Trigger_channel.csv");
        if (!triggerFile.is_open()) {
            throw std::runtime_error("Cannot open output CSV: Trigger_channel.csv");
        }

        triggerFile << "eventID";
        for (Int_t ch : sortedTriggerChannels) {
            triggerFile << ",ch" << ch;
        }
        triggerFile << "\n";

        for (Int_t evtID : allEventIDs) {
            triggerFile << evtID;

            const auto& channels = triggerChannels[evtID];
            for (Int_t ch : sortedTriggerChannels) {
                auto it = channels.find(ch);
                if (it != channels.end()) {
                    triggerFile << "," << it->second;
                } else {
                    triggerFile << ",0";
                }
            }
            triggerFile << "\n";
        }
        triggerFile.close();
    }

    if (!sortedVetoChannels.empty()) {
        std::ofstream vetoFile(fs::path(opticDir) / "Veto_channel.csv");
        if (!vetoFile.is_open()) {
            throw std::runtime_error("Cannot open output CSV: Veto_channel.csv");
        }

        vetoFile << "eventID";
        for (Int_t ch : sortedVetoChannels) {
            vetoFile << ",ch" << ch;
        }
        vetoFile << "\n";

        for (Int_t evtID : allEventIDs) {
            vetoFile << evtID;

            const auto& channels = vetoChannels[evtID];
            for (Int_t ch : sortedVetoChannels) {
                auto it = channels.find(ch);
                if (it != channels.end()) {
                    vetoFile << "," << it->second;
                } else {
                    vetoFile << ",0";
                }
            }
            vetoFile << "\n";
        }
        vetoFile.close();
    }

    if (!sortedBottomVetoChannels.empty()) {
        std::ofstream bottomFile(fs::path(opticDir) / "BottomVeto_channel.csv");
        if (!bottomFile.is_open()) {
            throw std::runtime_error("Cannot open output CSV: BottomVeto_channel.csv");
        }

        bottomFile << "eventID";
        for (Int_t ch : sortedBottomVetoChannels) {
            bottomFile << ",ch" << ch;
        }
        bottomFile << "\n";

        for (Int_t evtID : allEventIDs) {
            bottomFile << evtID;

            const auto& channels = bottomVetoChannels[evtID];
            for (Int_t ch : sortedBottomVetoChannels) {
                auto it = channels.find(ch);
                if (it != channels.end()) {
                    bottomFile << "," << it->second;
                } else {
                    bottomFile << ",0";
                }
            }
            bottomFile << "\n";
        }
        bottomFile.close();
    }
}
