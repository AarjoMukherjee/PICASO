#include "Common.h"
#include "RootOutputManager.h"

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCutG.h>
#include <TObject.h>

#include <array>
#include <cstdio>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <string>
#include <utility>

namespace {
constexpr int kChannels = 96;
constexpr int kEnergyPairs = 6;
constexpr int kClovers = 18;

constexpr int pairCh1[kEnergyPairs] = {80, 82, 84, 86, 88, 90};
constexpr int pairCh2[kEnergyPairs] = {81, 83, 85, 87, 89, 91};

const char* cloverNames[kClovers] = {
    "Clover1_Dummy", "Clover2", "Clover3", "Clover4", "Clover5",
    "Clover6", "Clover7", "Clover8", "Clover9", "Clover11",
    "Clover12_Dummy", "Clover14", "Clover15", "Clover16",
    "Clover17", "Clover19", "Clover20", "Clover21"
};

std::unique_ptr<TCutG> makeGate(
    const char* name, const char* varX, const char* varY,
    std::initializer_list<std::pair<double,double>> points)
{
    auto gate = std::make_unique<TCutG>(name, static_cast<int>(points.size()));
    gate->SetVarX(varX);
    gate->SetVarY(varY);
    int i = 0;
    for (const auto& point : points)
        gate->SetPoint(i++, point.first, point.second);
    return gate;
}

} // namespace

struct RootOutputManager::Impl {
    struct ThreadData {
        TFile* file = nullptr;
        TTree* tree = nullptr;
        TTree* traceTree = nullptr;

        int run = 0, channel = 0, finish_code = 0, cfd_source = 0;
        unsigned long long coincID = 0, coarseT = 0;
        double fineT = 0.0, energy = 0.0;

        int trace_channel = 0, trace_length = 0, trace_run = 0;
        double trace_energy = 0.0, trace_fineT = 0.0;
        unsigned long long trace_coincID = 0, trace_coarseT = 0;
        std::vector<uint16_t> trace_data;

        TH2D* hE1vsE2[kEnergyPairs] = {};
        TH1D* hAddback[kClovers] = {};
        TH2D* hMatrix = nullptr;
        TH1D* hAlphaGamma = nullptr;
        TH1D* hProtonGamma = nullptr;

        // Each OpenMP thread owns its own gate objects. They are read-only
        // during filling and are destroyed with ThreadData.
        std::array<std::unique_ptr<TCutG>, kEnergyPairs> alphaGates;
        std::array<std::unique_ptr<TCutG>, kEnergyPairs> protonGates;
    };

    std::vector<ThreadData> threads;
    bool writeAnalysis = false;
    bool writeTree = false;
    bool saveTrace = false;
    PixieType type = PixieType::PIXIE_100MHZ;

    void createGates(ThreadData& d, int tid);
    void createHistograms(ThreadData& d);
    void fillRootHistograms(ThreadData& d, const CoincidenceBlock& block);
};

void RootOutputManager::Impl::createGates(ThreadData& d, int tid)
{
    const std::string alphaName1 = "Alpha_82_83_thread_" + std::to_string(tid);
    d.alphaGates[1] = makeGate(alphaName1.c_str(), "Channels 82 and 83", "E1 vs E2", {
        {273.31396414346477, 1867.0539312700876},
        {363.74874809974381, 1849.7254996393247},
        {625.00479064010551, 1681.3921637976291},
        {1484.135238224756, 1339.7745116483056},
        {2800.4637602550392, 1020.4362715956772},
        {4393.1207888183972, 807.54411156059155},
        {6618.8213050757086, 644.16175618482816},
        {6608.7729957472329, 597.12744175847206},
        {4910.6087192348832, 663.96567804855704},
        {3272.7342986933854, 800.11764086169319},
        {2222.6859738677013, 926.36764274296502},
        {1448.9661555750918, 1146.6862734769488},
        {851.09175053080298, 1379.3823553757634},
        {343.65212944279301, 1661.5882419339002},
        {258.24150015075156, 1738.3284391558495},
        {218.04826283684997, 1820.0196168437315},
        {258.24150015075156, 1864.5784410371214},
        {273.31396414346477, 1867.0539312700876}
    });
    const std::string alphaName2 = "Alpha_84_85_thread_" + std::to_string(tid);
    d.alphaGates[2] = makeGate(alphaName2.c_str(), "Channels 84 and 85", "E1 vs E2", {
        {512.07243426489492, 1582.6139766174535},
        {637.83330570409544, 1561.8235351311812},
        {766.81881487250621, 1485.058828104945},
        {1305.333315650621, 1272.3566190530821},
        {2130.8405743284502, 1035.6654390555204},
        {3095.0072553623208, 861.34558351677549},
        {3749.6087143920054, 784.58087649053914},
        {4446.1304639014234, 730.2058756802885},
        {4500.9493052979979, 671.03308068089814},
        {3943.0869781446213, 687.02572797803077},
        {3188.5217495094189, 754.19484662598734},
        {2259.826083496861, 894.93014284075389},
        {1814.8260768658438, 990.88602662354913},
        {1244.0651987956262, 1149.2132348651614},
        {902.25359949933761, 1283.5514721610748},
        {637.83330570409544, 1408.2941210787087},
        {508.84779653568467, 1493.0551517535114},
        {476.60141924358197, 1560.2242704014679},
        {518.52170972331544, 1582.6139766174535},
        {512.07243426489492, 1582.6139766174535}
    });
    const std::string alphaName3 = "Alpha_86_87_thread_" + std::to_string(tid);
    d.alphaGates[3] = makeGate(alphaName3.c_str(), "Channels 86 and 87", "E1 vs E2", {
        {269.26262804069165, 2079.838245752323},
        {280.63340115698929, 1863.6924091981655},
        {820.74512418112681, 1571.2598068013642},
        {1508.6768977171337, 1250.8553902622602},
        {2139.7548056716523, 1083.0245054084439},
        {3038.0458818591651, 935.53675811266589},
        {4112.5839413492922, 816.02082495919058},
        {6295.7723796784376, 671.07596985816735},
        {7404.4227585174558, 630.3896947420908},
        {7375.9958257267117, 686.33332302669623},
        {6051.3007576780383, 767.70587325884958},
        {4669.7518240478748, 897.39337519134392},
        {3225.6636382780757, 1100.8247507717274},
        {2230.7209906020335, 1304.2561263521109},
        {1491.6207380426872, 1545.8308848538163},
        {877.59898976261502, 1792.4914277450312},
        {337.4872667384775, 2084.9240301418326},
        {274.94801459884047, 2074.7524613628134},
        {269.26262804069165, 2079.838245752323}
    });
    const std::string alphaName4 = "Alpha_88_89_thread_" + std::to_string(tid);
    d.alphaGates[4] = makeGate(alphaName4.c_str(), "Channels 88 and 89", "E1 vs E2", {
        {174.01384555651657, 1762.8235389309777},
        {304.29162527558469, 1729.5514796116561},
        {695.12496443278883, 1501.0833389523132},
        {1472.1388649000876, 1170.5808830470505},
        {2179.3610976607433, 975.38480170702951},
        {3156.4444455537532, 784.62499494291819},
        {4682.5555794056936, 620.48283563426412},
        {6013.2500436790315, 514.01224581243468},
        {5994.6389322905934, 440.81371530992669},
        {4691.8611350999126, 518.44852038834415},
        {3500.7500062398622, 629.35538478608328},
        {2416.6527678633311, 797.93381867064682},
        {1686.1666458671282, 959.85784069134615},
        {1323.2499737925812, 1075.200979664995},
        {964.98607956514411, 1199.4166677904627},
        {662.55551950302197, 1345.8137287954785},
        {388.04162652355694, 1492.2107898004942},
        {169.36106770940705, 1647.4803999573289},
        {178.6666234036262, 1765.0416762189325},
        {174.01384555651657, 1762.8235389309777}
    });
    const std::string alphaName5 = "Alpha_90_91_thread_" + std::to_string(tid);
    d.alphaGates[5] = makeGate(alphaName5.c_str(), "Channels 90 and 91", "E1 vs E2", {
        {276.22580512125126, 1696.5453528311364},
        {775.59175459136975, 1420.9448585282903},
        {1517.8924902901945, 1145.3443642254445},
        {2800.048306497255, 840.90195772811467},
        {3740.2959050490995, 699.89705366619341},
        {3668.3152276480018, 600.55268944074896},
        {2773.0555524718434, 639.00857236672755},
        {1859.8007079454105, 818.46935935462716},
        {1072.5120488708994, 1052.4093138209964},
        {591.14126875105569, 1289.5539251978639},
        {253.73184343340802, 1510.6752520222403},
        {208.74392005772177, 1600.4056455161901},
        {262.72942810854533, 1693.3406959206382},
        {276.22580512125126, 1696.5453528311364}
    });
    const std::string protonName1 = "Proton_82_83_thread_" + std::to_string(tid);
    d.protonGates[1] = makeGate(protonName1.c_str(), "Channels 82 and 83", "E1 vs E2", {
        {141.44926793333369, 450.82352850016423},
        {327.92270066369531, 391.21568447468326},
        {582.99516823268743, 273.568623898076},
        {573.33333233992255, 246.901956834045},
        {340.48308732428961, 246.901956834045},
        {260.28984941434135, 248.47058430839979},
        {164.63767407596933, 293.96078106468792},
        {120.19322896925101, 367.68627235936174},
        {122.12559614780398, 444.54901860274515},
        {147.24636946899258, 450.82352850016423},
        {141.44926793333369, 450.82352850016423}
    });
    const std::string protonName2 = "Proton_84_85_thread_" + std::to_string(tid);
    d.protonGates[2] = makeGate(protonName2.c_str(), "Channels 84 and 85", "E1 vs E2", {
        {247.33332742253899, 346.19607798725946},
        {301.777772678269, 359.23529386783343},
        {481.63888646951966, 255.85293938613978},
        {456.36110831507358, 223.25489968470481},
        {291.08332807446482, 254.92156682324162},
        {245.38888294912007, 302.42156753104683},
        {250.24999413266741, 344.33333286146319},
        {247.33332742253899, 346.19607798725946}
    });
    const std::string protonName3 = "Proton_86_87_thread_" + std::to_string(tid);
    d.protonGates[3] = makeGate(protonName3.c_str(), "Channels 86 and 87", "E1 vs E2", {
        {179.37800262510171, 463.22794012188467},
        {401.14732960121995, 385.37867425595778},
        {846.64854396917417, 272.14337845097327},
        {1013.4661793051923, 238.17278970947791},
        {938.88888350791376, 222.60293653629253},
        {707.30675445320617, 238.17278970947791},
        {481.61230664565215, 272.14337845097327},
        {316.75723172535191, 323.09926156321637},
        {193.11592553512679, 385.37867425595778},
        {159.75239846792311, 413.68749820720393},
        {155.82727763648739, 451.90441054138626},
        {187.22824428797315, 467.47426371457163},
        {179.37800262510171, 463.22794012188467}
    });
    const std::string protonName4 = "Proton_88_89_thread_" + std::to_string(tid);
    d.protonGates[4] = makeGate(protonName4.c_str(), "Channels 88 and 89", "E1 vs E2", {
        {108.2240271738452, 402.84191070436771},
        {151.71436598615298, 415.91176384030018},
        {325.67572123538412, 324.42279188877296},
        {348.25724331100548, 277.95220296101309},
        {330.69383725218887, 267.78676163306557},
        {261.27656568638992, 264.88234982508061},
        {181.82306208698145, 282.30882067299058},
        {114.07849586011739, 311.35293875284049},
        {94.842384462365885, 336.0404391207129},
        {90.66062111502859, 367.98896900854777},
        {93.169679123430967, 395.58088118440526},
        {114.07849586011739, 405.74632251235266},
        {108.2240271738452, 402.84191070436771}
    });
}

void RootOutputManager::Impl::createHistograms(ThreadData& d)
{
    for (int p = 0; p < kEnergyPairs; ++p) {
        const std::string name = "hE1vsE2_ch" + std::to_string(pairCh1[p]) +
                                 "_ch" + std::to_string(pairCh2[p]);
        const std::string title = "E1 vs E2: Channels " + std::to_string(pairCh1[p]) +
            " and " + std::to_string(pairCh2[p]) + ";Channel " +
            std::to_string(pairCh1[p]) + " energy;Channel " +
            std::to_string(pairCh2[p]) + " energy";
        d.hE1vsE2[p] = new TH2D(name.c_str(), title.c_str(),
                                2500, 0.0, 10000.0, 2500, 0.0, 10000.0);
        d.hE1vsE2[p]->SetDirectory(d.file);
    }

    for (int c = 0; c < kClovers; ++c) {
        const std::string name = std::string("AddBackSpectrum_") + cloverNames[c];
        const std::string title = std::string("Add-back spectrum ") + cloverNames[c] +
                                  ";Summed calibrated energy;Counts";
        d.hAddback[c] = new TH1D(name.c_str(), title.c_str(), 10000, 0.0, 10000.0);
        d.hAddback[c]->SetDirectory(d.file);
    }

    d.hMatrix = new TH2D(
        "Matrix",
        "Gamma-Gamma Matrix (Clover Addback);E_{#gamma 1} (keV);E_{#gamma 2} (keV)",
        4000, 0.0, 12000.0, 4000, 0.0, 12000.0);
    d.hMatrix->SetDirectory(d.file);

    d.hAlphaGamma = new TH1D(
        "AlphaGamma",
        "Gamma Spectrum with Alpha Gate 90Degree;Gamma Energy (keV);Counts",
        10000, 0.0, 10000.0);
    d.hAlphaGamma->SetDirectory(d.file);

    d.hProtonGamma = new TH1D(
        "ProtonGamma",
        "Gamma Spectrum with Proton Gate 90Degree;Gamma Energy (keV);Counts",
        10000, 0.0, 10000.0);
    d.hProtonGamma->SetDirectory(d.file);
}

void RootOutputManager::Impl::fillRootHistograms(
    ThreadData& d, const CoincidenceBlock& block)
{
    double channelEnergy[kChannels] = {};
    bool channelValid[kChannels] = {};
    for(const auto& member : block.members) {
        if(member.channel < 0 || member.channel >= kChannels || !member.energyValid)
            continue;
        channelEnergy[member.channel] = member.calibratedEnergy;
        channelValid[member.channel] = true;
    }

    bool alphaEvent = false;
    bool protonEvent = false;

    // 1. Channel-pair matrices and particle gates.
    for(int p = 0; p < kEnergyPairs; ++p) {
        const int c1 = pairCh1[p];
        const int c2 = pairCh2[p];
        if(!channelValid[c1] || !channelValid[c2])
            continue;

        const double e1 = channelEnergy[c1];
        const double e2 = channelEnergy[c2];
        if(e1 > 1.0 && e1 < 40000.0 && e2 > 1.0 && e2 < 40000.0) {
            if(d.alphaGates[p] && d.alphaGates[p]->IsInside(e2, e1))
                alphaEvent = true;
            if(d.protonGates[p] && d.protonGates[p]->IsInside(e2, e1))
                protonEvent = true;
            if(d.hE1vsE2[p])
                d.hE1vsE2[p]->Fill(e2, e1);
        }
    }

    // 2. Fill addback spectra from the values calculated by CoincidenceProcessor.
    for(int c = 0; c < kClovers; ++c) {
        if(block.addbackValid[c] && d.hAddback[c])
            d.hAddback[c]->Fill(block.addbackEnergy[c]);
    }

    // 3. Preserve the original 7Li.marcos particle-gated gamma condition.
    const bool requiredCloverPresent =
        block.addbackValid[9] || block.addbackValid[11] || block.addbackValid[12];
    if(requiredCloverPresent) {
        for(int c = 0; c < kClovers; ++c) {
            if(!block.addbackValid[c])
                continue;
            const double gammaEnergy = block.addbackEnergy[c];
            if(gammaEnergy <= 1.0 || gammaEnergy >= 10000.0)
                continue;
            if(alphaEvent && d.hAlphaGamma)
                d.hAlphaGamma->Fill(gammaEnergy);
            if(protonEvent && d.hProtonGamma)
                d.hProtonGamma->Fill(gammaEnergy);
        }
    }

    // 4. Gamma-gamma matrix: each distinct Clover pair is filled once.
    // if(!d.hMatrix)
    //     return;
    // for(int i = 0; i < kClovers; ++i) {
    //     if(!block.addbackValid[i])
    //         continue;
    //     const double e1 = block.addbackEnergy[i];
    //     if(e1 <= 1.0 || e1 >= 10000.0)
    //         continue;
    //     for(int j = i + 1; j < kClovers; ++j) {
    //         if(!block.addbackValid[j])
    //             continue;
    //         const double e2 = block.addbackEnergy[j];
    //         if(e2 <= 1.0 || e2 >= 10000.0)
    //             continue;
    //         d.hMatrix->Fill(e1, e2);
    //     }
    // }
    if(!d.hMatrix)
        return;

    const int selectedChannels[] = {9, 11, 12};

    for(int i = 0; i < 3; ++i) {
        const int ch1 = selectedChannels[i];

        if(!block.addbackValid[ch1])
            continue;

        const double e1 = block.addbackEnergy[ch1];
        if(e1 <= 1.0 || e1 >= 10000.0)
            continue;

        for(int j = i + 1; j < 3; ++j) {
            const int ch2 = selectedChannels[j];

            if(!block.addbackValid[ch2])
                continue;

            const double e2 = block.addbackEnergy[ch2];
            if(e2 <= 1.0 || e2 >= 10000.0)
                continue;

            d.hMatrix->Fill(e1, e2);
        }
    }
}

RootOutputManager::RootOutputManager() = default;
RootOutputManager::~RootOutputManager() { finalize(); }

bool RootOutputManager::initialize(int nthreads, PixieType type,
                                   bool writeAnalysis, bool writeTree,
                                   bool saveTrace)
{
    finalize();
    impl = std::make_unique<Impl>();
    impl->writeAnalysis = writeAnalysis;
    impl->writeTree = writeTree;
    impl->saveTrace = writeTree && saveTrace;
    impl->type = type;

    // No ROOT output requested.
    if (!writeAnalysis && !writeTree)
        return true;
    if (nthreads <= 0) {
        std::cerr << "RootOutputManager: invalid thread count.\n";
        return false;
    }

    ROOT::EnableThreadSafety();
    TH1::AddDirectory(false);
    impl->threads.resize(static_cast<size_t>(nthreads));

    for (int tid = 0; tid < nthreads; ++tid) {
        auto& d = impl->threads[static_cast<size_t>(tid)];
        char filename[128];
        std::snprintf(filename, sizeof(filename), "thread_%d.root", tid);
        d.file = TFile::Open(filename, "RECREATE");
        if (!d.file || d.file->IsZombie()) {
            std::cerr << "Cannot create ROOT file " << filename << "\n";
            finalize();
            return false;
        }

        d.file->cd();

        // ROOT Tree: event tree and optional trace tree.
        if (writeTree) {
            d.tree = new TTree("tree", "Tree");
            d.tree->Branch("run", &d.run, "run/I");
            d.tree->Branch("coincID", &d.coincID, "coincID/l");
            d.tree->Branch("channel", &d.channel, "channel/I");
            d.tree->Branch("energy", &d.energy, "energy/D");
            d.tree->Branch("coarseT", &d.coarseT, "coarseT/l");
            d.tree->Branch("fineT", &d.fineT, "fineT/D");
            d.tree->Branch("finish_code", &d.finish_code, "finish_code/I");
            if (type != PixieType::PIXIE_500MHZ)
                d.tree->Branch("cfd_forced", &d.cfd_source, "cfd_forced/I");
            else
                d.tree->Branch("cfd_source", &d.cfd_source, "cfd_source/I");

            if (impl->saveTrace) {
                d.traceTree = new TTree("trace_tree", "Trace Tree");
                d.traceTree->Branch("channel", &d.trace_channel, "channel/I");
                d.traceTree->Branch("run", &d.trace_run, "run/I");
                d.traceTree->Branch("coincID", &d.trace_coincID, "coincID/l");
                d.traceTree->Branch("coarseT", &d.trace_coarseT, "coarseT/l");
                d.traceTree->Branch("fineT", &d.trace_fineT, "fineT/D");
                d.traceTree->Branch("trace_length", &d.trace_length, "trace_length/I");
                d.traceTree->Branch("trace", &d.trace_data);
                d.traceTree->Branch("energy", &d.trace_energy, "energy/D");
            }
        }

        // ROOT Analysis: gates, histograms and matrices.
        if (writeAnalysis) {
            impl->createGates(d, tid);
            impl->createHistograms(d);
        }
    }
    return true;
}

void RootOutputManager::fillCoincidence(
    int tid, int run, const CoincidenceBlock& block)
{
    if(!impl || (!impl->writeAnalysis && !impl->writeTree) || tid < 0 ||
       tid >= static_cast<int>(impl->threads.size()))
        return;

    auto& d = impl->threads[static_cast<size_t>(tid)];

    // ROOT Tree: write accepted event rows only when that option is enabled.
    if (impl->writeTree) {
        for(const auto& member : block.members) {
            if(!member.energyValid)
                continue;
            const PixieEvent& ev = member.event;
            d.run = run;
            d.coincID = block.coincID;
            d.channel = member.channel;
            d.finish_code = ev.finish_code;
            d.coarseT = ev.coarse_time;
            d.fineT = ev.fine_time;
            d.cfd_source = ev.cfd_source;
            d.energy = member.calibratedEnergy;
            if(d.tree)
                d.tree->Fill();
        }
    }

    // if (impl->writeTree) {

    //     const int channelPairs[6][2] = {
    //         {80, 81},
    //         {82, 83},
    //         {84, 85},
    //         {86, 87},
    //         {88, 89},
    //         {90, 91}
    //     };

    //     bool validPair = false;

    //     bool channelPresent[92] = {false};

    //     for (const auto& member : block.members) {
    //         if (member.energyValid && member.channel >= 0 && member.channel < 92)
    //             channelPresent[member.channel] = true;
    //     }

    //     for (const auto& pair : channelPairs) {
    //         if (channelPresent[pair[0]] && channelPresent[pair[1]]) {
    //             validPair = true;
    //             break;
    //         }
    //     }

    //     // Fill tree only if at least one complete pair exists
    //     if (validPair) {
    //         for (const auto& member : block.members) {
    //             if (!member.energyValid)
    //                 continue;

    //             const PixieEvent& ev = member.event;

    //             d.run = run;
    //             d.coincID = block.coincID;
    //             d.channel = member.channel;
    //             d.finish_code = ev.finish_code;
    //             d.coarseT = ev.coarse_time;
    //             d.fineT = ev.fine_time;
    //             d.cfd_source = ev.cfd_source;
    //             d.energy = member.calibratedEnergy;

    //             if (d.tree)
    //                 d.tree->Fill();
    //         }
    //     }
    // }

    // ROOT Analysis: fill histograms/matrices only when requested.
    if (impl->writeAnalysis)
        impl->fillRootHistograms(d, block);
}

void RootOutputManager::fillTrace(
    int tid, int run, unsigned long long coincID, int channel,
    int traceLength, double energy, unsigned long long coarseT,
    double fineT, const std::vector<uint16_t>& trace)
{
    if (!impl || !impl->writeTree || !impl->saveTrace || tid < 0 ||
        tid >= static_cast<int>(impl->threads.size()))
        return;

    auto& d = impl->threads[static_cast<size_t>(tid)];
    if (!d.traceTree || trace.empty())
        return;

    d.trace_run = run;
    d.trace_coincID = coincID;
    d.trace_channel = channel;
    d.trace_length = traceLength;
    d.trace_energy = energy;
    d.trace_coarseT = coarseT;
    d.trace_fineT = fineT;
    d.trace_data = trace;
    d.traceTree->Fill();
}

void RootOutputManager::finalize()
{
    if (!impl)
        return;

    if (impl->writeAnalysis || impl->writeTree) {
        for (auto& d : impl->threads) {
            if (d.file) {
                d.file->cd();
                d.file->Write();
                d.file->Close();
                delete d.file;
                d.file = nullptr;
                d.tree = nullptr;
                d.traceTree = nullptr;
                d.hMatrix = nullptr;
                d.hAlphaGamma = nullptr;
                d.hProtonGamma = nullptr;
                for (auto& h : d.hE1vsE2) h = nullptr;
                for (auto& h : d.hAddback) h = nullptr;
            }
        }
    }
    impl->threads.clear();
    impl.reset();
}
