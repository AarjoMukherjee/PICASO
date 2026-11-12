#include "CoincidenceProcessor.h"
#include "PixieReader.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <filesystem>
#include <limits>
#include <utility>
#include <unistd.h>

void process_file(const std::vector<std::string>& filenames,
                  double timewin,
                  PixieType type,
                  int run_number,
                  bool timeSorting,
                  unsigned int *rng_state,
                  const CoincidenceCallback& onCoincidence)
{
    // If TimeSorting is enabled, sort each module independently on disk.
    // If disabled, read the original module files directly without creating
    // sorted copies. In either mode, the module streams are combined using
    // the same one-look-ahead k-way merge.
    namespace fs = std::filesystem;
    std::vector<std::string> sortedFilenames;
    std::vector<std::string> processingFilenames;
    processingFilenames.reserve(filenames.size());

    const std::string tempBase = (fs::temp_directory_path() /
        ("marcos_run_" + std::to_string(run_number) + "_pid_" +
         std::to_string((long long)getpid()))).string();

    if(timeSorting) {
        sortedFilenames.reserve(filenames.size());
        for(size_t i = 0; i < filenames.size(); ++i) {
            const std::string sortedPath = tempBase + "_module_" + std::to_string(i) + ".sorted.bin";
            if(!sort_module_file(filenames[i], sortedPath, type)) {
                std::cerr << "Aborting run because module sorting failed: " << filenames[i] << "\n";
                for(const auto& path : sortedFilenames) {
                    std::error_code ec;
                    fs::remove(path, ec);
                }
                return;
            }
            sortedFilenames.push_back(sortedPath);
        }
        processingFilenames = sortedFilenames;
        std::cout << "TimeSorting enabled: using sorted temporary module streams for run "
                  << run_number << "\n";
    } else {
        // Use the raw input files directly, preserving their on-disk event order.
        processingFilenames = filenames;
        std::cout << "TimeSorting disabled: reading raw module files as-is for run "
                  << run_number << "\n";
    }

    auto removeSortedFiles = [&]() {
        for(const auto& path : sortedFilenames) {
            std::error_code ec;
            fs::remove(path, ec);
        }
    };

    std::vector<FILE*> files(processingFilenames.size(), nullptr);
    std::vector<PixieEvent> pending(processingFilenames.size());
    std::vector<unsigned char> valid(processingFilenames.size(), 0);

    size_t active = 0;
    for(size_t i = 0; i < processingFilenames.size(); ++i) {
        files[i] = fopen(processingFilenames[i].c_str(), "rb");
        if(!files[i]) {
            std::cerr << "Cannot open module file: " << processingFilenames[i] << "\n";
            continue;
        }
        if(read_event(files[i], type, &pending[i])) {
            valid[i] = 1;
            ++active;
        } else {
            fclose(files[i]);
            files[i] = nullptr;
        }
    }

    if(active == 0) {
        for(FILE *fp : files) if(fp) fclose(fp);
        removeSortedFiles();
        std::cerr << "No readable events for run " << run_number << "\n";
        return;
    }

    double CTime0 = -1.0;
    double FTime0 = -1.0;
    bool first_event = true;
    unsigned long long coincID = 0;
    bool present[96] = {};
    PixieEvent blockEvents[96];

    // Clover membership is defined centrally here. Output managers must not
    // independently reconstruct addback from the individual detector energies.
    static constexpr int cloverChannels[MARCOS_CLOVER_COUNT][4] = {
        {100,100,100,100}, {16,17,18,19}, {0,1,2,3}, {44,45,46,47},
        {4,5,6,7}, {20,21,22,23}, {24,25,26,27}, {28,29,30,31},
        {32,33,34,35}, {36,37,38,39}, {100,100,100,100}, {48,49,50,51},
        {68,69,70,71}, {72,73,74,75}, {60,61,62,63}, {76,77,78,79},
        {64,65,66,67}, {56,57,58,59}
    };

    auto emitBlock = [&](bool finalBlock) {
        CoincidenceBlock block;
        block.coincID = coincID;
        block.finalBlock = finalBlock;
        for(int i = 0; i < 96; ++i) {
            if(!present[i]) continue;

            CoincidenceMember member;
            member.channel = i;
            member.event = blockEvents[i];

            // Keep the legacy accepted raw-energy range and calibration formula.
            // The random dither is applied once per accepted event, in channel order,
            // preserving the previous BinaryOutputManager calibration order.
            if(member.event.energy >= 0 && member.event.energy < 100000) {
                const double E = EnergyScale * member.event.energy;
                member.calibratedEnergy = Offset[i] + Slope[i] * E + Quad[i] * E * E +
                    Cube[i] * E * E * E + SqrtCoeff[i] * std::sqrt(E) +
                    (rng_state ? rand_r(rng_state) * (1.0 / RAND_MAX) : 0.0);
                member.energyValid = true;
            }
            block.members.push_back(std::move(member));
        }

        // Calculate addback from calibrated, accepted member energies. Keep the
        // legacy ROOT selection: each crystal energy must be >1 and <40000 keV.
        double channelEnergy[96] = {};
        bool channelValid[96] = {};
        for(const auto& member : block.members) {
            if(member.channel < 0 || member.channel >= 96 || !member.energyValid)
                continue;
            channelEnergy[member.channel] = member.calibratedEnergy;
            channelValid[member.channel] = true;
        }

        for(int c = 0; c < MARCOS_CLOVER_COUNT; ++c) {
            double sum = 0.0;
            int hits = 0;
            for(int j = 0; j < 4; ++j) {
                const int ch = cloverChannels[c][j];
                if(ch < 0 || ch >= 96 || !channelValid[ch])
                    continue;
                const double energy = channelEnergy[ch];
                if(energy > 1.0 && energy < 40000.0) {
                    sum += energy;
                    ++hits;
                }
            }
            if(hits > 0) {
                block.addbackEnergy[c] = sum;
                block.addbackValid[c] = true;
            }
        }

        if(!block.members.empty() && onCoincidence) onCoincidence(block);
        for(int i = 0; i < 96; ++i) present[i] = false;
    };

    while(active > 0) {
        size_t earliest = files.size();
        long double earliest_time = std::numeric_limits<long double>::infinity();
        for(size_t i = 0; i < files.size(); ++i) {
            if(!valid[i]) continue;
            const long double event_time =
                static_cast<long double>(pending[i].coarse_time) +
                static_cast<long double>(pending[i].fine_time);
            if(event_time < earliest_time) {
                earliest_time = event_time;
                earliest = i;
            }
        }
        if(earliest == files.size()) break;

        PixieEvent ev = std::move(pending[earliest]);
        if(read_event(files[earliest], type, &pending[earliest])) {
            valid[earliest] = 1;
        } else {
            valid[earliest] = 0;
            --active;
            fclose(files[earliest]);
            files[earliest] = nullptr;
        }

        const int ch = ev.channel + (ev.slot - 2) * 16;
        if(ch < 0 || ch >= 96) continue;

        double dt = 0.0;
        if(first_event) {
            CTime0 = static_cast<double>(ev.coarse_time);
            FTime0 = ev.fine_time;
            first_event = false;
        } else {
            dt = static_cast<double>(ev.coarse_time) - CTime0 + ev.fine_time - FTime0;
        }

        // Do not discard a raw event solely because its timestamp moves backward.
        // With TimeSorting disabled, the original stream order is retained.
        if(std::fabs(dt) < timewin) {
            present[ch] = true;
            blockEvents[ch] = std::move(ev); // Keep one event per channel, matching previous behavior.
        } else {
            emitBlock(false);
            ++coincID;
            CTime0 = static_cast<double>(ev.coarse_time);
            FTime0 = ev.fine_time;
            present[ch] = true;
            blockEvents[ch] = std::move(ev);
        }
    }

    emitBlock(true);
    for(FILE *fp : files) if(fp) fclose(fp);
    removeSortedFiles();
}
