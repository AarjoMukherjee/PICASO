#pragma once
#include "Common.h"
#include <functional>
#include <array>

struct CoincidenceMember {
    int channel = -1;       // Global channel index: channel + (slot - 2) * 16
    PixieEvent event;
    double calibratedEnergy = 0.0;
    bool energyValid = false; // Mirrors the legacy energy acceptance range.
};

constexpr int MARCOS_CLOVER_COUNT = 18;
constexpr int MARCOS_ADDBACK_BINS = 10000;

struct CoincidenceBlock {
    unsigned long long coincID = 0;
    bool finalBlock = false; // True for the last block in a run; preserves legacy tdiff threshold.
    std::vector<CoincidenceMember> members;

    // Calculated exactly once by CoincidenceProcessor for this coincidence.
    // The same values are consumed by BinaryOutputManager and RootOutputManager.
    std::array<double, MARCOS_CLOVER_COUNT> addbackEnergy{};
    std::array<bool, MARCOS_CLOVER_COUNT> addbackValid{};
};

using CoincidenceCallback = std::function<void(const CoincidenceBlock&)>;

// Reads and chronologically merges module files, builds coincidence blocks, calibrates
// accepted event energies, and calculates Clover addback energies. Histogram updates and
// output-file writing belong to output managers.
void process_file(const std::vector<std::string>& filenames,
                  double timewin,
                  PixieType type,
                  int run_number,
                  bool timeSorting,
                  unsigned int *rng_state,
                  const CoincidenceCallback& onCoincidence);
