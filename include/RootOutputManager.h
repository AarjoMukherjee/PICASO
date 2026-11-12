#pragma once

#include "Common.h"
#include "CoincidenceProcessor.h"
#include <memory>
#include <vector>
#include <cstdint>

// ROOT output manager. Receives complete coincidence blocks, including the
// precomputed Clover addback energies from CoincidenceProcessor.
class RootOutputManager {
public:
    RootOutputManager();
    ~RootOutputManager();

    RootOutputManager(const RootOutputManager&) = delete;
    RootOutputManager& operator=(const RootOutputManager&) = delete;

    // ROOT Analysis controls histograms/matrices/gates; ROOT Tree controls
    // the event tree and optional trace tree. They can be enabled independently.
    bool initialize(int nthreads, PixieType type,
                    bool writeAnalysis, bool writeTree, bool saveTrace);

    // Writes all event rows for one coincidence and fills ROOT histograms using
    // the addback energies already calculated by CoincidenceProcessor.
    void fillCoincidence(int tid, int run, const CoincidenceBlock& block);

    void fillTrace(int tid, int run, unsigned long long coincID, int channel,
                   int traceLength, double energy,
                   unsigned long long coarseT, double fineT,
                   const std::vector<uint16_t>& trace);

    void finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
