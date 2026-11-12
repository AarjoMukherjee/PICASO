#pragma once
#include <cstdint>
#include <cstdio>
#include <vector>
#include <string>

#define OneDHisSize 4194304
#define TwoDHisDim 4096
#define OneDHisDim 8192
#define MAX_FILES 500

struct Histograms {
    int *sp = nullptr;
    int *tsp = nullptr;
    int *Twod = nullptr;
};

enum class PixieType { PIXIE_100MHZ, PIXIE_250MHZ, PIXIE_500MHZ };

struct PixieEvent {
    int channel = 0, slot = 0, crate = 0;
    int header_length = 0, event_length = 0, finish_code = 0;
    unsigned long long coarse_time = 0;
    double fine_time = 0, cfd_time = 0;
    int cfd_fraction = 0, cfd_source = 0;
    bool cfd_forced = false;
    int energy = 0;
    int trace_length = 0;
    bool trace_out_of_range = false;
    std::vector<uint16_t> trace;
    unsigned int qdc[8] = {};
    unsigned long long external_ts = 0;
};

extern float Offset[96], Slope[96], Quad[96], Cube[96], SqrtCoeff[96];
extern float TimeBinWidthPs, TimeScale, TimeOffsetBin, EnergyScale;
extern int ECal, rootWrite, storeTrace;
