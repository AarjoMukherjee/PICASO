#include "Common.h"
#include "BinaryOutputManager.h"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <new>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <string>

void allocate_histograms(Histograms *H)
{
    H->sp = (int*)calloc(OneDHisSize/4, sizeof(int));
    H->tsp = (int*)calloc(OneDHisSize/4, sizeof(int));
    H->Twod = (int*)calloc((size_t)TwoDHisDim*TwoDHisDim, sizeof(int));
    if(!H->sp || !H->tsp || !H->Twod) {
        free_histograms(H);
        throw std::bad_alloc();
    }
}

static void merge_arrays(int *g, const int *l, size_t n)
{
    #pragma omp simd
    for(size_t i = 0; i < n; ++i) g[i] += l[i];
}

void merge_histograms(Histograms *G, Histograms *L)
{
    merge_arrays(G->sp, L->sp, OneDHisSize/4);
    merge_arrays(G->tsp, L->tsp, OneDHisSize/4);
    merge_arrays(G->Twod, L->Twod, (size_t)TwoDHisDim*TwoDHisDim);
}

void free_histograms(Histograms *H)
{
    free(H->sp); free(H->tsp); free(H->Twod);
    H->sp = H->tsp = H->Twod = nullptr;
}

void add_coincidence_to_histograms(const CoincidenceBlock& block, Histograms *H)
{
    double rawEnergy[96] = {};
    double cfdTime[96] = {};
    bool present[96] = {};

    for(const auto& member : block.members) {
        const int ch = member.channel;
        if(ch < 0 || ch >= 96) continue;
        rawEnergy[ch] = member.event.energy;
        cfdTime[ch] = member.event.cfd_time;
        present[ch] = true;

        // Calibration is performed by CoincidenceProcessor. This manager only
        // increments the calibrated spectrum histogram.
        if(!member.energyValid) continue;
        const double GEn = member.calibratedEnergy;
        const int base = ch << 13;
        const int idx = base + static_cast<int>(GEn);
        // Keep each channel within its own 8192-bin region. This also protects
        // the reserved addback spectra stored after the 96 channel regions.
        if(idx >= base && idx < base + 8192 && idx < OneDHisSize / 4)
            H->sp[idx]++;
    }

    // Clover addback spectra occupy the unused tail of output_1d.his after the
    // 96 x 8192 individual-channel spectra. Each Clover has 10000 bins [0,10000),
    // matching the ROOT addback histogram range. The 18 spectra fit in the
    // remaining OneDHisSize/4 storage.
    // constexpr size_t channelSpectrumBins = 96u * 8192u;
    // constexpr size_t addbackBase = channelSpectrumBins;
    // for(int ch1 = 0; ch1 < MARCOS_CLOVER_COUNT; ++ch1) {
    //     if(!block.addbackValid[ch1])
    //         continue;
    //     const double energy = block.addbackEnergy[ch1];
    //     if(energy < 0.0 || energy >= MARCOS_ADDBACK_BINS)
    //         continue;
    //     const size_t idx = addbackBase + static_cast<size_t>(ch1) * MARCOS_ADDBACK_BINS
    //                      + static_cast<size_t>(energy);
    //     if(idx < OneDHisSize / 4)
    //         H->sp[idx]++;
        // for(int ch2 = ch1 + 1; ch2 < MARCOS_CLOVER_COUNT; ++ch2) {
        //     if(!block.addbackValid[ch2])
        //         continue;
        //     const double energy2 = block.addbackEnergy[ch2];
        //     if(energy2 < 0.0 || energy2 >= MARCOS_ADDBACK_BINS)
        //         continue;
        //     const size_t matrixIdx = static_cast<size_t>(ch1) * TwoDHisDim
        //                     + static_cast<size_t>(ch2);
        //     if(matrixIdx < static_cast<size_t>(TwoDHisDim) * TwoDHisDim)
        //         H->Twod[matrixIdx]++;
        // }
    // }

    constexpr size_t channelSpectrumBins = 96u * 8192u;
    constexpr size_t addbackBase = channelSpectrumBins;

    const int selectedChannels[] = {9, 11, 12};

    for(int ch1 = 0; ch1 < MARCOS_CLOVER_COUNT; ++ch1) {

        if(!block.addbackValid[ch1])
            continue;

        const double energy1 = block.addbackEnergy[ch1];

        if(energy1 < 0.0 || energy1 >= MARCOS_ADDBACK_BINS)
            continue;

        // 1D Addback spectrum: ALL channels
        const size_t idx = addbackBase
                        + static_cast<size_t>(ch1) * MARCOS_ADDBACK_BINS
                        + static_cast<size_t>(energy1);

        if(idx < OneDHisSize / 4)
            H->sp[idx]++;

        // 2D Matrix: ONLY channels 9, 11, 12
        if(ch1 != 9 && ch1 != 11 && ch1 != 12)
            continue;

        for(int i = 0; i < 3; ++i) {

            const int ch2 = selectedChannels[i];

            if(ch2 <= ch1)
                continue;

            if(!block.addbackValid[ch2])
                continue;

            const double energy2 = block.addbackEnergy[ch2];

            if(energy2 < 0.0 || energy2 >= MARCOS_ADDBACK_BINS)
                continue;

            const size_t matrixIdx = static_cast<size_t>(energy1) * TwoDHisDim
                                + static_cast<size_t>(energy2);

            if(matrixIdx < static_cast<size_t>(TwoDHisDim) * TwoDHisDim)
                H->Twod[matrixIdx]++;
        }
    }

    // Time-difference histogram increments remain exclusively in this manager.
    const double emin = 10;
    const double emax = 100000.0;
    if(present[0] && rawEnergy[0] >= emin && rawEnergy[0] <= emax) {
        for(int ch = 1; ch < 96; ++ch) {
            if(!present[ch] || rawEnergy[ch] < emin || rawEnergy[ch] > emax) continue;
            const double dt = TimeScale * (cfdTime[0] - cfdTime[ch]);
            const int idx = OneDHisDim * ch + TimeOffsetBin + static_cast<int>(dt);
            if(idx >= 0 && idx < OneDHisSize / 4) H->tsp[idx]++;
        }
    }
}

static bool write_one(const char *name, const int *data, size_t n)
{
    FILE *fp = fopen(name, "wb");
    if(!fp) { std::cerr << "Cannot create " << name << "\n"; return false; }
    const size_t written = fwrite(data, sizeof(int), n, fp);
    fclose(fp);
    if(written != n) { std::cerr << "Incomplete write to " << name << "\n"; return false; }
    return true;
}

namespace {
constexpr int SPK_DIRECTORY_BYTES = 2048;
constexpr int SPK_HEADER_BYTES = 128;
constexpr int SPK_HEADER_HALFWORDS = 64;
constexpr int SPK_INITIAL_NXWD = 1025;
constexpr int SPK_MAX_SPECTRA = 254;
constexpr int SPK_CHANNEL_SPECTRA = 96;
constexpr int SPK_CHANNEL_BINS = 8192;
constexpr int SPK_ADDBACK_BINS = MARCOS_ADDBACK_BINS;
constexpr int SPK_TOTAL_SPECTRA = SPK_CHANNEL_SPECTRA + MARCOS_CLOVER_COUNT;
constexpr size_t SPK_ADDBACK_BASE =
    static_cast<size_t>(SPK_CHANNEL_SPECTRA) * SPK_CHANNEL_BINS;

void spk_set_int32_le(uint8_t *buffer, size_t offset, int32_t value)
{
    const uint32_t v = static_cast<uint32_t>(value);
    buffer[offset]     = static_cast<uint8_t>(v & 0xffu);
    buffer[offset + 1] = static_cast<uint8_t>((v >> 8) & 0xffu);
    buffer[offset + 2] = static_cast<uint8_t>((v >> 16) & 0xffu);
    buffer[offset + 3] = static_cast<uint8_t>((v >> 24) & 0xffu);
}

bool spk_write_int32_le(std::ofstream& out, int32_t value)
{
    uint8_t bytes[4];
    spk_set_int32_le(bytes, 0, value);
    out.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
    return static_cast<bool>(out);
}

void spk_create_header(uint8_t *header, int spectrumID,
                       int logicalChannels, int minChannel, int maxChannel)
{
    std::memset(header, 0, SPK_HEADER_BYTES);
    spk_set_int32_le(header, 0, spectrumID);

    char label[13];
    std::snprintf(label, sizeof(label), "%-11d*", spectrumID);
    std::memcpy(header + 4, label, 12);

    // DAMM/SPK header fields. The logical histogram length is retained here;
    // only the non-zero MINC..MAXC interval is physically written below.
    spk_set_int32_le(header, 7 * 4, -4);                    // bytes/channel
    spk_set_int32_le(header, 8 * 4, SPK_HEADER_HALFWORDS); // header length
    spk_set_int32_le(header, 9 * 4, 2 * logicalChannels);
    spk_set_int32_le(header, 10 * 4, 1);                   // dimensionality
    spk_set_int32_le(header, 11 * 4, logicalChannels);
    spk_set_int32_le(header, 12 * 4, logicalChannels);
    spk_set_int32_le(header, 13 * 4, logicalChannels);
    spk_set_int32_le(header, 14 * 4, 0);
    spk_set_int32_le(header, 15 * 4, minChannel);
    spk_set_int32_le(header, 16 * 4, maxChannel);

    std::memset(header + 22 * 4, ' ', 40);
    std::memcpy(header + 22 * 4 + 36, "****", 4);
}

bool spk_write_spectrum(std::ofstream& out, const std::vector<int32_t>& his,
                        size_t base, int logicalChannels, int spectrumID,
                        int32_t& nextHalfWord, uint8_t *directory)
{
    int minChannel = logicalChannels;
    int maxChannel = -1;
    for(int ch = 0; ch < logicalChannels; ++ch) {
        if(his[base + static_cast<size_t>(ch)] != 0) {
            minChannel = std::min(minChannel, ch);
            maxChannel = std::max(maxChannel, ch);
        }
    }

    // SPKIO convention for an empty spectrum: store one zero channel.
    if(maxChannel < 0) minChannel = maxChannel = 0;
    const int storedChannels = maxChannel - minChannel + 1;

    // Four directory words precede the two-word entries for spectra.
    const size_t directoryOffset = 4u * (4u + 2u * (spectrumID - 1u));
    spk_set_int32_le(directory, directoryOffset, spectrumID);
    spk_set_int32_le(directory, directoryOffset + 4, nextHalfWord);

    uint8_t header[SPK_HEADER_BYTES];
    spk_create_header(header, spectrumID, logicalChannels,
                      minChannel, maxChannel);
    out.write(reinterpret_cast<const char*>(header), sizeof(header));
    if(!out) return false;

    for(int ch = minChannel; ch <= maxChannel; ++ch) {
        if(!spk_write_int32_le(out, his[base + static_cast<size_t>(ch)]))
            return false;
    }

    nextHalfWord += SPK_HEADER_HALFWORDS + 2 * storedChannels;
    return true;
}
} // namespace

bool write_spk_file(const char *hisFilename, const char *spkFilename)
{
    if(!hisFilename || !spkFilename) {
        std::cerr << "SPK conversion: invalid filename.\n";
        return false;
    }
    if(SPK_TOTAL_SPECTRA > SPK_MAX_SPECTRA) {
        std::cerr << "SPK conversion: too many spectra for directory.\n";
        return false;
    }

    std::ifstream in(hisFilename, std::ios::binary | std::ios::ate);
    if(!in) {
        std::cerr << "SPK conversion: cannot open HIS file " << hisFilename << "\n";
        return false;
    }
    const std::streamoff fileSize = in.tellg();
    const size_t expectedBytes = static_cast<size_t>(OneDHisSize);
    if(fileSize < 0 || static_cast<size_t>(fileSize) != expectedBytes) {
        std::cerr << "SPK conversion: unexpected HIS size. Expected "
                  << expectedBytes << " bytes, got " << fileSize << ".\n";
        return false;
    }
    in.seekg(0, std::ios::beg);

    std::vector<uint8_t> raw(expectedBytes);
    in.read(reinterpret_cast<char*>(raw.data()),
            static_cast<std::streamsize>(raw.size()));
    if(!in) {
        std::cerr << "SPK conversion: failed reading HIS file.\n";
        return false;
    }

    // output_1d.his is written as native 32-bit ints. Decode explicitly as
    // little-endian to match the SPK format and the existing HIS writer.
    std::vector<int32_t> his(OneDHisSize / 4);
    for(size_t i = 0; i < his.size(); ++i) {
        const size_t j = 4 * i;
        const uint32_t value =
              static_cast<uint32_t>(raw[j])
            | (static_cast<uint32_t>(raw[j + 1]) << 8)
            | (static_cast<uint32_t>(raw[j + 2]) << 16)
            | (static_cast<uint32_t>(raw[j + 3]) << 24);
        his[i] = static_cast<int32_t>(value);
    }

    std::ofstream out(spkFilename, std::ios::binary | std::ios::trunc);
    if(!out) {
        std::cerr << "SPK conversion: cannot create " << spkFilename << "\n";
        return false;
    }

    uint8_t directory[SPK_DIRECTORY_BYTES];
    std::memset(directory, 0, sizeof(directory));
    std::memcpy(directory, "HHIRFSPK", 8);
    spk_set_int32_le(directory, 8, 0); // number of spectra, finalized later
    spk_set_int32_le(directory, 12, SPK_INITIAL_NXWD);
    out.write(reinterpret_cast<const char*>(directory), sizeof(directory));
    if(!out) return false;

    int32_t nextHalfWord = SPK_INITIAL_NXWD;
    int spectrumID = 1;

    // Spectra 1..96: one 8192-channel spectrum for each global MARCOS channel.
    for(int ch = 0; ch < SPK_CHANNEL_SPECTRA; ++ch, ++spectrumID) {
        const size_t base = static_cast<size_t>(ch) * SPK_CHANNEL_BINS;
        if(!spk_write_spectrum(out, his, base, SPK_CHANNEL_BINS,
                               spectrumID, nextHalfWord, directory)) {
            std::cerr << "SPK conversion: failed writing channel spectrum "
                      << ch << ".\n";
            return false;
        }
    }

    // Spectra 97..114: Clover addback spectra, each with 10000 channels.
    // These are stored contiguously after the 96 channel spectra in output_1d.his.
    for(int c = 0; c < MARCOS_CLOVER_COUNT; ++c, ++spectrumID) {
        const size_t base = SPK_ADDBACK_BASE
                          + static_cast<size_t>(c) * SPK_ADDBACK_BINS;
        if(base + SPK_ADDBACK_BINS > his.size()) {
            std::cerr << "SPK conversion: addback spectrum exceeds HIS buffer.\n";
            return false;
        }
        if(!spk_write_spectrum(out, his, base, SPK_ADDBACK_BINS,
                               spectrumID, nextHalfWord, directory)) {
            std::cerr << "SPK conversion: failed writing Clover addback spectrum "
                      << c << ".\n";
            return false;
        }
    }

    spk_set_int32_le(directory, 8, SPK_TOTAL_SPECTRA);
    spk_set_int32_le(directory, 12, nextHalfWord);
    out.seekp(0, std::ios::beg);
    out.write(reinterpret_cast<const char*>(directory), sizeof(directory));
    out.close();
    if(!out) {
        std::cerr << "SPK conversion: error finalizing output file.\n";
        return false;
    }

    std::cout << "SPK file created: " << spkFilename << " ("
              << SPK_TOTAL_SPECTRA << " spectra: 96 channel + 18 addback)\n";
    return true;
}

bool write_binary_histograms(const Histograms& H)
{
    const bool a = write_one("output_1d.his", H.sp, OneDHisSize/4);
    const bool b = write_one("output_tdiff.his", H.tsp, OneDHisSize/4);

    if(!a || !b)
        return false;

    // Convert completed 1D HIS files to DAMM SPK.
    const bool c = write_spk_file("output_1d.his", "output_1d.spk");
    const bool e = write_spk_file("output_tdiff.his", "output_tdiff.spk");

    // Write 2D coincidence matrix.
    const bool d = write_one(
        "output_2d.spn",
        H.Twod,
        (size_t)TwoDHisDim * TwoDHisDim
    );

    return c && d && e;
}