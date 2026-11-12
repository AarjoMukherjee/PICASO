#pragma once
#include "Common.h"
#include "CoincidenceProcessor.h"

// Non-ROOT histogram storage, histogram additions, merging, and binary .his output.
// Addback energies are supplied precomputed in CoincidenceBlock by CoincidenceProcessor.
void allocate_histograms(Histograms *H);
void merge_histograms(Histograms *G, Histograms *L);
void free_histograms(Histograms *H);
void add_coincidence_to_histograms(const CoincidenceBlock& block, Histograms *H);
bool write_binary_histograms(const Histograms& H);

// Convert the completed output_1d.his file into a DAMM-compatible SPK file.
// The converter writes 96 channel spectra (8192 bins each) followed by
// 18 Clover addback spectra (10000 bins each).
bool write_spk_file(const char *hisFilename = "output_1d.his",
                    const char *spkFilename = "output_1d.spk");
