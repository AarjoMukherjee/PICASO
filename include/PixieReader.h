#pragma once
#include "Common.h"
bool valid_header(uint32_t w0);
bool read_event(FILE *f, PixieType type, PixieEvent *ev);

// Sort one PIXIE module file by event timestamp using bounded-memory external sorting.
// Returns false on failure; sorted output is written to outputPath.
bool sort_module_file(const std::string& inputPath,
                      const std::string& outputPath,
                      PixieType type,
                      std::size_t chunkLimitBytes = 64u * 1024u * 1024u);
