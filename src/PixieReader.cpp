#include "PixieReader.h"
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <queue>
#include <string>
#include <unistd.h>
#include <utility>

bool valid_header(uint32_t w0)
{
    return ((w0 & 0x0FFFFFF0) == 0x00084020);
}

// ===================== EVENT READ =====================

bool read_event(FILE *f, PixieType type, PixieEvent *ev)
{
    uint32_t header[4];
    if(fread(header, sizeof(uint32_t), 4, f) != 4) return false;

    uint32_t w0 = header[0], w1 = header[1], w2 = header[2], w3 = header[3];
    ev->channel       =  w0 & 0xF;
    ev->slot          = (w0 >> 4)  & 0xF;
    ev->crate         = (w0 >> 8)  & 0xF;
    ev->header_length = (w0 >> 12) & 0x1F;
    ev->event_length  = (w0 >> 17) & 0x3FFF;
    ev->finish_code   = (w0 >> 31) & 0x1;

    uint64_t ts_lo = w1;
    uint64_t ts_hi = (uint64_t)(w2 & 0xFFFF);

    switch(type) {
        case PixieType::PIXIE_100MHZ: {
            uint32_t cfd_frac = (w2 >> 16) & 0x7FFF;
            ev->cfd_fraction = cfd_frac;
            ev->cfd_source = 0;
            ev->cfd_forced = false;
            ev->coarse_time = ((ts_hi << 32) | ts_lo) * 10;
            ev->fine_time = cfd_frac * (10.0 / 32768.0);
            break;
        }
        case PixieType::PIXIE_250MHZ: {
            uint32_t forced = (w2 >> 31) & 0x1;
            uint32_t source = (w2 >> 30) & 0x1;
            uint32_t cfd_frac = (w2 >> 16) & 0x3FFF;
            ev->coarse_time = ((ts_hi << 32) | ts_lo) * 8;
            ev->cfd_fraction = cfd_frac;
            ev->cfd_source = source;
            ev->cfd_forced = forced;
            ev->fine_time = forced ? 0.0 : (-source + cfd_frac / 16384.0) * 4.0;
            break;
        }
        case PixieType::PIXIE_500MHZ: {
            uint32_t cfd_frac = (w2 >> 16) & 0x1FFF;
            uint32_t source = (w2 >> 29) & 0x7;
            ev->coarse_time = ((ts_hi << 32) | ts_lo) * 10;
            ev->cfd_fraction = cfd_frac;
            ev->cfd_source = source;
            ev->cfd_forced = false;
            ev->fine_time = source == 7 ? 0.0 : 2.0 * source + cfd_frac * (2.0 / 8192.0);
            break;
        }
    }
    ev->cfd_time = (double)ev->coarse_time + ev->fine_time;
    ev->energy = w3 & 0xFFFF;
    ev->trace_length = (w3 >> 16) & 0x7FFF;
    ev->trace_out_of_range = (w3 >> 31) & 0x1;

    int extra_header_words = ev->header_length - 4;
    if(extra_header_words > 0 && fseek(f, (long)extra_header_words * 4, SEEK_CUR) != 0)
        return false;

    int event_data_words = ev->event_length - ev->header_length;
    if(event_data_words < 0) return false;

    ev->trace.clear();
    const int trace_words = (ev->trace_length + 1) / 2;
    if(storeTrace && ev->trace_length > 0) {
        ev->trace.resize(ev->trace_length);
        for(int i = 0; i < trace_words; ++i) {
            uint32_t word;
            if(fread(&word, sizeof(uint32_t), 1, f) != 1) return false;
            const int sample = 2 * i;
            if(sample < ev->trace_length) ev->trace[sample] = (uint16_t)(word & 0xFFFF);
            if(sample + 1 < ev->trace_length) ev->trace[sample + 1] = (uint16_t)((word >> 16) & 0xFFFF);
        }
    } else if(trace_words > 0) {
        if(fseek(f, (long)trace_words * 4, SEEK_CUR) != 0) return false;
    }

    // Preserve any event words beyond the trace, if present.
    const int remaining = event_data_words - trace_words;
    if(remaining > 0 && fseek(f, (long)remaining * 4, SEEK_CUR) != 0) return false;
    return true;
}

namespace {
struct RawRecord {
    long double time = 0;
    uint64_t sequence = 0;
    std::vector<uint32_t> words;
};

long double timestamp_from_header(const uint32_t* h, PixieType type)
{
    const uint32_t w1 = h[1], w2 = h[2];
    const uint64_t ticks = (((uint64_t)(w2 & 0xFFFF)) << 32) | w1;
    uint64_t coarse = 0;
    double fine = 0.0;
    switch(type) {
        case PixieType::PIXIE_100MHZ: {
            const uint32_t frac = (w2 >> 16) & 0x7FFF;
            coarse = ticks * 10ULL;
            fine = frac * (10.0 / 32768.0);
            break;
        }
        case PixieType::PIXIE_250MHZ: {
            const uint32_t forced = (w2 >> 31) & 1;
            const uint32_t source = (w2 >> 30) & 1;
            const uint32_t frac = (w2 >> 16) & 0x3FFF;
            coarse = ticks * 8ULL;
            fine = forced ? 0.0 : (-static_cast<double>(source) + frac / 16384.0) * 4.0;
            break;
        }
        case PixieType::PIXIE_500MHZ: {
            const uint32_t frac = (w2 >> 16) & 0x1FFF;
            const uint32_t source = (w2 >> 29) & 0x7;
            coarse = ticks * 10ULL;
            fine = source == 7 ? 0.0 : 2.0 * source + frac * (2.0 / 8192.0);
            break;
        }
    }
    return static_cast<long double>(coarse) + static_cast<long double>(fine);
}

// Read one complete raw record. On the source file, resynchronise by one 32-bit
// word when a plausible event-length/header combination is not found.
bool read_raw_record(FILE* fp, PixieType type, RawRecord& rec, bool recover)
{
    uint32_t h[4];
    while(true) {
        const size_t n = fread(h, sizeof(uint32_t), 4, fp);
        if(n == 0) return false;
        if(n != 4) return false;

        const unsigned headerLength = (h[0] >> 12) & 0x1F;
        const unsigned eventLength = (h[0] >> 17) & 0x3FFF;
        const bool plausible = headerLength >= 4 && headerLength <= eventLength &&
                               eventLength >= 4 && eventLength <= 16383;
        if(!plausible) {
            if(!recover || fseek(fp, -3L * (long)sizeof(uint32_t), SEEK_CUR) != 0)
                return false;
            continue;
        }

        rec.words.resize(eventLength);
        std::copy(h, h + 4, rec.words.begin());
        if(eventLength > 4 && fread(rec.words.data() + 4, sizeof(uint32_t), eventLength - 4, fp) != eventLength - 4)
            return false;
        rec.time = timestamp_from_header(h, type);
        return true;
    }
}

bool write_record(FILE* fp, const RawRecord& rec)
{
    return fwrite(rec.words.data(), sizeof(uint32_t), rec.words.size(), fp) == rec.words.size();
}

bool write_sorted_chunk(const std::string& path, std::vector<RawRecord>& records)
{
    std::stable_sort(records.begin(), records.end(), [](const RawRecord& a, const RawRecord& b) {
        if(a.time != b.time) return a.time < b.time;
        return a.sequence < b.sequence;
    });
    FILE* out = fopen(path.c_str(), "wb");
    if(!out) return false;
    bool ok = true;
    for(const auto& r : records) if(!write_record(out, r)) { ok = false; break; }
    if(fclose(out) != 0) ok = false;
    return ok;
}
}

bool sort_module_file(const std::string& inputPath,
                      const std::string& outputPath,
                      PixieType type,
                      std::size_t chunkLimitBytes)
{
    namespace fs = std::filesystem;
    const auto tempDir = fs::temp_directory_path();
    const std::string base = (tempDir / ("marcos_sort_" + std::to_string((long long)getpid()) +
                                  "_" + std::to_string(std::hash<std::string>{}(inputPath)))).string();
    std::vector<std::string> chunks;
    FILE* in = fopen(inputPath.c_str(), "rb");
    if(!in) {
        std::cerr << "Cannot open module file for sorting: " << inputPath << "\n";
        return false;
    }

    bool ok = true;
    uint64_t seq = 0;
    std::size_t chunkBytes = 0;
    std::vector<RawRecord> records;
    records.reserve(65536);
    while(true) {
        RawRecord rec;
        if(!read_raw_record(in, type, rec, true)) break;
        rec.sequence = seq++;
        chunkBytes += rec.words.size() * sizeof(uint32_t);
        records.push_back(std::move(rec));
        if(chunkBytes >= chunkLimitBytes) {
            const std::string path = base + "_chunk_" + std::to_string(chunks.size()) + ".bin";
            if(!write_sorted_chunk(path, records)) { ok = false; break; }
            chunks.push_back(path);
            records.clear();
            chunkBytes = 0;
        }
    }
    if(ferror(in)) ok = false;
    fclose(in);

    if(ok && !records.empty()) {
        const std::string path = base + "_chunk_" + std::to_string(chunks.size()) + ".bin";
        if(!write_sorted_chunk(path, records)) ok = false;
        else chunks.push_back(path);
    }

    if(ok) {
        FILE* out = fopen(outputPath.c_str(), "wb");
        if(!out) ok = false;
        else {
            struct Node { RawRecord rec; std::size_t index; };
            struct Later {
                bool operator()(const Node& a, const Node& b) const {
                    if(a.rec.time != b.rec.time) return a.rec.time > b.rec.time;
                    if(a.rec.sequence != b.rec.sequence) return a.rec.sequence > b.rec.sequence;
                    return a.index > b.index;
                }
            };
            std::vector<FILE*> inputs(chunks.size(), nullptr);
            std::priority_queue<Node, std::vector<Node>, Later> queue;
            for(std::size_t i = 0; i < chunks.size(); ++i) {
                inputs[i] = fopen(chunks[i].c_str(), "rb");
                if(!inputs[i]) { ok = false; break; }
                RawRecord r;
                if(read_raw_record(inputs[i], type, r, false)) {
                    // sequence within chunks is not persisted; ordering of equal timestamps
                    // is immaterial to the chronological merge.
                    queue.push(Node{std::move(r), i});
                }
            }
            while(ok && !queue.empty()) {
                Node node = queue.top(); queue.pop();
                if(!write_record(out, node.rec)) { ok = false; break; }
                RawRecord next;
                if(read_raw_record(inputs[node.index], type, next, false))
                    queue.push(Node{std::move(next), node.index});
            }
            for(FILE* fp : inputs) if(fp) fclose(fp);
            if(fclose(out) != 0) ok = false;
        }
    }

    for(const auto& path : chunks) {
        std::error_code ec;
        fs::remove(path, ec);
    }
    if(!ok) {
        std::error_code ec;
        fs::remove(outputPath, ec);
        std::cerr << "Failed sorting module file: " << inputPath << "\n";
    }
    return ok;
}
