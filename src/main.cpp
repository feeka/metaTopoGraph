#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _WIN32
#include <direct.h>
#define PATH_SEP '\\'
static bool MakeDir(const std::string& p) { return _mkdir(p.c_str()) == 0; }
static void RemoveDir(const std::string& p) {
    std::string cmd = "rmdir /S /Q \"" + p + "\"";
    system(cmd.c_str());
}
#else
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#define PATH_SEP '/'
static bool MakeDir(const std::string& p) { return mkdir(p.c_str(), 0755) == 0; }
static void RemoveDir(const std::string& p) {
    std::string cmd = "rm -rf \"" + p + "\"";
    system(cmd.c_str());
}
// Returns 90 % of total physical RAM in bytes, matching MEGAHIT's default.
static uint64_t AutoMemBytes() {
    struct sysinfo si;
    if (sysinfo(&si) == 0)
        return static_cast<uint64_t>(si.totalram * si.mem_unit * 0.9);
    return static_cast<uint64_t>(8ULL * 1073741824ULL); // fallback: 8 GB
}
#endif

#include "sdbg/sdbg.h"
#include "topo_extractor.h"
#include "topo_json_writer.h"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static std::string GetBinaryDir(const char* argv0) {
    std::string p(argv0);
    size_t sep = p.find_last_of("/\\");
    return sep == std::string::npos ? std::string(".") : p.substr(0, sep);
}

// Find megahit_core_no_hw_accel or megahit_core.
// Search order: next to this binary, then every directory in PATH.
static std::string FindMegahitCore(const std::string& bin_dir) {
#ifdef _WIN32
    const char* ext = ".exe";
    const char  path_sep = ';';
#else
    const char* ext = "";
    const char  path_sep = ':';
#endif
    const char* names[] = {"megahit_core_no_hw_accel", "megahit_core", nullptr};

    // Build search directories: binary dir first, then PATH entries.
    std::vector<std::string> dirs;
    dirs.push_back(bin_dir);
    const char* env_path = getenv("PATH");
    if (env_path) {
        std::string p(env_path);
        std::string::size_type start = 0, end;
        while ((end = p.find(path_sep, start)) != std::string::npos) {
            dirs.push_back(p.substr(start, end - start));
            start = end + 1;
        }
        dirs.push_back(p.substr(start));
    }

    for (const std::string& dir : dirs) {
        for (int i = 0; names[i]; ++i) {
            std::string path = dir + PATH_SEP + names[i] + ext;
            FILE* f = fopen(path.c_str(), "rb");
            if (f) { fclose(f); return path; }
        }
    }
    return "";
}

// Write a MEGAHIT reads-lib text file.  Format per MEGAHIT source:
//   <metadata line>   <- any string (used for logging)
//   <type> <file1> [<file2>]
// Then megahit_core buildlib converts it to .bin + .lib_info,
// and megahit_core read2sdbg consumes those binary files.
static void WriteLibFile(const std::string& reads1,
                         const std::string& reads2,   // empty for SE
                         const std::string& lib_path) {
    FILE* f = fopen(lib_path.c_str(), "w");
    if (!f)
        throw std::runtime_error("Cannot create lib file: " + lib_path);
    if (reads2.empty()) {
        fprintf(f, "%s\nse %s\n", reads1.c_str(), reads1.c_str());
    } else {
        fprintf(f, "%s,%s\npe %s %s\n",
                reads1.c_str(), reads2.c_str(),
                reads1.c_str(), reads2.c_str());
    }
    fclose(f);
}

// Invoke megahit_core {buildlib, read2sdbg} and return the SDBG prefix.
static std::string BuildSdbgFromReads(const std::string& reads1,
                                      const std::string& reads2,  // empty = SE
                                      const std::string& bin_dir,
                                      const std::string& tmp_dir,
                                      int n_threads,
                                      double host_mem_gb,
                                      int min_count) {
    const std::string core = FindMegahitCore(bin_dir);
    if (core.empty())
        throw std::runtime_error(
            "Cannot find megahit_core_no_hw_accel or megahit_core in: " + bin_dir +
            "\nBuild MEGAHIT first (cmake && make in libs/megahit/build).");

    if (!MakeDir(tmp_dir))
        throw std::runtime_error("Cannot create temp directory: " + tmp_dir);

    // ---- Step 1: write lib text file and convert FASTQ to binary ----
    const std::string lib_prefix = tmp_dir + PATH_SEP + "reads_lib";
    WriteLibFile(reads1, reads2, lib_prefix + ".txt");

    {
        std::ostringstream cmd;
        cmd << "\"" << core << "\" buildlib"
            << " \"" << lib_prefix << ".txt\""
            << " \"" << lib_prefix << "\"";
        std::cerr << "buildlib: " << cmd.str() << "\n";
        int ret = system(cmd.str().c_str());
        if (ret != 0)
            throw std::runtime_error(
                "megahit_core buildlib failed (exit " + std::to_string(ret) + ")");
    }

    // ---- Step 2: build SDBG from the binary lib ----
    const std::string sdbg_prefix = tmp_dir + PATH_SEP + "k21";
    const uint64_t    mem_bytes   = host_mem_gb > 0.0
        ? static_cast<uint64_t>(host_mem_gb * 1073741824.0)
        : AutoMemBytes();

    {
        std::ostringstream cmd;
        cmd << "\"" << core << "\" read2sdbg"
            << " --kmer_k 21"
            << " --min_kmer_frequency " << min_count
            << " --read_lib_file \""    << lib_prefix   << "\""
            << " --output_prefix \""    << sdbg_prefix  << "\""
            << " --num_cpu_threads "    << n_threads
            << " --host_mem "           << mem_bytes
            << " --need_mercy 1";
        std::cerr << "read2sdbg: " << cmd.str() << "\n";
        int ret = system(cmd.str().c_str());
        if (ret != 0)
            throw std::runtime_error(
                "megahit_core read2sdbg failed (exit " + std::to_string(ret) + ")");
    }

    return sdbg_prefix;
}

// ---------------------------------------------------------------------------
// Usage
// ---------------------------------------------------------------------------

static void PrintUsage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  From pre-built SDBG:\n"
        << "    " << prog << " --graph <sdbg_prefix> --output <out.json> [options]\n\n"
        << "  From single-end reads:\n"
        << "    " << prog << " --reads <reads.fa/fq[.gz]> --output <out.json> [options]\n\n"
        << "  From paired-end reads:\n"
        << "    " << prog << " --reads <R1.fq.gz> --reads2 <R2.fq.gz> --output <out.json> [options]\n\n"
        << "Options:\n"
        << "  --graph      SDBG file prefix\n"
        << "  --reads      R1 (or SE) FASTA/FASTQ reads file\n"
        << "  --reads2     R2 FASTA/FASTQ file (paired-end only)\n"
        << "  --output     Output JSON file path\n"
        << "  --sample     Edges to sample for node/walk/bubble phases (default 100000)\n"
        << "  --threads    OpenMP threads (default: all available)\n"
        << "  --mem        Memory for SDBG build in GB, reads mode only (default: 90% of RAM)\n"
        << "  --min-count  Min k-mer frequency for SDBG build (default 2)\n";
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    std::string graph_prefix;
    std::string reads_file;
    std::string reads2_file;
    std::string output_path;
    ExtractionOptions opts;
    int     n_threads  = 0;
    double  mem_gb     = 0.0;  // 0 = auto (90% of system RAM, matching MEGAHIT default)
    int     min_count  = 2;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if      (arg == "--graph"     && i + 1 < argc) graph_prefix = argv[++i];
        else if (arg == "--reads"     && i + 1 < argc) reads_file   = argv[++i];
        else if (arg == "--reads2"    && i + 1 < argc) reads2_file  = argv[++i];
        else if (arg == "--output"    && i + 1 < argc) output_path  = argv[++i];
        else if (arg == "--sample"    && i + 1 < argc) opts.sample_size = static_cast<uint64_t>(std::atoll(argv[++i]));
        else if (arg == "--threads"   && i + 1 < argc) n_threads    = std::atoi(argv[++i]);
        else if (arg == "--mem"       && i + 1 < argc) mem_gb       = std::atof(argv[++i]);
        else if (arg == "--min-count" && i + 1 < argc) min_count    = std::atoi(argv[++i]);
        else { std::cerr << "Unknown argument: " << arg << "\n"; PrintUsage(argv[0]); return 1; }
    }

    if (!reads2_file.empty() && reads_file.empty()) {
        std::cerr << "Error: --reads2 requires --reads (R1 file)\n";
        PrintUsage(argv[0]);
        return 1;
    }
    if (output_path.empty() || (graph_prefix.empty() == reads_file.empty())) {
        // both empty or both set
        PrintUsage(argv[0]);
        return 1;
    }

#ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
    std::cerr << "Threads: " << omp_get_max_threads() << "\n";
#else
    std::cerr << "Threads: 1 (OpenMP not available)\n";
#endif

    // ----- reads mode: build SDBG, remember to clean up -----
    bool owns_tmp = false;
    std::string tmp_dir;

    if (!reads_file.empty()) {
        const std::string bin_dir = GetBinaryDir(argv[0]);
        tmp_dir   = bin_dir + PATH_SEP + "megahit_topo_tmp_" + std::to_string(
#ifdef _WIN32
                        static_cast<unsigned>(GetCurrentProcessId())
#else
                        static_cast<unsigned>(getpid())
#endif
                    );
        owns_tmp = true;
        try {
            graph_prefix = BuildSdbgFromReads(reads_file, reads2_file,
                                              bin_dir, tmp_dir,
                                              n_threads > 0 ? n_threads : 1,
                                              mem_gb, min_count);
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << "\n";
            if (owns_tmp) RemoveDir(tmp_dir);
            return 1;
        }
    }

    // ----- load SDBG -----
    SDBG dbg;
    try {
        std::cerr << "Loading SDBG from: " << graph_prefix << " ...\n";
        dbg.LoadFromFile(graph_prefix.c_str());
        std::cerr << "  edges: " << dbg.size() << "  k=" << dbg.k() << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error loading SDBG: " << e.what() << "\n";
        if (owns_tmp) RemoveDir(tmp_dir);
        return 1;
    }

    // ----- extract -----
    TopoFeatures features;
    try {
        std::cerr << "Extracting features (sample_size=" << opts.sample_size << ") ...\n";
        features = RunExtraction(dbg, opts);
    } catch (const std::exception& e) {
        std::cerr << "Error during extraction: " << e.what() << "\n";
        if (owns_tmp) RemoveDir(tmp_dir);
        return 1;
    }

    std::cerr << "  histogram:  " << features.timing.histogram_ms << " ms\n";
    std::cerr << "  node:       " << features.timing.node_ms      << " ms\n";
    std::cerr << "  tip walk:   " << features.timing.walk_ms      << " ms\n";
    std::cerr << "  bubbles:    " << features.timing.bubble_ms    << " ms\n";

    // ----- write output -----
    try {
        WriteFeatures(features, output_path);
        std::cerr << "Features written to: " << output_path << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error writing output: " << e.what() << "\n";
        if (owns_tmp) RemoveDir(tmp_dir);
        return 1;
    }

    // ----- clean up temp SDBG -----
    if (owns_tmp) {
        std::cerr << "Removing temp SDBG: " << tmp_dir << "\n";
        RemoveDir(tmp_dir);
    }

    return 0;
}
