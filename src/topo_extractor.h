#pragma once

#include <cstdint>
#include <string>
#include "sdbg/sdbg.h"
#include "topo_features.h"

// Options passed to RunExtraction.
struct ExtractionOptions {
    uint64_t    total_reads = 0;   // total read count (R1+R2); 0 = not provided

    // Jellyfish k-mer feature extraction (all required together; leave empty to skip).
    std::string ref_fasta;   // reference FASTA — used as ground truth for error/true labels
    std::string reads1;      // R1 (or SE) reads file for jellyfish count
    std::string reads2;      // R2 reads file; empty for single-end
    std::string tmp_dir;     // directory for jellyfish intermediate files (must already exist)
    int         kmer_k      = 21;
    int         n_threads   = 1;
};

// Extract graph topology features in a single parallel pass over all SDBG edges.
NodeFeatures ExtractNodeFeatures(SDBG& dbg, uint64_t total_reads);

// Extract k-mer features by running jellyfish on the reference and reads.
// Writes jellyfish files inside opts.tmp_dir, optimises thresholds for labels_out.
KmerFeatures ExtractKmerFeatures(const ExtractionOptions& opts, Labels& labels_out);

// Top-level driver: runs both extractions and returns a combined result.
TopoFeatures RunExtraction(SDBG& dbg, const ExtractionOptions& opts);
