#ifndef FILTER_BAM_H
#define FILTER_BAM_H

#include <unordered_map>

typedef struct FASTA_ENZYME_POSES FASTA_ENZYME_POSES;

std::unordered_map<uint64_t, size_t> filter_bam_statistic(const char *source, const FASTA_ENZYME_POSES &poses, int mapq, int threads);

void filter_bam_dump_edge(FILE *graph_file, const std::unordered_map<uint64_t, size_t> &raw_edges, const FASTA_ENZYME_POSES &poses);

#endif // FILTER_BAM_H
