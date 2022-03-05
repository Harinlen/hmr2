#ifndef FILTER_BAM_H
#define FILTER_BAM_H

#include <unordered_map>
#include <vector>

typedef struct FASTA_ENZYME_POSES FASTA_ENZYME_POSES;

typedef struct READ_EDGE
{
    int32_t start;
    int32_t end;
    double weight;
} READ_EDGE;

void filter_bam_statistic(const char *source, const FASTA_ENZYME_POSES &poses, int mapq, int threads,
                          std::vector<READ_EDGE> &graph_edges, const char *reads_path);

void filter_bam_dump_edge(const char *filepath, const std::vector<READ_EDGE> &edges);

#endif // FILTER_BAM_H
