#ifndef FILTER_BAM_H
#define FILTER_BAM_H

#include <vector>

typedef struct FASTA_ENZYME_POSES FASTA_ENZYME_POSES;

typedef struct READ_EDGE
{
    int32_t start;
    int32_t end;
    double weight;
} READ_EDGE;

std::vector<READ_EDGE> filter_bam_statistic(const char *source, const FASTA_ENZYME_POSES &poses, int mapq, int threads);

void filter_bam_dump_edge(const char *filepath, const std::vector<READ_EDGE> &edges);

#endif // FILTER_BAM_H
