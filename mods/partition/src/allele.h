#ifndef ALLELE_H
#define ALLELE_H

#include <vector>

#include "hmr_contig_graph_type.h"

typedef std::vector<int32_t> ALLELE_IDS;
typedef std::vector<ALLELE_IDS> ALLELE_TABLE;

ALLELE_TABLE allele_load(const char* filepath, const HMR_CONTIGS& contigs, int32_t allele_groups);

#endif // ALLELE_H