#ifndef PARTITION_H
#define PARTITION_H

#include "partition_type.h"
#include "hmr_contig_graph_type.h"

void partition_load_edges(const char* filepath, CONTIG_EDGES& edges);

std::vector<CONTIG_ID_SET> partition_run(const HMR_CONTIGS& contigs, CONTIG_EDGES& edges, size_t num_of_group, const int32_t& threads);

#endif // PARTITION_H