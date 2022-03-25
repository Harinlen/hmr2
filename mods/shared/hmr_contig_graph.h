#ifndef HMR_CONTIG_GRAPH_H
#define HMR_CONTIG_GRAPH_H

#include <string>

#include "hmr_contig_graph_type.h"

HMR_EDGE hmr_graph_edge(int32_t a, int32_t b);

std::string hmr_graph_path_contig(const char* prefix);
bool hmr_graph_load_contig(const char* filepath, HMR_CONTIGS& contigs);
bool hmr_graph_save_contig(const char* filepath, const HMR_CONTIGS& contigs);

std::string hmr_graph_path_reads(const char* prefix);

std::string hmr_graph_path_edge(const char* prefix);
bool hmr_graph_save_edge(const char* filepath, const HMR_EDGE_WEIGHTS& edges);

#endif // HMR_CONTIG_GRAPH_H