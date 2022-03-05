//Harerun ID: 826917359
#ifndef HMR_NET_PARSE_H
#define HMR_NET_PARSE_H

#include "hmr_cgd_type.h"

void hmr_node_read(const char *node_path, CONTIG_NODE **graph_nodes, size_t *graph_node_size);

int hmr_net_read(const char *node_path, const char *edge_path, CONTIG_NODE **graph_nodes, size_t *graph_node_size);

#endif // HMR_NET_PARSE_H
