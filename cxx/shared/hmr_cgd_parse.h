//Harerun ID: 826917359
#ifndef HMR_NET_PARSE_H
#define HMR_NET_PARSE_H

#include <cstdint>
#include <unordered_map>

typedef struct CONTIG_EDGE
{
    int32_t id;
    double weight;
} CONTIG_EDGE;

typedef struct CONTIG_NODE
{
    size_t l_name;
    char *name;
    std::vector<CONTIG_EDGE> links;
} CONTIG_NODE;

int hmr_net_read(const char *node_path, const char *edge_path, CONTIG_NODE **graph_nodes, size_t *graph_node_size);

#endif // HMR_NET_PARSE_H
