#ifndef HMR_CGD_TYPE_H
#define HMR_CGD_TYPE_H

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
    size_t length;
    std::vector<CONTIG_EDGE> links;
} CONTIG_NODE;

#endif // HMR_CGD_TYPE_H
