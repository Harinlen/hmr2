#ifndef PARTITION_TYPE_H
#define PARTITION_TYPE_H

#include <vector>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

typedef struct CONTIG_EDGE
{
    int32_t id;
    double weight;
} CONTIG_EDGE;

typedef std::vector<CONTIG_EDGE> CONTIG_EDGE_LIST;
typedef std::vector<CONTIG_EDGE_LIST> CONTIG_EDGES;
typedef std::unordered_map<int32_t, double> CONTIG_EDGE_MAP;

typedef std::unordered_set<int32_t> CONTIG_ID_SET;

#endif // PARTITION_TYPE_H