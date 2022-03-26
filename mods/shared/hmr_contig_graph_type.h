#ifndef HMR_CONTIG_GRAPH_TYPE_H
#define HMR_CONTIG_GRAPH_TYPE_H

#include <cstdint>
#include <vector>

typedef struct HMR_CONTIG
{
    int32_t name_size;
    char* name;
    int32_t length;
} HMR_CONTIG;

typedef std::vector<HMR_CONTIG> HMR_CONTIGS;

typedef union HMR_EDGE
{
    struct {
        int32_t start;
        int32_t end;
    } pos;
    uint64_t data;
} HMR_EDGE;

typedef struct HMR_EDGE_WEIGHT
{
    HMR_EDGE edge;
    double weight;
} HMR_EDGE_WEIGHT;

typedef std::vector<HMR_EDGE_WEIGHT> HMR_EDGE_WEIGHTS;

typedef struct HMR_MAPPING
{
    int32_t refID;
    int32_t pos;
    int32_t next_refID;
    int32_t next_pos;
} HMR_MAPPING;

#endif // HMR_CONTIG_GRAPH_TYPE_H