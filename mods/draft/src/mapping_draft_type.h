#ifndef MAPPING_DRAFT_TYPE_H
#define MAPPING_DRAFT_TYPE_H

#include <unordered_map>

#include "hmr_contig_graph_type.h"

typedef std::unordered_map<std::string, int32_t> CONTIG_ID_MAP;
typedef std::unordered_map<uint64_t, uint64_t> RAW_EDGE_MAP;

typedef union MAPPING_READ
{
    struct {
        int32_t id;
        int32_t pos;
    } read;
    uint64_t data;
} MAPPING_READ;

typedef std::unordered_map<uint64_t, MAPPING_READ> READ_RECORD;

typedef struct MAPPING_DRAFT_USER
{
    READ_RECORD records;
    const CONTIG_ID_MAP& contig_ids;
    const HMR_CONTIG_INVALID_SET& invalid_ids;
    ENZYME_RANGES* contig_ranges;
    int32_t *contig_id_map;
    int32_t contig_idx;
    RAW_EDGE_MAP edges;
    FILE* reads_file;
    uint8_t mapq;
    char* output_buffer;
    size_t output_offset, output_size;
} MAPPING_DRAFT_USER;

#endif // MAPPING_DRAFT_H