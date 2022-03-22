#ifndef MAPPING_CORRECT_TYPE_H
#define MAPPING_CORRECT_TYPE_H

#include "contig_correct_type.h"

typedef struct BAM_CONTIG_MAP
{
    uint32_t size;
    int32_t* id;
} BAM_CONTIG_MAP;

typedef union POS_PAIR
{
    struct {
        int32_t a, b;
    } pos;
    uint64_t data;
} POS_PAIR;

typedef std::unordered_map<uint64_t, int32_t> HIC_DB;

typedef struct BAM_CORRECT_MAP
{
    HIC_DB *wide_db, *narrow_db;
    CONTIG_MAP contig_map;
    BAM_CONTIG_MAP bam_id_map;
    int32_t bam_contig_id;
    int32_t wide, narrow;
    uint8_t mapq;
} BAM_CORRECT_MAP;

#endif // MAPPING_CORRECT_H