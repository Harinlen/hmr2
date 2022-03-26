#ifndef FASTA_DRAFT_TYPE_H
#define FASTA_DRAFT_TYPE_H

#include <cstdint>

typedef struct ENZYME_RANGE
{
    int32_t start, end;
} ENZYME_RANGE;

typedef struct ENZYME_RANGES
{
    ENZYME_RANGE* ranges;
    size_t length, counter;
} ENZYME_RANGES;

#endif // FASTA_DRAFT_TYPE_H