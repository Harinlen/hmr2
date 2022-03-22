#ifndef HMR_ENZYME_TYPE_H
#define HMR_ENZYME_TYPE_H

#include <cstdint>

typedef struct ENZYME_RANGE
{
    size_t start, end;
} ENZYME_RANGE;

typedef struct CONTIG_ENZYME_RANGE
{
    size_t range_size;
    ENZYME_RANGE* ranges;
};

#endif // HMR_ENZYME_TYPE_H