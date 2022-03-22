#ifndef HMR_FASTA_TYPE_H
#define HMR_FASTA_TYPE_H

#include <cstdint>

typedef struct HMR_CONTIG
{
    char* seq_name;
    size_t seq_name_length, seq_length;
} HMR_CONTIG;

#endif // HMR_FASTA_TYPE_H