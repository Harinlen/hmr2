#ifndef FILTER_FASTA_TYPE_H
#define FILTER_FASTA_TYPE_H

#include <vector>

#include "hmr_thread_pool.h"

typedef struct ENZYME_RANGE
{
    size_t start;
    size_t end;
} ENZYME_RANGE;

typedef struct FASTA_ENZYME_POSES
{
    size_t seq_total;
    char **seq_names;
    size_t *seq_length;
    size_t *enz_range_size;
    ENZYME_RANGE **enz_ranges;
} FASTA_ENZYME_POSES;

typedef struct FASTA_FILTER_WORK_ARGS
{
    char *seq_name;
    size_t seq_name_length;
    char *seq;
    size_t seq_size;
    FASTA_ENZYME_POSES *enzyme_poses;
} FASTA_FILTER_WORK_ARGS;

typedef thread_pool<void (const FASTA_FILTER_WORK_ARGS &), FASTA_FILTER_WORK_ARGS> FILTER_WORKERS;

typedef struct
{
    FILTER_WORKERS *workers;
    FASTA_ENZYME_POSES enzyme_poses;
} FASTA_FILTER_USER;

#endif // FILTER_FASTA_TYPE_H
