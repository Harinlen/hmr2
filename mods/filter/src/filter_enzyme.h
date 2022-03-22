#ifndef FILTER_ENZYME_H
#define FILTER_ENZYME_H

#include <mutex>
#include <vector>

#include "hmr_fasta_type.h"
#include "hmr_fastq_type.h"
#include "hmr_thread_pool.h"
#include "filter_enzyme_type.h"

typedef struct FILTER_ENZYME_RANGE_DATA
{
    std::vector<FILTER_ENZYME_RANGE> ranges;
    char* seq;
    size_t seq_size;
} FILTER_ENZYME_RANGE_DATA;

typedef struct FILTER_ENZYME_RANGES
{
    std::vector<HMR_CONTIG>* contigs;
    std::vector<FILTER_ENZYME_RANGE_DATA>* enzyme_ranges;
    std::mutex* access_mutex;
} FILTER_ENZYME_RANGES;

typedef struct FILTER_ENZYME_SEARCH
{
    const char* enzyme;
    size_t enzyme_size;
} FILTER_ENZYME_SEARCH;

typedef struct FILTER_ENZYME_ARGS
{
    char* seq_name;
    size_t seq_name_length;
    char* seq;
    size_t seq_size;
    FILTER_ENZYME_SEARCH* search;
    FILTER_ENZYME_RANGES* ranges;
} FILTER_ENZYME_ARGS;

typedef hmr::thread_pool<FILTER_ENZYME_ARGS> FILTER_ENZYME_WORKERS;

typedef struct FILTER_ENZYME_USER
{
    FILTER_ENZYME_WORKERS* workers;
    FILTER_ENZYME_SEARCH* search;
    FILTER_ENZYME_RANGES* ranges;
} FILTER_ENZYME_USER;

void filter_enzyme_search_prepare(const char *enzyme, size_t enzyme_size, FILTER_ENZYME_SEARCH *search);

void filter_enzyme_search_submit(char* seq_name, size_t seq_name_length, char* seq, size_t seq_size, void* user);
void filter_enzyme_search(const FILTER_ENZYME_ARGS& args);

typedef struct FILTER_READ_PAIR_ARGS
{
    FASTQ_READ_PAIR read_pair;
} FILTER_READ_PAIR_ARGS;

typedef hmr::thread_pool<FILTER_READ_PAIR_ARGS> FILTER_READ_PAIR_WORKERS;

typedef struct FILTER_READ_PAIR_USER
{
    FILTER_READ_PAIR_WORKERS *workers;
} FILTER_READ_PAIR_USER;

void filter_read_matching_submit(const FASTQ_READ_PAIR& read_pair, void *user);
void filter_read_matching(const FILTER_READ_PAIR_ARGS& args);

#endif // FILTER_ENZYME_H