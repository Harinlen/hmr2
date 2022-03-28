#ifndef FASTA_DRAFT_H
#define FASTA_DRAFT_H

#include "hmr_contig_graph_type.h"
#include "hmr_thread_pool.h"

#include "fasta_draft_type.h"

typedef struct ENZYME_RANGE_CHAIN
{
    ENZYME_RANGES data;
    struct ENZYME_RANGE_CHAIN* next;
} ENZYME_RANGE_CHAIN;

typedef struct ENZYME_SEARCH
{
    const char* enzyme;
    int32_t enzyme_length;
    int32_t* kmpNext;
} ENZYME_SEARCH;

typedef struct ENZYME_RANGE_SEARCH
{
    ENZYME_SEARCH* search;
    ENZYME_RANGE_CHAIN* chain_node;
    char* seq;
    int32_t seq_size, range;
} ENZYME_RANGE_SEARCH;

void contig_range_search(const ENZYME_RANGE_SEARCH& param);

typedef hmr::thread_pool<ENZYME_RANGE_SEARCH> RANGE_SEARCH_POOL;

typedef struct DRAFT_NODES_USER
{
    const int32_t range;
    HMR_CONTIGS* nodes;
    ENZYME_SEARCH* search;
    RANGE_SEARCH_POOL* pool;
    ENZYME_RANGE_CHAIN* chain_head;
    ENZYME_RANGE_CHAIN* chain_tail;
} DRAFT_NODES_USER;

void contig_draft_search_start(const char* enzyme, int32_t enzyme_length, ENZYME_SEARCH& search);
void contig_draft_search_end(ENZYME_SEARCH& search);
int32_t contig_draft_search(const char *seq, size_t seq_size, ENZYME_SEARCH* search);

void contig_draft_build(int32_t index, char* seq_name, size_t seq_name_size, char* seq, size_t seq_size, void* user);

#endif // FASTA_DRAFT_H