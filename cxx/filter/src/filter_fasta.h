#ifndef FILTER_FASTA_INDEX_H
#define FILTER_FASTA_INDEX_H

#include <vector>

typedef struct FASTA_ENZYME_POSES FASTA_ENZYME_POSES;
typedef struct FASTA_FILTER_WORK_ARGS FASTA_FILTER_WORK_ARGS;

void filter_fasta_search_enzyme(const FASTA_FILTER_WORK_ARGS &args);

void filter_fasta_push_work(char *seq_name, size_t seq_name_length, char *seq, size_t seq_size, void *user);

void filter_fasta_dump_node(FILE *graph_file, const FASTA_ENZYME_POSES &poses);

#endif // FILTER_FASTA_INDEX_H
