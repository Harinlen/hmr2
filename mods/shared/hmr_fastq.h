#ifndef HMR_FASTQ_H
#define HMR_FASTQ_H

#include "hmr_fastq_type.h"

typedef void (*FASTQ_PAIR_PROC)(const FASTQ_READ_PAIR &read_pair, void* user);

void hmr_fastq_pair_read(const char* left_filepath, const char *right_filepath, FASTQ_PAIR_PROC parser, void *user);
void hmr_fastq_read_free(const FASTQ_READ& read);

#endif // HMR_FASTQ_H