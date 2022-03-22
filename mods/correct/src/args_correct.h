#ifndef ARGS_CORRECT_H
#define ARGS_CORRECT_H

#include <vector>

typedef struct HMR_ARGS
{
    const char *fasta = NULL;
    const char *output = NULL;
    std::vector<char *> mappings;
    double percent = 0.95, sensitive = 0.5;
    int mapq = 1, wide = 25000, narrow = 1000, depletion = 100000, threads = 1;
} HMR_ARGS;

#endif // ARGS_CORRECT_H
