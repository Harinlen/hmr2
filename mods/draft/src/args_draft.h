#ifndef ARGS_DRAFT_H
#define ARGS_DRAFT_H

#include <vector>

typedef struct HMR_ARGS
{
    const char *fasta = NULL;
    const char *output = NULL;
    std::vector<char *> mappings;
    char* enzyme = nullptr;
    const char* enzyme_nuc = nullptr;
    int enzyme_nuc_length = 0, mapq = 40, threads = 1;
} HMR_ARGS;

#endif // ARGS_DRAFT_H
