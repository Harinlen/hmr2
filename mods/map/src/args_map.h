#ifndef ARGS_MAP_H
#define ARGS_MAP_H

#include <vector>

typedef struct HMR_ARGS
{
    const char *fasta = NULL;
    const char *output = NULL;
    char *enzyme = NULL;
    std::vector<char *> left, right;
    int threads = 1;
} HMR_ARGS;

#endif // ARGS_MAP_H
