#ifndef ARGS_CORRECT_H
#define ARGS_CORRECT_H

#include <vector>

typedef struct HMR_ARGS
{
    const char *fasta = NULL;
    const char *output = NULL;
    char *enzyme = NULL;
    std::vector<char *> left, right;
    int threads = 1;
} HMR_ARGS;

#endif // ARGS_CORRECT_H
