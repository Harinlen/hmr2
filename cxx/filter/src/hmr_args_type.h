#ifndef HMR_ARGS_TYPE_H
#define HMR_ARGS_TYPE_H

typedef struct HMR_ARGS {
    char *mapping = nullptr;
    char *fasta = nullptr;
    char *output = nullptr;
    char *enzyme = nullptr;
    const char *nuc_seq = nullptr;
    int nuc_seq_length = 0;
    int mapq = 40;
    int threads = 1;
} HMR_ARGS;

#endif // HMR_ARGS_TYPE_H
