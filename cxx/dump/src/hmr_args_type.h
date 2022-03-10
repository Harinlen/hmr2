#ifndef HMR_ARGS_TYPE_H
#define HMR_ARGS_TYPE_H

typedef struct HMR_ARGS {
    char *nodes = nullptr;
    char *node_rename = nullptr;
    char *cluster = nullptr;
    char *tour = nullptr;
    char *fasta = nullptr;
    char *output = nullptr;
} HMR_ARGS;

#endif // HMR_ARGS_TYPE_H
