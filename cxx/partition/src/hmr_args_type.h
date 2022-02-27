#ifndef HMR_ARGS_TYPE_H
#define HMR_ARGS_TYPE_H

typedef struct HMR_ARGS {
    char *node = nullptr;
    char *edge = nullptr;
    char *output = nullptr;
    int group = -1;
    int threads = 1;
} HMR_ARGS;

#endif // HMR_ARGS_TYPE_H
