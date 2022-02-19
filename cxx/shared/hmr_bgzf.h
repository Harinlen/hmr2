#ifndef HMR_BGZF_H
#define HMR_BGZF_H

#include <cstdint>
#include <condition_variable>

typedef struct HMR_BGZF_SLICE
{
    char *data;
    size_t data_size;
} HMR_BGZF_SLICE;

typedef struct HMR_BGZF_QUEUE
{
    bool finish;
    HMR_BGZF_SLICE *slices;
    size_t head, tail, size;
    std::mutex mutex;
    std::condition_variable pop_cv, push_cv;
} HMR_BGZF_QUEUE;

void hmr_bgzf_queue_create(HMR_BGZF_QUEUE *queue, size_t size);
void hmr_bgzf_queue_free(HMR_BGZF_QUEUE *queue);
HMR_BGZF_SLICE hmr_bgzf_pop(HMR_BGZF_QUEUE *queue);

void hmr_bgzf_parse(const char *filepath, HMR_BGZF_QUEUE *queue, int threads);

#endif // HMR_BGZF_H
