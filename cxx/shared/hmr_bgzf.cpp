#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <thread>

#include <libdeflate.h>
#include "hmr_ui.h"

#include "hmr_bgzf.h"

#define WORK_PER_THREAD (512)

typedef struct BGZF_HEADER
{
    uint8_t ID1;
    uint8_t ID2;
    uint8_t CM;
    uint8_t FLG;
    uint32_t MTIME;
    uint8_t XFL;
    uint8_t OS;
    uint16_t XLEN;
} BGZF_HEADER;

typedef struct BGZF_SUB_HEADER
{
    uint8_t SI1;
    uint8_t SI2;
    uint16_t SLEN;
} BGZF_SUB_HEADER;

typedef struct BGZF_FOOTER
{
    uint32_t CRC32;
    uint32_t ISIZE;
} BGZF_FOOTER;

typedef struct HMR_BGZF_DECOMPRESS
{
    char *cdata;
    uint16_t cdata_size;
    size_t offset;
    size_t raw_size;
} HMR_BGZF_DECOMPRESS;

void hmr_bgzf_queue_create(HMR_BGZF_QUEUE *queue, size_t size)
{
    //Create the slices buffer.
    queue->slices = static_cast<HMR_BGZF_SLICE *>(malloc(sizeof(HMR_BGZF_SLICE) * size));
    queue->size = size;
    queue->head = 0;
    queue->tail = 0;
    queue->finish = false;
}

void hmr_bgzf_queue_free(HMR_BGZF_QUEUE *queue)
{
    //Clear the slices.
    free(queue->slices);
}

HMR_BGZF_SLICE hmr_bgzf_pop(HMR_BGZF_QUEUE *queue)
{
    //Wait until the queue is not empty.
    std::unique_lock<std::mutex> pop_lock(queue->mutex);
    queue->pop_cv.wait(pop_lock, [queue]
    {
        return (queue->head != queue->tail) || queue->finish;
    });
    //Extract the data.
    HMR_BGZF_SLICE slice = queue->slices[queue->head];
    if(queue->finish)
    {
        slice.data = NULL;
        slice.data_size = 0;
    }
    else
    {
        queue->head = (queue->head+1 == queue->size) ? 0 : (queue->head+1);
    }
    //Notify the push variable.
    queue->push_cv.notify_one();
    return slice;
}

void hmr_bgzf_push(HMR_BGZF_QUEUE *queue, char *raw_data, size_t raw_data_size)
{
    //Check whether the queue is full.
    std::unique_lock<std::mutex> push_lock(queue->mutex);
    queue->push_cv.wait(push_lock, [queue]
    {
        return !((queue->tail+1==queue->head) || (queue->head==0&&queue->tail==queue->size - 1));
    });
    //Push the data to the queue.
    queue->slices[queue->tail] = HMR_BGZF_SLICE {raw_data, raw_data_size};
    queue->tail = (queue->tail+1 == queue->size) ? 0 : (queue->tail+1);
    //Notify the conditional variable.
    queue->pop_cv.notify_one();
}

void hmr_bgzf_decompress(int id, HMR_BGZF_DECOMPRESS *pool, char *bgzf_raw)
{
    struct libdeflate_decompressor *z = libdeflate_alloc_decompressor();
    size_t data_size;
    //Loop and decompress the data.
    for(int i=id*WORK_PER_THREAD, target=(id+1)*WORK_PER_THREAD; i<target; ++i)
    {
        //Decompress the data.
        libdeflate_deflate_decompress(z,
                                      pool[i].cdata, pool[i].cdata_size,
                                      bgzf_raw + pool[i].offset, pool[i].raw_size,
                                      &data_size);
        //Recover the compress data memory.
        free(pool[i].cdata);
    }
    libdeflate_free_decompressor(z);
}

void hmr_bgzf_parse(const char *filepath, HMR_BGZF_QUEUE *queue, int threads)
{
    //Read the BAM file.
#ifdef _MSC_VER
    FILE *bgzf_file = NULL;
    fopen_s(&bgzf_file, filepath, "rb");
#else
    FILE *bgzf_file = fopen(filepath, "rb");
#endif
    if(!bgzf_file)
    {
        time_error_str(1, "Failed to open bgzf file %s", filepath);
    }
    //Get the total file size.
    fseek(bgzf_file, 0L, SEEK_END);
#ifdef _MSC_VER
    size_t total_size = _ftelli64(bgzf_file);
#else
    size_t total_size = ftello64(bgzf_file);
#endif
    fseek(bgzf_file, 0L, SEEK_SET);
    size_t report_size = total_size / 10, report_pos = report_size;
    //Prepare the decompression buffer.
    size_t block_offset = 0;
    BGZF_HEADER header_buf;
    BGZF_FOOTER footer_buf;
    uint32_t buf_size = WORK_PER_THREAD * threads, buf_used = 0;
    HMR_BGZF_DECOMPRESS *bgzf_buf = static_cast<HMR_BGZF_DECOMPRESS *>(malloc(sizeof(HMR_BGZF_DECOMPRESS) * buf_size));
    //Read while to the end of the file.
    while(fread(&header_buf, sizeof(BGZF_HEADER), 1, bgzf_file) > 0)
    {
        //Read the Xlen data.
        char *subfield_data = static_cast<char *>(malloc(header_buf.XLEN));
        //Read the data.
        fread(subfield_data, header_buf.XLEN, 1, bgzf_file);
        //Go through the header.
        uint16_t subfield_left = header_buf.XLEN, bsize = 0;
        char *subfield_pos = subfield_data;
        while(subfield_left > 0)
        {
            BGZF_SUB_HEADER *subfield = reinterpret_cast<BGZF_SUB_HEADER *>(subfield_pos);
            //Check the ID matches the bsize.
            if(subfield->SI1 == 66 && subfield->SI2 == 67 && subfield->SLEN == 2)
            {
                bsize = *(reinterpret_cast<uint16_t *>(subfield_pos + sizeof(BGZF_SUB_HEADER)));
            }
            subfield_pos += subfield->SLEN + sizeof(BGZF_SUB_HEADER);
            subfield_left -= subfield->SLEN + sizeof(BGZF_SUB_HEADER);
        }
        free(subfield_data);
        uint16_t cdata_size = bsize - header_buf.XLEN - 19;
        char *cdata = static_cast<char *>(malloc(cdata_size));
        //Reading the compressed data.
        fread(cdata, cdata_size, 1, bgzf_file);
        //Fetch the footer data.
        fread(&footer_buf, sizeof(BGZF_FOOTER), 1, bgzf_file);
        //Buffer until reach the decompress limit.
        bgzf_buf[buf_used] = HMR_BGZF_DECOMPRESS { cdata, cdata_size, block_offset, footer_buf.ISIZE };
        ++buf_used;
        //Increase the block offset, and keep going.
        block_offset += footer_buf.ISIZE;
        //Check whether we are reaching the decompression limitation.
        if(buf_used == buf_size)
        {
            //Create the pool.
            char *bgzf_raw = static_cast<char *>(malloc(block_offset));
            //Decompress the data.
            std::thread *workers = new std::thread[threads];
            for(int i=0; i<threads; ++i)
            {
                workers[i] = std::thread(hmr_bgzf_decompress, i, bgzf_buf, bgzf_raw);
            }
            for(int i=0; i<threads; ++i)
            {
                workers[i].join();
            }
            delete[] workers;
            //Push the data to the parsing queue.
            hmr_bgzf_push(queue, bgzf_raw, block_offset);
            //Reset the buffer used.
            buf_used = 0;
            block_offset = 0;
        }
        //Check should we report the position.
#ifdef _MSC_VER
        size_t bgzf_pos = _ftelli64(bgzf_file);
#else
        size_t bgzf_pos = ftello64(bgzf_file);
#endif
        if(bgzf_pos >= report_pos)
        {
            float percent = static_cast<float>(bgzf_pos) / static_cast<float>(total_size) * 100.0f;
            time_print_float("BGZF parsed %.1f%%", percent);
            report_pos += report_size;
        }
    }
    //Mark BGZF complete.
    {
        std::unique_lock<std::mutex> finish_lock(queue->mutex);
        queue->finish = true;
        queue->pop_cv.notify_one();
    }
    fclose(bgzf_file);
}
