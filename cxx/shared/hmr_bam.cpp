#include "hmr_ui.h"
#include "hmr_bgzf.h"

#include "hmr_bam.h"

typedef struct BAM_DATA_BUF
{
    char *data;
    size_t size, reserve, offset;
} BAM_DATA_BUF;

char *hmr_bam_fetch(BAM_DATA_BUF *buf, HMR_BGZF_QUEUE *queue, size_t size)
{
    //Check the left data is enough.
    size_t residual = buf->size - buf->offset;
    if(residual >= size)
    {
        //Enough to hold the data.
        char *data = buf->data + buf->offset;
        buf->offset += size;
        return data;
    }
    //Well, we have to fill the buffer with the queue data size.
    while(residual < size)
    {
        //Fetch the data from the queue.
        HMR_BGZF_SLICE bam_slice = hmr_bgzf_pop(queue);
        if(bam_slice.data == NULL)
        {
            //Reach the end of the data.
            return NULL;
        }
        //Increase the bam slice.
        if(residual == 0)
        {
            //Directly use the slice memory.
            if(buf->data)
            {
                free(buf->data);
            }
            buf->data = bam_slice.data;
            buf->size = bam_slice.data_size;
            buf->reserve = bam_slice.data_size;
        }
        else
        {
            //Copy the residual data.
            memcpy(buf->data, buf->data + buf->offset, residual);
            //Calculate the new buffer size.
            size_t data_size = residual + bam_slice.data_size;
            //Update the buffer when necessary.
            if(data_size > buf->reserve)
            {
                buf->data = static_cast<char *>(realloc(buf->data, data_size));
                buf->reserve = data_size;
            }
            //Copy the slice data to buffer.
            memcpy(buf->data + residual, bam_slice.data, bam_slice.data_size);
            //Update the data size.
            buf->size = data_size;
            //Release the slice data.
            free(bam_slice.data);
        }
        //Update the residual data.
        buf->offset = 0;
        residual = buf->size;
    }
    //Now the data should be enough.
    char *data = buf->data;
    buf->offset = size;
    return data;
}

inline uint32_t hmr_fetch_uint32(BAM_DATA_BUF *buf, HMR_BGZF_QUEUE *queue)
{
    return *(reinterpret_cast<uint32_t *>(hmr_bam_fetch(buf, queue, 4)));
}

void hmr_bam_load(const char *filepath, BAM_PARSER *parser, int threads, void *user)
{
    //Prepare the BGZF parsing queue.
    HMR_BGZF_QUEUE queue;
    hmr_bgzf_queue_create(&queue, threads);
    //Prepare the BAM buffer.
    BAM_DATA_BUF buf{ NULL, 0, 0, 0 };
    //Start the BGZF parsing thread.
    std::thread bgzf_parse(hmr_bgzf_parse, filepath, &queue, threads);
    //Fetch and check the magic number.
    char *magic = hmr_bam_fetch(&buf, &queue, 4);
    if(strncmp(magic, "BAM\1", 4))
    {
        time_error(1, "BAM header magic string incorrect.");
    }
    //Read the header text.
    uint32_t l_text = hmr_fetch_uint32(&buf, &queue);
    char *text = hmr_bam_fetch(&buf, &queue, l_text);
    parser->header(text, l_text, user);
    //Fetch the n_ref.
    uint32_t n_ref = hmr_fetch_uint32(&buf, &queue);
    parser->n_ref(n_ref, user);
    //Loop until all the reference are parsed.
    while(n_ref--)
    {
        uint32_t l_name = hmr_fetch_uint32(&buf, &queue);
        char *name = hmr_bam_fetch(&buf, &queue, l_name);
        uint32_t l_ref = hmr_fetch_uint32(&buf, &queue);
        parser->ref(name, l_name, l_ref, user);
    }
    //Fetch the rest of the data.
    char *block_size_data = hmr_bam_fetch(&buf, &queue, 4);
    size_t block_id = 0;
    while(block_size_data)
    {
        uint32_t block_size = *(reinterpret_cast<uint32_t *>(block_size_data));
        //Fetch the data of the fetch.
        char *block_data = hmr_bam_fetch(&buf, &queue, block_size);
        //Parse the block data.
        parser->block(block_id, block_data, block_size, user);
        //Increase the block id.
        ++block_id;
        //Fetch the next block.
        block_size_data = hmr_bam_fetch(&buf, &queue, 4);
    }
    //Release the buf.
    if(buf.data)
    {
        free(buf.data);
    }
    //Wait the parse thread.
    bgzf_parse.join();
    hmr_bgzf_queue_free(&queue);
}
