#ifndef HMR_BAM_H
#define HMR_BAM_H

#include <cstdint>

typedef struct BAM_BLOCK_HEADER
{
    int32_t refID;
    int32_t pos;
    uint8_t l_read_name;
    uint8_t mapq;
    uint16_t bin;
    uint16_t n_cigar_op;
    uint16_t flag;
    uint32_t l_seq;
    int32_t next_refID;
    int32_t next_pos;
    int32_t tlen;
} BAM_BLOCK_HEADER;

typedef void (*BAM_PARSE_HEADER)(char *text, uint32_t l_text, void *user);
typedef void (*BAM_PARSE_N_REF)(uint32_t n_ref, void *user);
typedef void (*BAM_PARSE_REF)(char *name, uint32_t l_name, uint32_t l_ref, void *user);
typedef void (*BAM_PARSE_BLOCK)(size_t block_id, char *block, uint32_t block_size, void *user);

typedef struct BAM_PARSER
{
    BAM_PARSE_HEADER header;
    BAM_PARSE_N_REF n_ref;
    BAM_PARSE_REF ref;
    BAM_PARSE_BLOCK block;
} BAM_PARSER;

void hmr_bam_load(const char *filepath, BAM_PARSER *parser, int threads, void *user);

#endif // HMR_BAM_H
