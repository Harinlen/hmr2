#ifndef HMR_FASTQ_TYPE_H
#define HMR_FASTQ_TYPE_H

typedef struct FASTQ_READ
{
    size_t name_length, seq_length, split_length;
    char* name, *seq, *quality, *split;
} FASTQ_READ;

typedef struct FASTQ_READ_PAIR
{
    FASTQ_READ left, right;
} FASTQ_READ_PAIR;

#endif // HMR_FASTQ_TYPE_H