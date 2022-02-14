#ifndef HMR_FASTA_H
#define HMR_FASTA_H

typedef void (*FASTA_SEQ_PROC)(char *seq_name, size_t seq_name_length, char *seq, size_t seq_size, void *user);

void fasta_parser(const char *filepath, FASTA_SEQ_PROC parser, void *user);

#endif // HMR_FASTA_H
