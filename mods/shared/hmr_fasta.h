#ifndef HMR_FASTA_H
#define HMR_FASTA_H

typedef void (*FASTA_PROC)(int32_t, char *, size_t , char *, size_t , void *);

void hmr_fasta_read(const char *filepath, FASTA_PROC parser, void *user);

#endif // HMR_FASTA_H
