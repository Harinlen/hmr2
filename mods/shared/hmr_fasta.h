#ifndef HMR_FASTA_H
#define HMR_FASTA_H

#include <string>
#include <vector>

#include "hmr_fasta_type.h"

typedef void (*FASTA_PROC)(int32_t, char *, size_t , char *, size_t , void *);

void hmr_fasta_read(const char *filepath, FASTA_PROC parser, void *user);

std::string hmr_fasta_path_contig(const char* filepath);
bool hmr_fasta_load_contig(const char* filepath, std::vector<HMR_CONTIG>& contigs);
bool hmr_fasta_save_contig(const char* filepath, const std::vector<HMR_CONTIG>& contigs);

#endif // HMR_FASTA_H
