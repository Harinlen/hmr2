#ifndef HMR_BIN_ENZYME_H
#define HMR_BIN_ENZYME_H

#include <vector>

#include "hmr_fasta_type.h"
#include "hmr_enzyme_type.h"

void hmr_enzyme_formalize(char* enzyme, const char** nuc_seq, int* nuc_seq_size);

bool hmr_enzyme_load_range(const char* filepath, std::vector<CONTIG_ENZYME_RANGE>& contig_ranges);
bool hmr_enzyme_save_range(const char* filepath, const std::vector<CONTIG_ENZYME_RANGE>& contig_ranges);

#endif // HMR_BIN_ENZYME