#ifndef CONTIG_CORRECT_H
#define CONTIG_CORRECT_H

#include <cstdint>

void contig_correct_build(int32_t index, char* seq_name, size_t seq_name_size, char* seq, size_t seq_size, void* user);

#endif // CONTIG_CORRECT_H