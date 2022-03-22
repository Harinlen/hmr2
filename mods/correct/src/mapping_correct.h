#ifndef MAPPING_CORRECT_H
#define MAPPING_CORRECT_H

#include "hmr_mapping_type.h"

void mapping_correct_n_contig(uint32_t n_ref, void* user);
void mapping_correct_contig(uint32_t name_length, char* name, uint32_t length, void* user);
void mapping_correct_read_align(size_t id, const MAPPING_INFO& mapping_info, void* user);

#endif // MAPPING_CORRECT_H