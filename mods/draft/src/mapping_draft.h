#ifndef MAPPING_DRAFT_H
#define MAPPING_DRAFT_H

#include "hmr_mapping_type.h"
#include "mapping_draft_type.h"

void mapping_draft_n_contig(uint32_t n_ref, void* user);
void mapping_draft_contig(uint32_t name_length, char* name, uint32_t length, void* user);
void mapping_draft_read_align(size_t id, const MAPPING_INFO& mapping_info, void* user);

std::vector<HMR_EDGE_WEIGHT> mapping_draft_get_edge_weights(const RAW_EDGE_MAP &edge_map, const HMR_CONTIGS &contigs);

#endif // MAPPING_DRAFT_H