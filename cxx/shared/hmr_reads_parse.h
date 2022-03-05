#ifndef HMR_READS_PARSE_H
#define HMR_READS_PARSE_H

#include "hmr_polar_type.h"

typedef struct CONTIG_NODE CONTIG_NODE;

POLAR_INFO hmr_polar_parse(const char *filepath, const std::vector<int> &sequence, CONTIG_NODE *nodes);

#endif // HMR_READS_PARSE_H
