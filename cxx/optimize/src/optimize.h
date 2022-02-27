#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <vector>

typedef struct CONTIG_NODE CONTIG_NODE;

void optimize_group(std::vector<int> &node_ids, CONTIG_NODE *nodes, size_t node_size);

#endif // HMR_OPTIMIZE_H
