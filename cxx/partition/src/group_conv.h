#ifndef GROUP_CONV_H
#define GROUP_CONV_H

#include <set>
#include <vector>

typedef struct CONTIG_NODE CONTIG_NODE;

std::vector<std::set<int> > group_hanamaru(CONTIG_NODE *nodes, size_t node_size, int groups, int threads);

#endif // GROUP_CONV_H
