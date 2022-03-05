#ifndef GROUP_CONV_H
#define GROUP_CONV_H

#include <unordered_set>
#include <list>
#include <vector>

typedef struct CONTIG_NODE CONTIG_NODE;

void group_hanamaru(CONTIG_NODE *nodes, size_t node_size, int no_of_group, int threads,
                    std::vector<std::vector<int> > &groups, std::list<int> &unknown_ids);

#endif // GROUP_CONV_H
