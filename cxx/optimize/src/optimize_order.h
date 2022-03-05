#ifndef OPTIMIZE_ORDER_H
#define OPTIMIZE_ORDER_H

#include <map>
#include <vector>

typedef std::vector<int> SEQUENCE;

typedef struct ORDER_EDGE
{
    double weight;
    int32_t order;
} ORDER_EDGE;

typedef std::map<int32_t, ORDER_EDGE> ORDER_EDGE_MAP;

typedef struct ORDER_NODE_INFO
{
    SEQUENCE *parent;
    double penalty;
    int32_t prev_order, next_order;
    ORDER_EDGE_MAP edge_map;
} ORDER_NODE_INFO;

typedef struct CONTIG_NODE CONTIG_NODE;

SEQUENCE optimize_order(const std::vector<int32_t> &node_ids, CONTIG_NODE *nodes, ORDER_NODE_INFO *node_infos);

#endif // HMR_OPTIMIZE_H
