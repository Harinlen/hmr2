//Harerun ID: 826917359
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <set>
#include <unordered_set>

#include "hmr_cgd_type.h"
#include "hmr_global.h"

#include "optimize_order.h"

typedef struct RELATION
{
    double weight;
    double mark;
} RELATION;

typedef union EDGE_DATA
{
    struct EDGE
    {
        int32_t start;
        int32_t end;
    };

    EDGE edge;
    uint64_t data;
} EDGE_DATA;

inline void get_edge(const uint64_t &data, int32_t &x, int32_t &y)
{
    EDGE_DATA edge_data;
    edge_data.data = data;
    x = edge_data.edge.start;
    y = edge_data.edge.end;
}

inline uint64_t edge_data(int32_t x, int32_t y)
{
    EDGE_DATA edge_data;
    if(x < y)
    {
        edge_data.edge.start = x;
        edge_data.edge.end = y;
    }
    else
    {
        edge_data.edge.start = y;
        edge_data.edge.end = x;
    }
    return edge_data.data;
}

typedef std::unordered_map<uint64_t, RELATION> GROUP_RELATION;

typedef struct EDGE_RELATION
{
    uint64_t edge;
    double weight;
    double mark;
} EDGE_RELATION;

inline double edge_penualty(uint8_t edge)
{
    if(edge == 0x12 || edge == 0x13)
    {
        return 18.0;
    }
//    if(edge == 0x24 || edge == 0x25 || edge == 0x34 || edge == 0x35)
//    {
//        return 7.0;
//    }
    return 3.0;
}

inline void add_relation(int32_t x, int32_t y, double weight, uint8_t edge, GROUP_RELATION &relation)
{
    if(weight == 0.0)
    {
        return;
    }
    //Decide the edge mark.
    uint64_t xy_edge = edge_data(x, y);
    auto finder = relation.find(xy_edge);
    if(finder == relation.end())
    {
        relation.insert(std::make_pair(xy_edge, RELATION {weight, edge_penualty(edge)}));
    }
    else
    {
        finder->second.mark += edge_penualty(edge);
    }
}

inline double edge_weight(int32_t x, int32_t y, ORDER_NODE_INFO *node_infos)
{
    const auto &edge_map = node_infos[x].edge_map;
    auto finder = edge_map.find(y);
    return (finder == edge_map.end()) ? 0.0 : finder->second.weight;
}

SEQUENCE optimize_order(const std::vector<int32_t> &node_ids, CONTIG_NODE *nodes,
                        ORDER_NODE_INFO *node_infos)
{
    //Build the node id set.
    std::unordered_set<int32_t> group_nodes;
    //Generate the node related.
    for(const int32_t node_id: node_ids)
    {
        group_nodes.insert(node_id);
    }
    GROUP_RELATION group_relation;
    for(const int32_t node_id: node_ids)
    {
        auto &node_info = node_infos[node_id];
        auto &node = nodes[node_id];
        std::sort(node.links.begin(), node.links.end(), [](const CONTIG_EDGE &lhs, const CONTIG_EDGE &rhs) {
            return lhs.weight > rhs.weight;
        });
        //Update the node info id.
        node_info.node_id = node_id;
        //Extract the edges within the node info.
        std::vector<CONTIG_EDGE> edges;
        edges.reserve(node.links.size());
        for(const auto &edge: node.links)
        {
            if(hInSet(edge.id, group_nodes))
            {
                edges.push_back(edge);
            }
        }
        edges.reserve(edges.size());
        node.links = edges;
        //Create the edge map.
        ORDER_EDGE_MAP edge_map;
        for(size_t i=0; i<edges.size(); ++i)
        {
            const auto &edge = edges[i];
            edge_map.insert(std::make_pair(edge.id, ORDER_EDGE {edge.weight, static_cast<int32_t>(i)}));
        }
        node_info.edge_map = edge_map;
    }

    ;
    exit(-1);
}
