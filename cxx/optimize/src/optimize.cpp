//Harerun ID: 826917359
#include <algorithm>
#include <set>

#include "hmr_cgd_parse.h"

#include "optimize.h"

typedef struct OPTIMIZE_NODE
{
    int id;
    std::set<int> closest;
    std::vector<CONTIG_EDGE> edges, trust_edges;
} OPTIMIZE_NODE;

void optimize_find_nearest()
{
    ;
}

void optimize_group(std::vector<int> &node_ids, CONTIG_NODE *nodes, size_t node_size)
{
    //Decide the convlution range, half of the size.
    int detect_range = static_cast<int>(node_ids.size()) >> 1;
    //Build the relation matrix.
    std::set<int> node_set(node_ids.begin(), node_ids.end());
    std::vector<OPTIMIZE_NODE> group_nodes;
    group_nodes.reserve(node_ids.size());
    for(int id: node_ids)
    {
        //Find the node that inside the group.
        const CONTIG_NODE &node_info = nodes[id];
        OPTIMIZE_NODE node;
        node.id = id;
        node.edges.reserve(node_ids.size());
        for(auto edge: node_info.links)
        {
            //Save the edge if we can find the node.
            if(node_set.find(edge.id) != node_set.end())
            {
                node.edges.push_back(edge);
            }
        }
        //Sort the edges.
        std::sort(node.edges.begin(), node.edges.end(), [](const auto &lhs, const auto &rhs){
            return lhs.weight > rhs.weight;
        });
        //Push the trusted edges.
        if(static_cast<int>(node.edges.size()) < detect_range)
        {
            node.trust_edges = node.edges;
        }
        else
        {
            node.trust_edges.insert(node.trust_edges.end(), node.edges.begin(), node.edges.begin() + detect_range);
        }
        int closest_size = 4;
        if(node.trust_edges.size() < 4)
        {
            closest_size = static_cast<int>(node.trust_edges.size());
        }
        node.closest.insert(id);
        for(int i=0; i<closest_size; ++i)
        {
            node.closest.insert(node.trust_edges[i].id);
        }
        group_nodes.push_back(node);
    }
    //Print the edge size.
    for(size_t i=0; i<group_nodes.size(); ++i)
    {
        printf("%s\n", nodes[group_nodes[i].id].name);
        for(size_t j=0; j<group_nodes[i].closest.size()-1; ++j)
        {
            printf("%s %lf\t", nodes[group_nodes[i].trust_edges[j].id].name, group_nodes[i].trust_edges[j].weight);
        }
        printf("\n");
    }
}
