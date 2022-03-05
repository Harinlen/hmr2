#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <vector>
#include <list>

#include "hmr_ui.h"
#include "hmr_cgd_parse.h"

typedef struct READ_EDGE
{
    int32_t start;
    int32_t end;
    double weight;
} READ_EDGE;


void hmr_node_read(const char *node_path, CONTIG_NODE **graph_nodes, size_t *graph_node_size)
{
    //Read the vertices.
#ifdef _MSC_VER
    FILE *node_file = NULL;
    fopen_s(&node_file, node_path, "rb");
#else
    FILE *node_file = fopen(node_path, "rb");
#endif
    if(node_file == NULL)
    {
        time_error(-1, "Failed to read node file %s", node_path);
    }
    size_t node_size;
    fread(&node_size, sizeof(size_t), 1, node_file);
    *graph_node_size = node_size;
    //Create the node size.
    CONTIG_NODE *nodes = new CONTIG_NODE[node_size];
    *graph_nodes = nodes;
    for(size_t i=0; i<node_size; ++i)
    {
        fread(&(nodes[i].l_name), sizeof(size_t), 1, node_file);
        nodes[i].name = new char[nodes[i].l_name + 1];
        fread(nodes[i].name, nodes[i].l_name, 1, node_file);
        nodes[i].name[nodes[i].l_name] = '\0';
        fread(&(nodes[i].length), sizeof(size_t), 1, node_file);
    }
    fclose(node_file);
}

int hmr_net_read(const char *node_path, const char *edge_path, CONTIG_NODE **graph_nodes, size_t *graph_node_size)
{
    //Read the nodes.
    hmr_node_read(node_path, graph_nodes, graph_node_size);
    CONTIG_NODE *nodes = *graph_nodes;
    //Read the edges.
#ifdef _MSC_VER
    FILE *edge_file = NULL;
    fopen_s(&edge_file, edge_path, "rb");
#else
    FILE *edge_file = fopen(edge_path, "rb");
#endif
    if(edge_file == NULL)
    {
        time_error(-1, "Failed to read edge file %s", edge_path);
    }
    size_t edge_size;
    READ_EDGE edge;
    fread(&edge_size, sizeof(size_t), 1, edge_file);
    while(edge_size--)
    {
        fread(&edge, sizeof(READ_EDGE), 1, edge_file);
        //Ignore the edge point to itself.
        if(edge.start == edge.end)
        {
            continue;
        }
        //Add the edge data to nodes.
        nodes[edge.start].links.push_back(CONTIG_EDGE{edge.end, edge.weight});
        nodes[edge.end].links.push_back(CONTIG_EDGE{edge.start, edge.weight});
    }
    fclose(edge_file);
    return 0;
}
