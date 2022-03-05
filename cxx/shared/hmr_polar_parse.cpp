#include <cstdio>
#include <cstdlib>

#include "hmr_ui.h"

#include "hmr_polar_parse.h"

void hmr_polar_parse_half(FILE *polar_file, HALF_EDGE &half_edge)
{
    //Read the node size.
    size_t half_size, edge_count;
    int32_t edge_id;
    fread(&half_size, sizeof(size_t), 1, polar_file);
    //Read until all the data are done.
    while(half_size--)
    {
        fread(&edge_id, sizeof(int32_t), 1, polar_file);
        fread(&edge_count, sizeof(size_t), 1, polar_file);
        //Insert the edge count info.
        half_edge.insert(std::make_pair(edge_id, edge_count));
    }
}

POLAR_INFO hmr_polar_parse(const char *filepath)
{
    POLAR_INFO graph_half_edges;
#ifdef _MSC_VER
    FILE *polar_file = NULL;
    fopen_s(&polar_file, filepath, "rb");
#else
    FILE *polar_file = fopen(filepath, "rb");
#endif
    if(polar_file == NULL)
    {
        time_error(-1, "Failed to read contig polarity file %s", filepath);
    }
    //Read the contig size of the file.
    size_t contig_size;
    int32_t node_id;
    fread(&contig_size, sizeof(size_t), 1, polar_file);
    while(contig_size--)
    {
        //Read the node id.
        fread(&node_id, sizeof(int32_t), 1, polar_file);
        //Read the left side info.
        BAM_CONTIG_HALF node_info;
        hmr_polar_parse_half(polar_file, node_info.left);
        hmr_polar_parse_half(polar_file, node_info.right);
        //Insert the node.
        graph_half_edges.insert(std::make_pair(node_id, node_info));
    }
    fclose(polar_file);
    return graph_half_edges;
}
