#include <cstdlib>
#include <cstdio>

#include "hmr_ui.h"
#include "hmr_cluster.h"

std::vector<int> hmr_cluster_read(const char *filepath)
{
#ifdef _MSC_VER
    FILE *cluster_file = NULL;
    fopen_s(&cluster_file, filepath, "rb");
#else
    FILE *cluster_file = fopen(filepath, "rb");
#endif
    if(cluster_file == NULL)
    {
        time_error(-1, "Failed to read cluster node file %s", filepath);
    }
    //Prepare the value set.
    std::vector<int> nodes;
    int32_t num_of_nodes, node_id;
    fread(&num_of_nodes, sizeof(int32_t), 1, cluster_file);
    nodes.reserve(num_of_nodes);
    for(int j=0; j<num_of_nodes; ++j)
    {
        fread(&node_id, sizeof(int32_t), 1, cluster_file);
        nodes.push_back(node_id);
    }
    fclose(cluster_file);
    return nodes;
}
