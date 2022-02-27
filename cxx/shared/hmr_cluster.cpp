#include "hmr_ui.h"
#include "hmr_cluster.h"


void hmr_cluster_read(const char *filepath, std::vector<std::vector<int> > &groups)
{
#ifdef _MSC_VER
    FILE *cluster_file = NULL;
    fopen_s(&cluster_file, filepath, "rb");
#else
    FILE *cluster_file = fopen(filepath, "rb");
#endif
    if(cluster_file == NULL)
    {
        time_error_str(-1, "Failed to cluster node file %s", filepath);
    }
    //Read the size of the group.
    int num_of_group, num_of_nodes, node_id;
    fread(&num_of_group, sizeof(int), 1, cluster_file);
    groups.reserve(num_of_group);
    for(int i=0; i<num_of_group; ++i)
    {
        //Prepare the value set.
        std::vector<int> nodes;
        fread(&num_of_nodes, sizeof(int), 1, cluster_file);
        nodes.reserve(num_of_nodes);
        for(int j=0; j<num_of_nodes; ++j)
        {
            fread(&node_id, sizeof(int), 1, cluster_file);
            nodes.push_back(node_id);
        }
        groups.push_back(nodes);
    }
    fclose(cluster_file);
}
