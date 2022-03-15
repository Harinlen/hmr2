#include <cstdlib>
#include <cstdio>
#include <algorithm>

#include "hmr_args.h"
#include "hmr_args_type.h"
#include "hmr_ui.h"
#include "partition.h"

#include "hmr_cgd_parse.h"

extern HMR_ARGS opts;

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Print the current working parameter.
    time_print("Execution configuration:");
    time_print("\tPartition group: %d", opts.group);
    time_print("\tThreads: %d", opts.threads);
    //Read and store the network.
    CONTIG_NODE *nodes;
    size_t node_size;
    time_print("Reading network node and edge files...");
    hmr_net_read(opts.node, opts.edge, &nodes, &node_size);
    time_print("%zu node(s) read.", node_size);
    //Run hana-maru algorithm.
    std::vector<std::vector<int> > groups;
    std::list<int> unknown_ids;
    group_hanamaru(nodes, node_size, opts.group, opts.threads, groups, unknown_ids);
    //Check defected ids.
    if(!unknown_ids.empty())
    {
        time_print("%zu unknown nodes detected.", unknown_ids.size());
        char defect_path[1024];
        sprintf(defect_path, "%s.defect", opts.output);
        time_print("Dump defect nodes to %s", defect_path);
#ifdef _MSC_VER
        FILE *defect_file = NULL;
        fopen_s(&defect_file, defect_path, "wb");
#else
        FILE *defect_file = fopen(defect_path, "wb");
#endif
        int defect_size = static_cast<int>(unknown_ids.size());
        fwrite(&defect_size, sizeof(int), 1, defect_file);
        for(int i: unknown_ids)
        {
            fwrite(&i, sizeof(int), 1, defect_file);
        }
        fclose(defect_file);
        time_print("Defected nodes dumped.");
    }
    std::sort(groups.begin(), groups.end(), [](const std::vector<int> &lhs, const std::vector<int> &rhs) {
        return lhs.size() > rhs.size();
    });
    //Write the summary file.
    char output_path[1024];
    sprintf(output_path, "%s_cluster.txt", opts.output);
#ifdef _MSC_VER
    FILE *summary_file = NULL;
    fopen_s(&summary_file, output_path, "w");
#else
    FILE *summary_file = fopen(output_path, "w");
#endif
    for(int i=0; i<opts.group; ++i)
    {
        if(i > 0)
        {
            fprintf(summary_file, "\n");
        }
        fprintf(summary_file, "Group %2d:\n", i);
        for(int node_id: groups[i])
        {
            fprintf(summary_file, "%s\n", nodes[node_id].name);
        }
    }
    if(!unknown_ids.empty())
    {
        fprintf(summary_file, "\nDefected:\n");
        for(int node_id: unknown_ids)
        {
            fprintf(summary_file, "%s\n", nodes[node_id].name);
        }
    }
    fclose(summary_file);
    //Write the partition result to the output file.
    for(int i=0; i<opts.group; ++i)
    {
        sprintf(output_path, "%s_%dg%d.cluster", opts.output, i+1, opts.group);
#ifdef _MSC_VER
        FILE *output_file = NULL;
        fopen_s(&output_file, output_path, "wb");
#else
        FILE *output_file = fopen(output_path, "wb");
#endif
        if(output_file == NULL)
        {
            time_error(-1, "Failed to open output file %s", opts.output);
        }
        //Write the size of the group first.
        const auto node_ids = groups[i];
        time_print("Dumping partition %2d (%zu nodes)...", i, node_ids.size());
        int group_size = static_cast<int>(node_ids.size());
        fwrite(&group_size, sizeof(int), 1, output_file);
        //Write the node id.
        for(int j: node_ids)
        {
            fwrite(&j, sizeof(int), 1, output_file);
        }
        fclose(output_file);
    }
    time_print("Partition info output complete.");
    return 0;
}
