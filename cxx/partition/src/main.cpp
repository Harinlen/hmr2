#include <cstdlib>
#include <cstdio>
#include <algorithm>

#include "hmr_args.h"
#include "hmr_args_type.h"
#include "hmr_ui.h"
#include "group_conv.h"

#include "hmr_cgd_parse.h"

extern HMR_ARGS opts;

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Print the current working parameter.
    time_print("Execution configuration:");
    time_print_int("\tPartition group: %d", opts.group);
    time_print_int("\tThreads: %d", opts.threads);
    //Read and store the network.
    CONTIG_NODE *nodes;
    size_t node_size;
    time_print("Reading network node and edge files...");
    hmr_net_read(opts.node, opts.edge, &nodes, &node_size);
    time_print_size("%zu node(s) read.", node_size);
    //Run hana-maru algorithm.
    std::vector<std::set<int> > groups = group_hanamaru(nodes, node_size, opts.group, opts.threads);
    //Write the partition result to the output file.
    time_print_str("Dumping partition result to %s", opts.output);
#ifdef _MSC_VER
    FILE *output_file = NULL;
    fopen_s(&output_file, opts.output, "wb");
#else
    FILE *output_file = fopen(opts.output, "wb");
#endif
    if(output_file == NULL)
    {
        time_error_str(-1, "Failed to open output file %s", opts.output);
    }
    //Write the total group number.
    int result_groups = static_cast<int>(groups.size());
    fwrite(&result_groups, sizeof(int), 1, output_file);
    for(int i=0; i<result_groups; ++i)
    {
        //Write the size of the group first.
        const auto node_ids = groups[i];
        int group_size = static_cast<int>(node_ids.size());
        fwrite(&group_size, sizeof(int), 1, output_file);
        //Write the node id.
        for(int j: node_ids)
        {
            fwrite(&j, sizeof(int), 1, output_file);
        }
    }
    fclose(output_file);
    time_print("Partition complete.");
    return 0;
}
