#include <cstdlib>
#include <cstdio>
#include <algorithm>

#include "hmr_args.h"
#include "hmr_args_type.h"
#include "hmr_cluster.h"
#include "hmr_ui.h"
#include "optimize.h"

#include "hmr_cgd_parse.h"

extern HMR_ARGS opts;

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Print the current working parameter.
    time_print("Execution configuration:");
    time_print_int("\tThreads: %d", opts.threads);
    //Read and store the network.
    CONTIG_NODE *nodes;
    size_t node_size;
    time_print("Reading network node and edge files...");
    hmr_net_read(opts.node, opts.edge, &nodes, &node_size);
    time_print_size("%zu node(s) read.", node_size);
    //Read the partition records.
    time_print("Reading clustering result...");
    std::vector<std::vector<int> > groups;
    hmr_cluster_read(opts.clusters, groups);
    time_print("Cluster result read.");
    //Optimize the first group.
    optimize_group(groups[0], nodes, node_size);
    return 0;
}
