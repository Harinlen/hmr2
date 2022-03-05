#include <cstdlib>
#include <cstdio>
#include <algorithm>

#include "hmr_args.h"
#include "hmr_args_type.h"
#include "hmr_cluster.h"
#include "hmr_ui.h"
#include "hmr_global.h"
#include "hmr_reads_parse.h"
#include "hmr_cgd_parse.h"
#include "optimize_order.h"
#include "optimize_direction.h"
#include "hmr_tour.h"

extern HMR_ARGS opts;

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Print the current working parameter.
    time_print("Execution configuration:");
    time_print("\tThreads: %d", opts.threads);
    //Read and store the network.
    CONTIG_NODE *nodes;
    size_t node_size;
    time_print("Reading network node and edge files...");
    hmr_net_read(opts.node, opts.edge, &nodes, &node_size);
    time_print("%zu node(s) read.", node_size);
    //Read the partition records.
    time_print("Reading clustering result...");
    std::vector<int> node_id_group = hmr_cluster_read(opts.cluster);
    time_print("Cluster result read, %zu nodes readed.", node_id_group.size());
    //Calculate the sequence.
    ORDER_NODE_INFO *node_infos = new ORDER_NODE_INFO[node_size];
    time_print("Optimize the order of the cluster...");
    SEQUENCE node_sequence = optimize_order(node_id_group, nodes, node_infos);
    delete[] node_infos;
    time_print("Sequence determined.");
    //Determine the direction of each group.
    time_print("Reading HiC reads file...");
    POLAR_INFO polar_infos = hmr_polar_parse(opts.reads, node_sequence, nodes);
    time_print("Group reads imported.");
    time_print("Deciding directions...");
    optimize_direction(polar_infos);
    time_print("Direction optimized complete.");
    //Write the tour file for the result.
    hmr_tour_dump(opts.output, nodes, polar_infos);
    time_print("Optimize result dumps to %s", opts.output);
    return 0;
}
