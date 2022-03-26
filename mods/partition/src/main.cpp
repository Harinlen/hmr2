#include <cstdlib>
#include <cstdio>

#include "hmr_args.h"
#include "hmr_ui.h"
#include "hmr_path.h"
#include "hmr_contig_graph.h"
#include "partition.h"

#include "args_partition.h"

extern HMR_ARGS opts;

int main(int argc, char* argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.nodes) { help_exit(-1, "Missing HMR graph contig information file path."); }
    if (!path_can_read(opts.nodes)) { time_error(-1, "Cannot read HMR graph contig file %s", opts.nodes); }
    if (!opts.edge) { help_exit(-1, "Missing HMR graph edge weight file path."); }
    if (!path_can_read(opts.edge)) { time_error(-1, "Cannot read HMR graph edge weight file %s", opts.edge); }
    if (opts.groups < 1) { time_error(-1, "Please specify the group to be separated."); }
    //Print the execution configuration.
    time_print("Execution configuration:");
    time_print("\tNumber of Partitions: %d", opts.groups);
    time_print("\tHas allele: %s", opts.allele == NULL ? "Yes" : "No");
    time_print("\tThreads: %d", opts.threads);
    //Load the contig node information.
    HMR_CONTIGS contigs;
    time_print("Loading contig information from %s", opts.nodes);
    hmr_graph_load_contig(opts.nodes, contigs);
    time_print("%zu contig(s) information loaded.", contigs.size());
    //Read the edge information.
    time_print("Loading edge information from %s", opts.edge);
    //Based on the allele table construct the nodes.
    HMR_NODES nodes;
    //Build 
    return 0;
}