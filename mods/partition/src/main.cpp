#include <cstdlib>
#include <cstdio>

#include "hmr_args.h"
#include "hmr_ui.h"
#include "hmr_path.h"
#include "hmr_contig_graph.h"
#include "partition.h"
#include "allele.h"

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
    if (opts.allele_groups > -1 && opts.allele_groups < 2) { time_error(-1, "Allele groups must be greater than 2."); }
    if (opts.allele_groups > 0 && (!opts.allele_table)) { time_error(-1, "Allele table must be provided for allele group division."); }
    //Print the execution configuration.
    time_print("Execution configuration:");
    time_print("\tNumber of Partitions: %d", opts.groups);
    time_print("\tThreads: %d", opts.threads);
    time_print("\tAllele mode: %s", opts.allele_groups > 0 ? "Yes" : "No");
    //Load the contig node information.
    HMR_CONTIGS contigs;
    time_print("Loading contig information from %s", opts.nodes);
    hmr_graph_load_contig(opts.nodes, contigs);
    time_print("%zu contig(s) information loaded.", contigs.size());
    //Loading the allele table if needed.
    ALLELE_TABLE allele_table;
    if (opts.allele_groups > 0)
    {
        time_print("Loading allele table %s", opts.allele_table);
        allele_table = allele_load(opts.allele_table, contigs, opts.allele_groups);
        time_print("%zu allele id rule(s) loaded.", allele_table.size());
    }
    //Read the edge information.
    time_print("Loading edge information from %s", opts.edge);
    CONTIG_EDGES edges;
    edges.resize(contigs.size());
    partition_load_edges(opts.edge, edges);
    time_print("Contig edges loaded.");
    //Partition mission start.
    time_print("Dividing contigs into %d groups...", opts.groups);
    auto partition_result = partition_run(contigs, edges, static_cast<size_t>(opts.groups), opts.threads);
    time_print("%zu group(s) of contigs generated.", partition_result.size());
    //Check whether we have to divide them into allele groups.
    if (opts.allele_groups > 0)
    {
        time_print("Divided into %d allele group(s)...", opts.allele_groups);
        ;
    }
    return 0;
}