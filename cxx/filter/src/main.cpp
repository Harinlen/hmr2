#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "hmr_args.h"
#include "hmr_args_type.h"
#include "hmr_ui.h"
#include "filter_enzyme.h"
#include "filter_fasta_type.h"
#include "filter_fasta.h"
#include "filter_bam.h"

#include "hmr_fasta.h"

extern HMR_ARGS opts;

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Print the current working parameter.
    time_print("Execution configuration:");
    time_print("\tMinimum Map Quality: %d", opts.mapq);
    time_print("\tThreads: %d", opts.threads);
    //Read the FASTA file and find the enzyme.
    const char *nuc_seq = NULL;
    int nuc_seq_size = 0;
    filter_enzyme_formalize(opts.enzyme, &nuc_seq, &nuc_seq_size);
    opts.nuc_seq = nuc_seq;
    opts.nuc_seq_length = nuc_seq_size;
    time_print("\tEnzyme: %s", opts.nuc_seq);
    time_print("Building enzyme index from FASTA file...");
    //Initial the FASTA filter user.
    FASTA_FILTER_USER fasta_user;
    FASTA_ENZYME_POSES &enzyme_poses = fasta_user.enzyme_poses;
    {
        //Prepare the workers.
        FILTER_WORKERS workers(filter_fasta_search_enzyme);
        fasta_user.workers = &workers;
        fasta_user.enzyme_poses = FASTA_ENZYME_POSES {0, NULL, NULL, NULL, NULL, NULL};
        //Start parsing FASTA.
        fasta_parser(opts.fasta, filter_fasta_push_work, &fasta_user);
        //Wait until all the workers are completed.
        workers.wait_for_tasks();
    }
    time_print("Enzyme index built, total sequenece: %d", enzyme_poses.seq_total);
    char node_path[1024];
    sprintf(node_path, "%s.node", opts.output);
    time_print("Dump node info to %s", node_path);
    filter_fasta_dump_node(node_path, enzyme_poses);
    time_print("Graph nodes dump completed.");
    //Read the filter the mapping file.
    time_print("Filtering mapping file %s", opts.mapping);
    char reads_path[1024];
    sprintf(reads_path, "%s.reads", opts.output);
    time_print("Dump reads position info to %s", reads_path);
    std::vector<READ_EDGE> edges;
    filter_bam_statistic(opts.mapping, fasta_user.enzyme_poses, opts.mapq, opts.threads, edges, reads_path);
    time_print("Filter complete, %zu edges generated.", edges.size());
    //Based on the node and raw edges, dump the weighted edges.
    char edge_path[1024];
    sprintf(edge_path, "%s.edge", opts.output);
    time_print("Dump edge info to %s", edge_path);
    filter_bam_dump_edge(edge_path, edges);
    time_print("Graph edges dump completed.");
    return 0;
}
