#include <cstdlib>
#include <cstdio>

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
    time_print_int("\tMinimum Map Quality: %d", opts.mapq);
    time_print_int("\tThreads: %d", opts.threads);
    //Read the FASTA file and find the enzyme.
    const char *nuc_seq = NULL;
    int nuc_seq_size = 0;
    filter_enzyme_formalize(opts.enzyme, &nuc_seq, &nuc_seq_size);
    opts.nuc_seq = nuc_seq;
    opts.nuc_seq_length = nuc_seq_size;
    time_print_str("\tEnzyme: %s", opts.nuc_seq);
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
    time_print_size("Enzyme index built, total sequenece: %d", enzyme_poses.seq_total);
    //Dump the enzyme count file.
#ifdef _MSC_VER
    FILE *output_file = NULL;
    fopen_s(&output_file, opts.output, "w");
#else
    FILE *output_file = fopen(opts.output, "w");
#endif
    time_print_str("Dump node info to %s", opts.output);
    filter_fasta_dump_node(output_file, enzyme_poses);
    time_print("Graph nodes dump completed.");
    //Read the filter the mapping file.
    time_print_str("Filtering mapping file %s", opts.mapping);
    std::unordered_map<uint64_t, size_t> raw_edges =
            filter_bam_statistic(opts.mapping, fasta_user.enzyme_poses, opts.mapq, opts.threads);
    time_print_size("Filter complete, %zu edges generated.", raw_edges.size());
    //Based on the node and raw edges, dump the weighted edges.
    time_print_str("Dump edge info to %s", opts.output);
    filter_bam_dump_edge(output_file, raw_edges, fasta_user.enzyme_poses);
    time_print("Graph edges dump completed.");
    fclose(output_file);
    return 0;
}
