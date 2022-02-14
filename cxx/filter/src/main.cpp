#include <cstdlib>
#include <cstdio>

#include "hmr_args.h"
#include "hmr_args_type.h"
#include "filter_enzyme.h"
#include "filter_fasta_type.h"
#include "filter_fasta.h"
#include "hmr_ui.h"

#include "hmr_fasta.h"

extern HMR_ARGS opts;

static const char *edge_suffix = ".edge";

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Print the current working parameter.
    time_print("Execution configuration:");
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
        fasta_user.enzyme_poses = FASTA_ENZYME_POSES {0, NULL, NULL, NULL, NULL};
        //Start parsing FASTA.
        fasta_parser(opts.fasta, filter_fasta_push_work, &fasta_user);
        //Wait until all the workers are completed.
        workers.wait_for_tasks();
    }
    time_print_size("Enzyme index built, total sequenece: %d", enzyme_poses.seq_total);
    //Dump the enzyme count file.
    {
        char enzyme_count_path[2048];
        sprintf(enzyme_count_path, "%s.enz_count", opts.fasta);
        time_print_str("Dump enzyme index to %s", enzyme_count_path);
        filter_fasta_dump_enzyme_count(enzyme_count_path, enzyme_poses);
        time_print("Enzyme index dumped.");
    }
    //Read the filter the mapping file.
    time_print_str("Filtering mapping file %s", opts.mapping);
    return 0;
}
