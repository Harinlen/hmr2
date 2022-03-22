#include <cstdlib>
#include <cstdio>

#include "hmr_args.h"
#include "hmr_enzyme.h"
#include "hmr_ui.h"

#include "args_map.h"
#include "map_fasta.h"

extern HMR_ARGS opts;

int main(int argc, char* argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.fasta) { help_exit(-1, "No FASTA file path."); }
    if (opts.left.empty() || opts.right.empty()) { help_exit(-1, "Missing at least one side pair FASTQ path."); }
    if (opts.left.size() != opts.right.size()) { help_exit(-1, "HiC reads FASTQ are not paired."); }
    if (opts.enzyme == NULL) { help_exit(-1, "Missing enzyme."); }
    //Print the working parameters.
    time_print("Execution configuration:");
    time_print("\tThreads: %d", opts.threads);
    //Formalize the enzyme.
    const char* enzyme_seq = NULL;
    int enzyme_seq_size = 0;
    hmr_enzyme_formalize(opts.enzyme, &enzyme_seq, &enzyme_seq_size);
    time_print("\tEnzyme: %s", enzyme_seq);
    //Load the FASTA data, including: enzyme ranges and BW table (bwt and sa).
    map_fasta_load(opts.fasta);
    return 0;
}