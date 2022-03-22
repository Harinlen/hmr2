#include <cstdlib>
#include <cstdio>

#include "hmr_args.h"
#include "hmr_fasta.h"
#include "hmr_fastq.h"
#include "hmr_enzyme.h"
#include "hmr_ui.h"

#include "args_filter.h"
#include "filter_enzyme.h"

extern HMR_ARGS opts;

int main(int argc, char *argv[])
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
    //Create the fasta information file and find the enzyme inside the sequences.
    time_print("Finding enzyme ranges in FASTA file %s", opts.fasta);
    std::vector<HMR_CONTIG> contigs;
    std::vector<FILTER_ENZYME_RANGE_DATA> enzyme_ranges;
    {
        //Do search preprocessing for enzyme.
        FILTER_ENZYME_SEARCH enzyme_search;
        filter_enzyme_search_prepare(enzyme_seq, enzyme_seq_size, &enzyme_search);
        //Build the search structures.
        std::mutex access_mutex;
        FILTER_ENZYME_RANGES search_ranges {&contigs, &enzyme_ranges, &access_mutex};
        FILTER_ENZYME_WORKERS workers(filter_enzyme_search, opts.threads << 4, opts.threads);
        //Start parsing the FASTA.
        FILTER_ENZYME_USER enzyme_search_user { &workers, &enzyme_search, &search_ranges };
        hmr_fasta_read(opts.fasta, filter_enzyme_search_submit, &enzyme_search_user);
    }
    for (size_t i = 0; i < enzyme_ranges.size(); ++i)
    {
        printf("%s\t%zu\t%zu\n", contigs[i].seq_name, enzyme_ranges[i].ranges.size() * 1000, contigs[i].seq_length);
    }
    time_print("%zu contig(s) processed.", contigs.size());
    //Read the pair-ends data.
    /*time_print("Filtering the pair-ends files:");
    for (size_t i = 0; i < opts.left.size(); ++i)
    {
        const char* left_path = opts.left[i], * right_path = opts.right[i];
        time_print("\t%s", left_path);
        time_print("\t%s", right_path);
        //Read and parse the left and right file.
        FILTER_READ_PAIR_WORKERS workers(filter_read_matching, opts.threads << 4, opts.threads);
        //Start parsing the FASTQ pairs.
        FILTER_READ_PAIR_USER read_pair_filter_user {&workers};
        hmr_fastq_pair_read(left_path, right_path, filter_read_matching_submit, &read_pair_filter_user);
    }*/

    return 0;
}
