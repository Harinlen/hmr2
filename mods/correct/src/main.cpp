#include <cstdlib>
#include <cstdio>

#include "hmr_args.h"
#include "hmr_mapping.h"
#include "hmr_fasta.h"
#include "hmr_path.h"
#include "hmr_ui.h"
#include "hmr_parallel.h"

#include "args_correct.h"
#include "contig_correct.h"
#include "mapping_correct.h"
#include "mismatch_correct.h"

extern HMR_ARGS opts;

int main(int argc, char* argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.fasta) { help_exit(-1, "Missing FASTA file path."); }
    if (!path_can_read(opts.fasta)) { time_error(-1, "Cannot read FASTA file %s", opts.fasta); }
    if (opts.mappings.empty()) { help_exit(-1, "Missing Hi-C mapping file path."); }
    if (!opts.output) { help_exit(-1, "Missing corrected FASTA file path."); }
    //Try to write to output FASTA file.
    MISMATCH_CORRECTING corrected_file;
    mismatch_correct_open(opts.output, &corrected_file);
    //Load the FASTA name and length.
    time_print("Execution configuration:");
    time_print("\tMinimum Map Quality: %d", opts.mapq);
    time_print("\tMismatch resolutions: %d, %d, %d", opts.narrow, opts.wide, opts.depletion);
    time_print("\tThreads: %d", opts.threads);
    //Build the contig map.
    time_print("Building contig map from FASTA %s", opts.fasta);
    BAM_CORRECT_MAP correct_map;
    hmr_fasta_read(opts.fasta, contig_correct_build, &correct_map.contig_map);
    time_print("Contig map built, %zu contig(s) read.", correct_map.contig_map.size());
    //Prepare the mapping quality and pos lists.
    correct_map.narrow = opts.narrow;
    correct_map.wide = opts.wide;
    correct_map.mapq = opts.mapq;
    correct_map.wide_db = new HIC_DB[correct_map.contig_map.size()];
    correct_map.narrow_db = new HIC_DB[correct_map.contig_map.size()];
    //Loop and parse the mapping information.
    time_print("Constructing Hi-C reads relations...");
    for (char* mapping_path : opts.mappings)
    {
        time_print("Loading reads from %s", mapping_path);
        //Build the reads mapping.
        hmr_mapping_read(mapping_path, 
            MAPPING_PROC {mapping_correct_n_contig, mapping_correct_contig, mapping_correct_read_align}, 
            &correct_map, opts.threads);
        //Recover the mapping array.
        delete[] correct_map.bam_id_map.id;
    }
    time_print("Read(s) positions loaded and filtered.");
    //Calculate all the mismatches.
    time_print("Calculating mismatches...");
    int32_t contig_size = correct_map.contig_map.size();
    corrected_file.mismatches = new RANGE_LIST[contig_size];
    hmr_parallel_for(opts.threads, contig_size,
        mismatch_calc, correct_map.wide_db, correct_map.narrow_db, opts.percent, opts.sensitive, opts.depletion, opts.wide, opts.narrow, corrected_file.mismatches);
    //Now we are safe to remove databases.
    delete[] correct_map.wide_db;
    delete[] correct_map.narrow_db;
    time_print("Mismatches found.");
    //Based on the mismatches, render the corrected FASTA.
    time_print("Building the corrected FASTA file...");
    hmr_fasta_read(opts.fasta, mismatch_corrected, &corrected_file);
    //Flush the data.
    fclose(corrected_file.fp);
    time_print("Corrected FASTA has been written to %s", opts.output);
    return 0;
}