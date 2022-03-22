#include <cstdlib>

#include "args_correct.h"

#include "hmr_args_types.h"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-f", "--fasta"}, "FASTA", "Contig FASTA file (.fasta/.fasta.gz)", LAMBDA_PARSE_ARG {opts.fasta = arg[0]; }},
    { {"-m", "--mapping"}, "MAPPING 1, MAPPING 2...", "Hi-C reads mapping files (.bam/.hmr_mapping)", LAMBDA_PARSE_ARG { opts.mappings = arg; }},
    { {"-o", "--output"}, "OUTPUT", "Corrected contig FASTA file (.fasta)", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-p", "--percent"}, "PERCENT", "Percent of the map to saturate (default: 0.95)", LAMBDA_PARSE_ARG {opts.percent = atof(arg[0]);}},
    { {"-s", "--sensitive"}, "SENSITIVE", "Sensitivity to depletion score (default: 0.5)", LAMBDA_PARSE_ARG {opts.sensitive = atof(arg[0]); }},
    { {"-q", "--mapq"}, "MAPQ", "MAPQ of mapping lower bound (default: 1)", LAMBDA_PARSE_ARG {opts.mapq = atoi(arg[0]); }},
    { {"-w", "--wide"}, "WIDE", "Resolution for first pass search of mismatches (default: 25000)", LAMBDA_PARSE_ARG {opts.wide = atoi(arg[0]); }},
    { {"-n", "--narrow"}, "NARROW", "Resolution for the precise mismatch localizaton, NARROW < WIDE (default: 1000)", LAMBDA_PARSE_ARG {opts.narrow = atoi(arg[0]); }},
    { {"-d", "--depletion"}, "DEPLETION", "The size of the region to aggregate the depletion score in the wide path, DEPLETION >= 2 * WIDE (default: 100000)", LAMBDA_PARSE_ARG {opts.depletion = atoi(arg[0]); }},
    { {"-t", "--threads"}, "THREAS", "Number of threads (default: 1)", LAMBDA_PARSE_ARG { opts.threads = atoi(arg[0]); }},
};
