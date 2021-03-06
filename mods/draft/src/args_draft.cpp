#include <cstdlib>

#include "args_draft.h"

#include "hmr_args_types.h"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-f", "--fasta"}, "FASTA", "Contig FASTA file (.fasta/.fasta.gz)", LAMBDA_PARSE_ARG {opts.fasta = arg[0]; }},
    { {"-m", "--mapping"}, "MAPPING 1, MAPPING 2...", "Hi-C reads mapping files (.bam/.hmr_mapping)", LAMBDA_PARSE_ARG { opts.mappings = arg; }},
    { {"-o", "--output"}, "OUTPUT", "Output graph prefix", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-e", "--enzyme"}, "ENZYME", "Enzyme to find in the sequence", LAMBDA_PARSE_ARG {opts.enzyme = arg[0];}},
    { {"-q", "--mapq"}, "MAPQ", "MAPQ of mapping lower bound (default: 1)", LAMBDA_PARSE_ARG {opts.mapq = atoi(arg[0]); }},
    { {"-r", "--range"}, "ENZYME_RANGE", "The enzyme position range size (default: 1000)", LAMBDA_PARSE_ARG {opts.range = atoi(arg[0]) >> 1; }},
    { {"-c", "--count"}, "ENZYME_COUNT", "The minimum enzyme count (default: 0)", LAMBDA_PARSE_ARG {opts.min_enzymes = atoi(arg[0]); }},
    { {"-t", "--threads"}, "THREAS", "Number of threads (default: 1)", LAMBDA_PARSE_ARG { opts.threads = atoi(arg[0]); }},
};
