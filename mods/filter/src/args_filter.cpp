#include <cstdlib>

#include "args_filter.h"

#include "hmr_args_types.h"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-f", "--fasta"}, "FASTA", "Contig FASTA file (.fasta)", LAMBDA_PARSE_ARG {opts.fasta = arg[0]; }},
    { {"-e", "--enzyme"}, "ENZYME", "Enzyme to be found in the sequence", LAMBDA_PARSE_ARG {opts.enzyme = arg[0]; }},
    { {"-1", "--r1"}, "READ1", "One side of reads FASTQ file (.fastq/.gz)", LAMBDA_PARSE_ARG { opts.left = arg; }},
    { {"-2", "--r2"}, "READ2", "Corrected FASTA file (.fastq/.gz)", LAMBDA_PARSE_ARG { opts.right = arg; }},
    { {"-o", "--output"}, "OUTPUT", "Generated FASTQ prefix", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-t", "--threads"}, "THREAS", "Number of threads (default: 1)", LAMBDA_PARSE_ARG { opts.threads = atoi(arg[0]); }},
};
