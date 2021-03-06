#include <cstdlib>

#include "args_partition.h"

#include "hmr_args_types.h"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-n", "--nodes"}, "NDOES", "HMR contig node file (.hmr_contig)", LAMBDA_PARSE_ARG {opts.nodes = arg[0]; }},
    { {"-e", "--edge"}, "EDGE", "HMR edge file (.hmr_edge)", LAMBDA_PARSE_ARG { opts.edge = arg[0];}},
    { {"-g", "--group"}, "GROUP", "Number of homologous chromosomes groups", LAMBDA_PARSE_ARG {opts.groups = atoi(arg[0]); }},
    { {"-a", "--allele"}, "ALLELE_GROUP", "Number of allele chromosomes groups", LAMBDA_PARSE_ARG {opts.allele_groups = atoi(arg[0]); }},
    { {"-x", "--table"}, "ALLELE_TABLE", "Allele contig table (.hmr_allele/.ctg.table)", LAMBDA_PARSE_ARG {opts.allele_table = arg[0]; }},
    { {"-o", "--output"}, "OUTPUT", "Output partition file prefix", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-t", "--threads"}, "THREAS", "Number of threads (default: 1)", LAMBDA_PARSE_ARG { opts.threads = atoi(arg[0]); }},
};
