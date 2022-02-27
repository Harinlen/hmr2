#include <cstdlib>

#include "hmr_args_parser.h"
#include "hmr_args_type.h"

using namespace std;

HMR_ARGS opts;

HMR_ARG_PARSER args_parser
({
     { "-f", {"FASTA", {"--fasta"}, "Contig fasta file (.fasta)", [](char *arg){opts.fasta = arg;}}},
     { "-m", {"MAPPING", {"--mapping"}, "Input mapping file (.bam)", [](char *arg){opts.mapping = arg;}}},
     { "-e", {"ENZYME", {"--enzyme"}, "Enzyme to find in the sequence", [](char *arg){opts.enzyme = arg;}}},
     { "-o", {"OUTPUT", {"--output"}, "Contig graph file", [](char *arg){opts.output = arg;}}},
     { "-q", {"MAPQ", {"--mapq"}, "MAPQ of mapping lower bound (default: 40)", [](char *arg){opts.mapq = atoi(arg);}}},
     { "-t", {"THREADS", {"--threads"}, "Number of threads (default: 1)", [](char *arg){opts.threads = atoi(arg);}}},
});

HMR_CONSTRAINTS args_constrains
{
    {[]{return opts.fasta != nullptr;}, "Failed to find fasta file."},
    {[]{return opts.mapping != nullptr;}, "Failed to find the mapping file."},
    {[]{return opts.output != nullptr;}, "Output graph file needs to be specified."},
    {[]{return opts.enzyme != nullptr;}, "Enzyme is not specified."},
};
