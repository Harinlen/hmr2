#include <cstdlib>

#include "hmr_args_parser.h"
#include "hmr_args_type.h"

using namespace std;

HMR_ARGS opts;

HMR_ARG_PARSER args_parser
({
     { "-n", {"NODE", {"--node"}, "HMR node binary file", [](char *arg){opts.nodes = arg;}}},
     { "-nr", {"NODE_RENAME", {"--node-rename"}, "Node rename rule file", [](char *arg){opts.node_rename = arg;}}},
     { "-c", {"CLUSTER", {"--cluster"}, "HMR node cluster binary file", [](char *arg){opts.cluster = arg;}}},
     { "-f", {"FASTA", {"--fasta"}, "Contig FASTA file", [](char *arg){opts.fasta = arg;}}},
     { "-t", {"TOUR", {"--tour"}, "Sequence order and direction description file", [](char *arg){opts.cluster = arg;}}},
     { "-o", {"OUTPUT", {"--output"}, "Output file name", [](char *arg){opts.output = arg;}}},
});

HMR_CONSTRAINTS args_constrains
{
    {[]{return opts.nodes != nullptr;}, "Please specify the HMR node binary file."},
};
