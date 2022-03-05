#include <cstdlib>

#include "hmr_args_parser.h"
#include "hmr_args_type.h"

using namespace std;

HMR_ARGS opts;

HMR_ARG_PARSER args_parser
({
     { "-n", {"NODE", {"--node"}, "Graph node record file (.node)", [](char *arg){opts.node = arg;}}},
     { "-e", {"EDGE", {"--edge"}, "Graph edge record file (.edge)", [](char *arg){opts.edge = arg;}}},
     { "-r", {"READS", {"--reads"}, "Contig reads index file (.reads)", [](char *arg){opts.reads = arg;}}},
     { "-c", {"CLUSTER", {"--cluster"}, "Group cluster record file (.clusters)", [](char *arg){opts.cluster = arg;}}},
     { "-o", {"OUTPUT", {"--output"}, "Group tour file (.tour)", [](char *arg){opts.output = arg;}}},
     { "-t", {"THREADS", {"--threads"}, "Number of threads (default: 1)", [](char *arg){opts.threads = atoi(arg);}}},
});

HMR_CONSTRAINTS args_constrains
{
    {[]{return opts.node != nullptr;}, "Failed to find the node record file."},
    {[]{return opts.edge != nullptr;}, "Failed to find the edge record file."},
    {[]{return opts.reads != nullptr;}, "Failed to find the reads index record file."},
    {[]{return opts.cluster != nullptr;}, "Failed to find the cluster record file."},
    {[]{return opts.output != nullptr;}, "Output tour file needs to be specified."},
};
