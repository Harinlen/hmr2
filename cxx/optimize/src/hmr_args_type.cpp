#include <cstdlib>

#include "hmr_args_parser.h"
#include "hmr_args_type.h"

using namespace std;

HMR_ARGS opts;

HMR_ARG_PARSER args_parser
({
     { "-n", {"NODE", {"--node"}, "Graph node record file (.node)", [](char *arg){opts.node = arg;}}},
     { "-e", {"EDGE", {"--edge"}, "Graph edge record file (.edge)", [](char *arg){opts.edge = arg;}}},
     { "-c", {"CLUSTER", {"--cluster"}, "Group cluster record file (.clusters)", [](char *arg){opts.clusters = arg;}}},
     { "-o", {"OUTPUT", {"--output"}, "Contig partition prefix", [](char *arg){opts.output = arg;}}},
     { "-t", {"THREADS", {"--threads"}, "Number of threads (default: 1)", [](char *arg){opts.threads = atoi(arg);}}},
});

HMR_CONSTRAINTS args_constrains
{
    {[]{return opts.node != nullptr;}, "Failed to find the node record file."},
    {[]{return opts.edge != nullptr;}, "Failed to find the edge record file."},
    {[]{return opts.clusters != nullptr;}, "Failed to find the cluster record file."},
    {[]{return opts.output != nullptr;}, "Output graph file needs to be specified."},
};