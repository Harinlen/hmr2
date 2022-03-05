#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "hmr_args.h"
#include "hmr_args_type.h"
#include "hmr_cgd_parse.h"
#include "hmr_ui.h"

#include "hmr_fasta.h"

extern HMR_ARGS opts;
static char output_path[1024];

bool suffix_match(const char *filepath, const char *suffix)
{
    size_t file_len = strlen(filepath);
    size_t suffix_len = strlen(suffix);
    if(file_len < suffix_len)
    {
        return false;
    }
    //Compare the last suffix size.
    size_t offset = file_len - suffix_len;
    for(size_t i=0; i<suffix_len; ++i)
    {
        if(filepath[offset + i] != suffix[i])
        {
            return false;
        }
    }
    return true;
}

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Print the current working parameter.
    printf("%s\n", opts.nodes);
    if(suffix_match(opts.nodes, ".node"))
    {
        //Read the node file.
        CONTIG_NODE *nodes;
        size_t node_size;
        hmr_node_read(opts.nodes, &nodes, &node_size);
        char *node_path = opts.output;
        if(node_path == nullptr)
        {
            sprintf(output_path, "%s.txt", opts.nodes);
            node_path = output_path;
        }
        //Write the node file.
#ifdef _MSC_VER
        FILE *node_file = NULL;
        fopen_s(&node_file, node_path, "w");
#else
        FILE *node_file = fopen(node_path, "w");
#endif
        for(size_t i=0; i<node_size; ++i)
        {
            fprintf(node_file, "%s\t%zu\n", nodes[i].name, nodes[i].length);
        }
        fclose(node_file);
        time_print_str("Node info dumps to %s.", node_path);
        return 0;
    }
    //No action takes.
    time_error_str(-1, "No action takes for file %s.", opts.nodes);
    return 0;
}
