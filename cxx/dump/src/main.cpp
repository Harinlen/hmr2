#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <unordered_map>
#include <string>

#include "hmr_args.h"
#include "hmr_args_type.h"
#include "hmr_cgd_parse.h"
#include "hmr_ui.h"

#include "hmr_fasta.h"

extern HMR_ARGS opts;
static char output_path[1024];

bool suffix_match(const char *filepath, const char *suffix)
{
    if(filepath==nullptr)
    {
        return false;
    }
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

bool provide_node()
{
    return suffix_match(opts.nodes, ".node");
}

bool provide_cluster()
{
    return suffix_match(opts.cluster, ".cluster");
}

bool provide_node_rename()
{
    return suffix_match(opts.node_rename, ".rename");
}

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    if(provide_node())
    {
        //Read the node file.
        CONTIG_NODE *nodes;
        size_t node_size;
        hmr_node_read(opts.nodes, &nodes, &node_size);
        //Check whether we provide the rename file.
        if(provide_node_rename())
        {
            //Read the rename rules.
#ifdef _MSC_VER
            FILE *rename_file = NULL;
            fopen_s(&rename_file, opts.node_rename, "r");
#else
            FILE *rename_file = fopen(opts.node_rename, "r");
#endif
            //Create the node matrix.
            std::unordered_map<std::string, int> name_map;
            for(size_t i=0; i<node_size; ++i)
            {
                name_map.insert(std::make_pair(std::string(nodes[i].name, nodes[i].l_name), static_cast<int>(i)));
            }
            int no_of_rules = 0;
            fscanf(rename_file, "%d", &no_of_rules);
            //Create the mapping rules.
            static char names[4096], u_names[4096];
            for(int i=0; i<no_of_rules; ++i)
            {
                fscanf(rename_file, "%s %s", names, u_names);
                size_t name_len = strlen(names), u_name_len = strlen(u_names);
                auto name_id = name_map.find(std::string(names, name_len));
                if(name_id == name_map.end())
                {
                    time_print("Rule ignored, failed to find %s", names);
                    continue;
                }
                //Replace the names.
                CONTIG_NODE &node = nodes[name_id->second];
                free(node.name);
                node.name = static_cast<char *>(malloc(u_name_len+1));
                node.l_name = u_name_len;
                strncpy(node.name, u_names, u_name_len);
                node.name[node.l_name] = '\0';
            }
            fclose(rename_file);
            //Write the renamed node file.
            char *node_path = opts.output;
            if(node_path == nullptr)
            {
                sprintf(output_path, "%s_rename.node", opts.nodes);
                node_path = output_path;
            }
#ifdef _MSC_VER
            FILE *dump_file = NULL;
            fopen_s(&dump_file, node_path, "wb");
#else
            FILE *dump_file = fopen(node_path, "wb");
#endif
            fwrite(&node_size, sizeof(size_t), 1, dump_file);
            for(size_t i=0; i<node_size; ++i)
            {
                const CONTIG_NODE &node = nodes[i];
                fwrite(&(node.l_name), sizeof(size_t), 1, dump_file);
                fwrite(node.name, node.l_name, 1, dump_file);
                fwrite(&(node.length), sizeof(size_t), 1, dump_file);
            }
            fclose(dump_file);
            time_print("Renamed node file dumps to %s", node_path);
            return 0;
        }
        //Check whether we provide cluster file.
        if(provide_cluster())
        {
            //Read the number of nodes.
            int no_of_nodes;
#ifdef _MSC_VER
            FILE *cluster_file = NULL;
            fopen_s(&cluster_file, opts.cluster, "rb");
#else
            FILE *cluster_file = fopen(opts.cluster, "rb");
#endif
            char *cluster_path = opts.output;
            if(cluster_path == nullptr)
            {
                sprintf(output_path, "%s.txt", opts.cluster);
                cluster_path = output_path;
            }
#ifdef _MSC_VER
            FILE *dump_file = NULL;
            fopen_s(&dump_file, cluster_path, "w");
#else
            FILE *dump_file = fopen(cluster_path, "w");
#endif
            fread(&no_of_nodes, sizeof(int), 1, cluster_file);
            int node_id;
            for(int i=0; i<no_of_nodes; ++i)
            {
                fread(&node_id, sizeof(int), 1, cluster_file);
                fprintf(dump_file, "%s\n", nodes[node_id].name);
            }
            fclose(dump_file);
            fclose(cluster_file);
            time_print("Cluster result dumps to %s", cluster_path);
            return 0;
        }
        char *node_path = opts.output;
        if(node_path == nullptr)
        {
            sprintf(output_path, "%s.txt", opts.nodes);
            node_path = output_path;
        }
        //Write the node file.
#ifdef _MSC_VER
        FILE *dump_file = NULL;
        fopen_s(&dump_file, node_path, "w");
#else
        FILE *dump_file = fopen(node_path, "w");
#endif
        for(size_t i=0; i<node_size; ++i)
        {
            fprintf(dump_file, "%s\t%zu\n", nodes[i].name, nodes[i].length);
        }
        fclose(dump_file);
        time_print("Node info dumps to %s", node_path);
        return 0;
    }
    //No action takes.
    time_error(-1, "No action takes for file %s.", opts.nodes);
    return 0;
}
