#include <cassert>

#include "hmr_bin_file.h"
#include "hmr_ui.h"
#include "hmr_path.h"

#include "hmr_contig_graph.h"

HMR_EDGE hmr_graph_edge(int32_t a, int32_t b)
{
    HMR_EDGE edge{};
    if (a < b)
    {
        edge.pos.start = a;
        edge.pos.end = b;
    }
    else
    {
        edge.pos.start = b;
        edge.pos.end = a;
    }
    return edge;
}

std::string hmr_graph_path_contig(const char* prefix)
{
    return std::string(prefix) + ".hmr_contig";
}

bool hmr_graph_load_contig(const char* filepath, HMR_CONTIGS& contigs)
{
    //Open the contig input file to read the data.
    FILE* contig_file;
    if (!bin_open(filepath, &contig_file, "rb"))
    {
        time_error(-1, "Failed to load contig information from %s", filepath);
        return false;
    }
    //Read the length of the contigs.
    size_t contig_sizes = 0;
    fread(&contig_sizes, sizeof(size_t), 1, contig_file);
    contigs = std::vector<HMR_CONTIG>();
    contigs.reserve(contig_sizes);
    for (size_t i = 0; i < contig_sizes; ++i)
    {
        //Read the contig from the file.
        HMR_CONTIG contig;
        fread(&contig.name_size, sizeof(int32_t), 1, contig_file);
        contig.name = static_cast<char*>(malloc(contig.name_size));
        assert(contig.name);
        fread(contig.name, sizeof(char), contig.name_size, contig_file);
        fread(&contig.length, sizeof(int32_t), 1, contig_file);
        contigs.push_back(contig);
    }
    fclose(contig_file);
    return true;
}

bool hmr_graph_save_contig(const char* filepath, const HMR_CONTIGS& contigs)
{
    //Open the contig output file to write the data.
    FILE* contig_file;
    if (!bin_open(filepath, &contig_file, "wb"))
    {
        time_error(-1, "Failed to save contig information from %s", filepath);
        return false;
    }
    //Write the length of the contigs.
    size_t contig_sizes = contigs.size();
    fwrite(&contig_sizes, sizeof(size_t), 1, contig_file);
    for (const auto& contig : contigs)
    {
        //Write the contig information.
        fwrite(&contig.name_size, sizeof(int32_t), 1, contig_file);
        fwrite(contig.name, sizeof(char), contig.name_size, contig_file);
        fwrite(&contig.length, sizeof(int32_t), 1, contig_file);
    }
    fclose(contig_file);
    return true;
}

std::string hmr_graph_path_reads(const char* prefix)
{
    return std::string(prefix) + ".hmr_reads";
}

std::string hmr_graph_path_edge(const char* prefix)
{
    return std::string(prefix) + ".hmr_edge";
}

bool hmr_graph_save_edge(const char* filepath, const HMR_EDGE_WEIGHTS& edges)
{
    //Open the contig output file to write the data.
    FILE* edge_file;
    if (!bin_open(filepath, &edge_file, "wb"))
    {
        time_error(-1, "Failed to save edge information from %s", filepath);
        return false;
    }
    //Write the edge information.
    size_t contig_sizes = edges.size();
    fwrite(&contig_sizes, sizeof(size_t), 1, edge_file);
    for (const auto& edge : edges)
    {
        fwrite(&edge, sizeof(HMR_EDGE_WEIGHT), 1, edge_file);
    }
    fclose(edge_file);
    return true;
}
