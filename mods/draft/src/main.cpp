#include <cstdlib>
#include <cstdio>

#include "hmr_args.h"
#include "hmr_ui.h"
#include "hmr_path.h"
#include "hmr_enzyme.h"
#include "hmr_fasta.h"
#include "hmr_mapping.h"
#include "hmr_contig_graph.h"
#include "hmr_bin_file.h"

#include "args_draft.h"
#include "fasta_draft.h"
#include "mapping_draft_type.h"
#include "mapping_draft.h"

extern HMR_ARGS opts;

int main(int argc, char* argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.fasta) { help_exit(-1, "Missing FASTA file path."); }
    if (!path_can_read(opts.fasta)) { time_error(-1, "Cannot read FASTA file %s", opts.fasta); }
    if (opts.mappings.empty()) { help_exit(-1, "Missing Hi-C mapping file path."); }
    if (!opts.enzyme) { help_exit(-1, "Missing restriction enzyme cutting site."); }
    if (!opts.output) { help_exit(-1, "Missing output file prefix."); }
    //Convert the enzyme into sequence.
    hmr_enzyme_formalize(opts.enzyme, &opts.enzyme_nuc, &opts.enzyme_nuc_length);
    time_print("Execution configuration:");
    time_print("\tMinimum map quality: %d", opts.mapq);
    time_print("\tRestriction enzyme: %s", opts.enzyme_nuc);
    time_print("\tMinimum restriction enzyme count: %d", opts.min_enzymes);
    time_print("\tHalf of enzyme range: %d", opts.range);
    time_print("\tThreads: %d", opts.threads);
    //Load the FASTA sequence and find the enzyme.
    HMR_CONTIGS contigs;
    HMR_CONTIG_INVALID_SET invalid_id_set;
    ENZYME_RANGES* contig_ranges;
    {
        //Prepare the enzyme for searching.
        ENZYME_SEARCH search;
        contig_draft_search_start(opts.enzyme_nuc, opts.enzyme_nuc_length, search);
        //Prepare the search user data.
        DRAFT_NODES_USER node_user{ opts.range, &contigs, &search, NULL, NULL, NULL };
        {
            //Prepare the thread pool for searching.
            RANGE_SEARCH_POOL search_pool(contig_range_search, opts.threads * 32, opts.threads);
            node_user.pool = &search_pool;
            time_print("Searching enzyme in %s", opts.fasta);
            hmr_fasta_read(opts.fasta, contig_draft_build, &node_user);
        }
        contig_draft_search_end(search);
        //Convert the search node information.
        contig_ranges = static_cast<ENZYME_RANGES*>(malloc(sizeof(ENZYME_RANGES) * contigs.size()));
        ENZYME_RANGE_CHAIN* chain_node = node_user.chain_head;
        int32_t range_index = 0;
        while (chain_node != NULL)
        {
            contig_ranges[range_index++] = chain_node->data;
            //Recover the chain node.
            ENZYME_RANGE_CHAIN* next = chain_node->next;
            free(chain_node);
            chain_node = next;
        }
        //Search complete.
        time_print("%zu contig(s) indexed.", contigs.size());
        //Dump the node data to target file.
        {
            std::string path_contig = hmr_graph_path_contig(opts.output);
            time_print("Save contig information to %s", path_contig.data());
            hmr_graph_save_contig(path_contig.data(), contigs);
            time_print("Done");
        }
        //Checking which contig is valid, if invalid, generate the invalid list.
        time_print("Checking invalid contig(s)...");
        HMR_CONTIG_INVALID_IDS invalid_ids;
        for (int32_t i = 0; i < range_index; ++i)
        {
            if (contig_ranges[i].counter < opts.min_enzymes)
            {
                invalid_ids.push_back(i);
            }
        }
        time_print("%zu invalid contig(s) detected.", invalid_ids.size());
        if (!invalid_ids.empty())
        {
            std::string path_invalid = hmr_graph_path_invalid(opts.output);
            time_print("Save invalid contig indices to %s", path_invalid.data());
            hmr_graph_save_invalid(path_invalid.data(), invalid_ids);
            time_print("Done");
            //Construct the invalid set.
            invalid_id_set = HMR_CONTIG_INVALID_SET(invalid_ids.begin(), invalid_ids.end());
        }
    }
    //Build the contig mapping.
    std::vector<HMR_EDGE_WEIGHT> edge_weights;
    {
        time_print("Building contig id map...");
        CONTIG_ID_MAP contig_ids;
        for (size_t i = 0; i < contigs.size(); ++i)
        {
            contig_ids.insert(std::make_pair(std::string(contigs[i].name, contigs[i].name_size),
                static_cast<int32_t>(i)));
            //Recover the memory.
            free(contigs[i].name);
        }
        time_print("Contig map has been built.");
        //Prepare the read-pair information output.
        std::string path_reads = hmr_graph_path_reads(opts.output);
        FILE* reads_file = NULL;
        if (!bin_open(path_reads.data(), &reads_file, "wb"))
        {
            time_error(-1, "Failed to create read information file %s", path_reads.data());
        }
        time_print("Writing reads summary information to %s", path_reads.data());
        //Loop and generate edge information.
        MAPPING_DRAFT_USER mapping_user{ READ_RECORD(), contig_ids, invalid_id_set, contig_ranges, NULL, 0, RAW_EDGE_MAP(), reads_file, static_cast<uint8_t>(opts.mapq), NULL, 0, 0 };
        time_print("Constructing Hi-C reads relations...");
        for (char* mapping_path : opts.mappings)
        {
            time_print("Loading reads from %s", mapping_path);
            //Build the reads mapping.
            hmr_mapping_read(mapping_path, MAPPING_PROC{ mapping_draft_n_contig, mapping_draft_contig, mapping_draft_read_align }, &mapping_user, opts.threads);
            //Recover the mapping array.
            delete[] mapping_user.contig_id_map;
        }
        time_print("Contig edges built from %zu file(s).", opts.mappings.size());
        //Check reads file buffer is complete.
        if (mapping_user.output_offset)
        {
            fwrite(mapping_user.output_buffer, mapping_user.output_offset, 1, reads_file);
        }
        fclose(reads_file);
        free(mapping_user.output_buffer);
        time_print("Reads summary information saved.");
        //Build the edge map.
        time_print("Calculating %zu edge weights...", mapping_user.edges.size());
        edge_weights = mapping_draft_get_edge_weights(mapping_user.edges, contig_ranges);
        time_print("Done");
    }
    //Dump the edge information into files.
    std::string path_edge = hmr_graph_path_edge(opts.output);
    time_print("Save contig edge information to %s", path_edge.data());
    hmr_graph_save_edge(path_edge.data(), edge_weights);
    time_print("Done");
    return 0;
}