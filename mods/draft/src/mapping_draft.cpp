#include "hmr_contig_graph.h"

#include "fasta_draft_type.h"

#include "mapping_draft.h"

void mapping_draft_n_contig(uint32_t n_ref, void* user)
{
    MAPPING_DRAFT_USER* mapping_user = reinterpret_cast<MAPPING_DRAFT_USER*>(user);
    //Reset the records.
    mapping_user->records = READ_RECORD();
    //Create the mapping user.
    mapping_user->contig_id_map = new int32_t[n_ref];
    mapping_user->contig_idx = 0;
}

void mapping_draft_contig(uint32_t name_length, char* name, uint32_t length, void* user)
{
    MAPPING_DRAFT_USER* mapping_user = reinterpret_cast<MAPPING_DRAFT_USER*>(user);
    //Find the name in the contig mapping.
    auto name_finder = mapping_user->contig_ids.find(std::string(name, name_length));
    mapping_user->contig_id_map[mapping_user->contig_idx] = name_finder == mapping_user->contig_ids.end() ? -1 : name_finder->second;
    ++mapping_user->contig_idx;
}

bool position_in_range(int32_t pos, const ENZYME_RANGES& ranges)
{
    int32_t start = 0, end = ranges.length;
    do
    {
        int32_t range_ptr = (start + end) >> 1;
        const ENZYME_RANGE& r = ranges.ranges[range_ptr];
        //If it is in the range.
        if (r.start <= pos && pos <= r.end)
        {
            return true;
        }
        //Check whether we only have two ranges.
        if (end == start + 1)
        {
            const ENZYME_RANGE& rn = ranges.ranges[range_ptr];
            return rn.start <= pos && pos <= rn.end;
        }
        //Change the start or end position.
        if (pos < r.start)
        {
            end = range_ptr;
        }
        else
        {
            start = range_ptr;
        }
    } while (start != end);
    //No position matched.
    return false;
}

void mapping_draft_read_align(size_t id, const MAPPING_INFO& mapping_info, void* user)
{
    MAPPING_DRAFT_USER* mapping_user = reinterpret_cast<MAPPING_DRAFT_USER*>(user);
    //Check whether the mapping meets the requirements.
    if (mapping_info.mapq < mapping_user->mapq || mapping_info.mapq == 255) //Map quality is invalid
    {
        return;
    }
    //Find the ranges.
    int32_t ref_index = mapping_user->contig_id_map[mapping_info.refID],
        next_ref_index = mapping_user->contig_id_map[mapping_info.next_refID];
    if (ref_index == -1 || next_ref_index == -1 || //Reference index invalid.
        !position_in_range(mapping_info.pos, mapping_user->contig_ranges[ref_index]) || //Not in range.
        !position_in_range(mapping_info.next_pos, mapping_user->contig_ranges[next_ref_index]))
    {
        return;
    }
    //Check whether the paired read is in the record.
    MAPPING_READ ref_read{ {mapping_info.refID, mapping_info.pos} },
        next_ref_read{ {mapping_info.next_refID, mapping_info.next_pos} };
    //Search the records inside the map.
    auto next_finder = mapping_user->records.find(next_ref_read.data);
    if (next_finder == mapping_user->records.end())
    {
        //Insert the current paired information into the records.
        mapping_user->records.insert(std::make_pair(ref_read.data, next_ref_read));
    }
    else
    {
        //Check whether the records are paired.
        if (next_finder->second.data == ref_read.data)
        {
            //Replace the next finder as an invalid data.
            next_finder->second.data = 0xFFFFFFFFFFFFFFFF;
            //Paired information are found.
            HMR_EDGE edge = hmr_graph_edge(mapping_info.refID, mapping_info.next_refID);
            auto edge_finder = mapping_user->edges.find(edge.data);
            if (edge_finder == mapping_user->edges.end())
            {
                //Insert a new record.
                mapping_user->edges.insert(std::make_pair(edge.data, 1));
            }
            else
            {
                //Increase the record.
                ++edge_finder->second;
            }
            //Construct and write the mapping info to the reads file.
            fwrite(&mapping_info, sizeof(MAPPING_INFO), 1, mapping_user->reads_file);
        }
    }
}

std::vector<HMR_EDGE_WEIGHT> mapping_draft_get_edge_weights(const RAW_EDGE_MAP& edge_map, const HMR_CONTIGS& contigs)
{
    std::vector<HMR_EDGE_WEIGHT> weights;
    weights.reserve(edge_map.size());
    //Loop and get edge counting.
    for (const auto& edge_counter : edge_map)
    {
        HMR_EDGE edge;
        edge.data = edge_counter.first;
        weights.push_back(HMR_EDGE_WEIGHT{ edge, static_cast<double>(edge_counter.second) / static_cast<double>(contigs[edge.pos.start].length + contigs[edge.pos.end].length) });
    }
    return weights;
}
