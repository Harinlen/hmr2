#include <cstdio>

#include "hmr_ui.h"
#include "mapping_correct_type.h"

#include "mapping_correct.h"

void mapping_correct_n_contig(uint32_t n_ref, void* user)
{
    BAM_CORRECT_MAP* bam_map = static_cast<BAM_CORRECT_MAP*>(user);
    //Initialize the contig id.
    bam_map->bam_id_map.size = n_ref;
    bam_map->bam_id_map.id = new int32_t[n_ref];
    for (int i = 0; i < n_ref; ++i)
    {
        bam_map->bam_id_map.id[i] = -1;
    }
    bam_map->bam_contig_id = 0;
}

void mapping_correct_contig(uint32_t name_length, char* name, uint32_t length, void* user)
{
    BAM_CORRECT_MAP* bam_map = static_cast<BAM_CORRECT_MAP*>(user);
    //Find the name in the contig map.
    auto contig_finder = bam_map->contig_map.find(std::string(name, name_length));
    if (contig_finder != bam_map->contig_map.end())
    {
        //Assign the contig id.
        bam_map->bam_id_map.id[bam_map->bam_contig_id] = (*contig_finder).second.id;
        ++bam_map->bam_contig_id;
    }
    else
    {
        time_print("Failed to find contig '%s'", name);
    }
}

inline void count_pair(int32_t a, int32_t b, HIC_DB& db)
{
    POS_PAIR reads_pair;
    reads_pair.pos = { a, b };
    auto pos_finder = db.find(reads_pair.data);
    if (pos_finder == db.end())
    {
        db.insert(std::make_pair(reads_pair.data, 1));
    }
    else
    {
        ++(*pos_finder).second;
    }
}

void mapping_correct_read_align(size_t id, const MAPPING_INFO& mapping_info, void* user)
{
    BAM_CORRECT_MAP* bam_map = static_cast<BAM_CORRECT_MAP*>(user);
    //Check the mapq reaches the limitation.
    if (mapping_info.mapq < bam_map->mapq || // Quality filter.
        mapping_info.pos == -1 || mapping_info.next_pos == -1 || // Mapped.
        mapping_info.refID != mapping_info.next_refID) //In the same contig.
    {
        return;
    }
    int32_t target_id = bam_map->bam_id_map.id[mapping_info.refID];
    if (target_id == -1)
    {
        return;
    }
    //Count at wide database.
    count_pair(mapping_info.pos / bam_map->wide * bam_map->wide,
        mapping_info.next_pos / bam_map->wide * bam_map->wide, bam_map->wide_db[target_id]);
    count_pair(mapping_info.pos / bam_map->narrow * bam_map->narrow,
        mapping_info.next_pos / bam_map->narrow * bam_map->narrow, bam_map->narrow_db[target_id]);
}
