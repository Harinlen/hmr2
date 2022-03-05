#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>

#include "hmr_bam.h"
#include "hmr_ui.h"
#include "filter_fasta_type.h"

#include "filter_bam.h"

#define PAIR_INVALID 0xFFFFFFFFFFFFFFFF

typedef struct BAM_PAIR_DETAIL
{
    int32_t ref_id;
    int32_t ref_pos;
} BAM_PAIR_DETAIL;

typedef union BAM_PAIR_INFO
{
    BAM_PAIR_DETAIL detail;
    uint64_t data;
} BAM_PAIR_INFO;

typedef struct BAM_PAIR_INDEX
{
    uint64_t paired_data;
    size_t key_id;
} BAM_PAIR_INDEX;

typedef struct BAM_EDGE_DETAIL
{
    int32_t edge_start;
    int32_t edge_end;
} BAM_EDGE_DETAIL;

typedef union BAM_EDGE_INFO
{
    BAM_EDGE_DETAIL detail;
    uint64_t data;
} BAM_EDGE_INFO;

typedef struct BAM_INDEX
{
    int mapq;
    uint32_t n_ref, ref_ptr;
    int32_t* ref_id_map;
    const FASTA_ENZYME_POSES *poses;
    std::unordered_map<uint64_t, size_t> *edges;
    std::unordered_map<uint64_t, BAM_PAIR_INDEX> pair_map;
    FILE *read_file;
    char *read_buf;
    size_t read_buf_size;
    size_t read_buf_offset;
} BAM_INDEX;

void index_header(char *text, uint32_t l_text, void *user)
{
    (void)text;
    (void)l_text;
    (void)user;
}

void index_n_ref(uint32_t n_ref, void *user)
{
    BAM_INDEX *index = static_cast<BAM_INDEX *>(user);
    //Assign the n_ref.
    index->n_ref = n_ref;
    index->ref_ptr = 0;
    //Prepare the index mapping array.
    index->ref_id_map = static_cast<int32_t *>(malloc(n_ref * sizeof(int32_t)));
}

int32_t search_ref(char *name, uint32_t l_name, const FASTA_ENZYME_POSES *poses)
{
    for(size_t i=0; i<poses->seq_total; ++i)
    {
        if(l_name-1 != poses->seq_name_length[i])
        {
            continue;
        }
        //Compare the name length.
        if(strncmp(name, poses->seq_names[i], l_name-1) == 0)
        {
            return static_cast<int32_t>(i);
        }
    }
    return -1;
}

void index_ref(char *name, uint32_t l_name, uint32_t l_ref, void *user)
{
    (void)l_ref;
    BAM_INDEX *index = static_cast<BAM_INDEX *>(user);
    //Search the name in the poses.
    index->ref_id_map[index->ref_ptr++] = search_ref(name, l_name, index->poses);
}

inline bool index_in_range(size_t pos, ENZYME_RANGE *ranges, size_t range_size)
{
    for(size_t i=0; i<range_size; ++i)
    {
        if(pos >= ranges[i].start && pos <= ranges[i].end)
        {
            return true;
        }
    }
    return false;
}

typedef struct READ_INFO
{
    int32_t ref_id;
    int32_t pos;
    int32_t next_ref_id;
    int32_t next_pos;
} READ_INFO;

void index_block(size_t block_id, char *block, uint32_t block_size, void *user)
{
    (void)block_size;
    BAM_INDEX *index = static_cast<BAM_INDEX *>(user);
    //Recast the block data into header.
    BAM_BLOCK_HEADER *header = reinterpret_cast<BAM_BLOCK_HEADER *>(block);
    //Check map quality.
    if(header->mapq < index->mapq || header->mapq == 255)
    {
        return;
    }
    //Check whether the header mapping id is valid, and in enzyme range.
    int32_t my_ref_id = index->ref_id_map[header->refID],
            paired_ref_id = index->ref_id_map[header->next_refID];
    const FASTA_ENZYME_POSES *poses = index->poses;
    if(my_ref_id == -1 || paired_ref_id == -1 ||
            !index_in_range(header->pos,
                            poses->enz_ranges[my_ref_id],
                            poses->enz_range_size[my_ref_id]))
    {
        return;
    }
    //Check whether this pair contains in the data.
    BAM_PAIR_INFO my_info, paired_info;
    my_info.detail = BAM_PAIR_DETAIL {my_ref_id, header->pos};
    paired_info.detail = BAM_PAIR_DETAIL {paired_ref_id, header->next_pos};
    auto pair_result = index->pair_map.find(paired_info.data);
    if(pair_result == index->pair_map.end())
    {
        //No paired result, insert mine data to the index.
        index->pair_map.insert(std::make_pair(
                                   my_info.data,
                                   BAM_PAIR_INDEX {paired_info.data, block_id}));
    }
    else
    {
        //Check whether the paired result is the current read.
        if(pair_result->second.paired_data == my_info.data)
        {
            //Matched, send the data to write list, and invalid the record.
            BAM_EDGE_INFO edge_info;
            if(my_ref_id < paired_ref_id)
            {
                edge_info.detail.edge_start = my_ref_id;
                edge_info.detail.edge_end = paired_ref_id;
            }
            else
            {
                edge_info.detail.edge_start = paired_ref_id;
                edge_info.detail.edge_end = my_ref_id;
            }
            //Check whether the data exist in the edge map.
            auto edge_result = index->edges->find(edge_info.data);
            if(edge_result == index->edges->end())
            {
                //Failed to find this edge, push a new edge.
                index->edges->insert(std::make_pair(edge_info.data, 1));
            }
            else
            {
                ++edge_result->second;
            }
            //Dump the reads result info.
            if(edge_info.detail.edge_start != edge_info.detail.edge_end)
            {
                //Write the data to buffer.
                READ_INFO *read_data = reinterpret_cast<READ_INFO *>(index->read_buf + index->read_buf_offset);
                //Set the data.
                *read_data = READ_INFO{edge_info.detail.edge_start, header->pos,
                        edge_info.detail.edge_end, header->next_pos};
                index->read_buf_offset += sizeof(READ_INFO);
                //Check whether the offset reach the limitation.
                if(index->read_buf_offset >= index->read_buf_size)
                {
                    //Write the buffer.
                    fwrite(index->read_buf, index->read_buf_offset, 1, index->read_file);
                    index->read_buf_offset = 0;
                }
            }
            //Replace the mapping id.
            pair_result->second.paired_data = PAIR_INVALID;
            pair_result->second.key_id = PAIR_INVALID;
        }
    }
}

void filter_bam_statistic(const char *source, const FASTA_ENZYME_POSES &poses,
        int mapq, int threads,
        std::vector<READ_EDGE> &graph_edges, const char *reads_path)
{
    std::unordered_map<uint64_t, size_t> raw_edges;
    BAM_PARSER filter {index_header, index_n_ref, index_ref, index_block};
    BAM_INDEX index;
    index.mapq = mapq;
    index.poses = &poses;
    index.edges = &raw_edges;
    index.read_buf_size = sizeof(READ_INFO) << 20;
    index.read_buf = static_cast<char *>(malloc(index.read_buf_size));
    index.read_buf_offset = 0;
#ifdef _MSC_VER
    index.read_file = NULL;
    fopen_s(&index.read_file, reads_path, "wb");
#else
    FILE *index.read_file = fopen(reads_path, "wb");
#endif
    if(index.read_file == NULL)
    {
        time_error(-1, "Failed to open reads info output file %s", reads_path);
    }
    //Load the BAM file with using the filter parser.
    hmr_bam_load(source, &filter, threads, &index);
    //Flush the buffer.
    if(index.read_buf_offset > 0)
    {
        fwrite(index.read_buf, index.read_buf_offset, 1, index.read_file);
    }
    //Complete the data dump.
    fclose(index.read_file);
    //Free the buffer.
    free(index.read_buf);
    //Sort the edges.
    std::vector<READ_EDGE> edges;
    edges.reserve(raw_edges.size());
    for(auto &edge : raw_edges)
    {
        BAM_EDGE_INFO edge_info;
        edge_info.data = edge.first;
        edges.push_back(READ_EDGE {edge_info.detail.edge_start, edge_info.detail.edge_end,
                                   static_cast<double>(edge.second) /
                                                       static_cast<double>(poses.enz_range_size[edge_info.detail.edge_start] +
                                                       poses.enz_range_size[edge_info.detail.edge_end])});
    }
    //Sort the edges.
    std::sort(edges.begin(), edges.end(), [](const READ_EDGE &left, const READ_EDGE &right)
    {
        return (left.start == right.start) ? (left.end < right.end) : (left.start < right.start);
    });
    //Return the edges and half edges.
    graph_edges = edges;
}

void filter_bam_dump_edge(const char *filepath, const std::vector<READ_EDGE> &edges)
{
#ifdef _MSC_VER
    FILE *graph_file = NULL;
    fopen_s(&graph_file, filepath, "wb");
#else
    FILE *graph_file = fopen(filepath, "wb");
#endif
    if(graph_file == NULL)
    {
        time_error(-1, "Failed to open edge output file %s", filepath);
    }
    size_t edge_count = edges.size();
    fwrite(&edge_count, sizeof(size_t), 1, graph_file);
    for(auto &edge : edges)
    {
        //Parse the edge information.
        fwrite(&edge, sizeof(READ_EDGE), 1, graph_file);
    }
    fclose(graph_file);
}
