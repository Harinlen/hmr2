#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>

#include "hmr_ui.h"
#include "hmr_global.h"
#include "hmr_polar_type.h"
#include "hmr_cgd_type.h"
#include "hmr_reads_parse.h"

typedef struct READ_INFO
{
    int32_t ref_id;
    int32_t pos;
    int32_t next_ref_id;
    int32_t next_pos;
} READ_INFO;

POLAR_INFO hmr_polar_parse(const char *filepath, const std::vector<int> &sequence, CONTIG_NODE *nodes)
{
#ifdef _MSC_VER
    FILE *reads_file = NULL;
    fopen_s(&reads_file, filepath, "rb");
#else
    FILE *reads_file = fopen(filepath, "rb");
#endif
    if(reads_file == NULL)
    {
        time_error(-1, "Failed to read contig polarity file %s", filepath);
    }
    POLAR_INFO node_polars;
    node_polars.resize(sequence.size());
    std::unordered_map<int, NODE_POLAR *> polar_map;
    for(size_t i=0; i<sequence.size(); ++i)
    {
        node_polars[i].node_id = sequence[i];
        node_polars[i].length = static_cast<int>(nodes[sequence[i]].length);
        node_polars[i].seq_pos = static_cast<int>(i);
        polar_map.insert(std::make_pair(sequence[i], &node_polars[i]));
    }
    //Read the contig size of the file.
    READ_INFO read_info;
    while(fread(&read_info, sizeof(READ_INFO), 1, reads_file))
    {
        //Read the node id.
        auto ref_finder = polar_map.find(read_info.ref_id),
                next_ref_finder = polar_map.find(read_info.next_ref_id);
        if(ref_finder != polar_map.end() && next_ref_finder != polar_map.end())
        {
            //Insert the polar info.
            auto ref = ref_finder->second, next_ref = next_ref_finder->second;
            ref->reads.push_back(READ_PAIR {read_info.pos, next_ref, read_info.next_pos});
            next_ref->reads.push_back(READ_PAIR {read_info.next_pos, ref, read_info.pos});
        }
    }
    fclose(reads_file);
    return node_polars;
}
