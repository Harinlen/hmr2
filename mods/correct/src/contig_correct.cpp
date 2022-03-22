#include "hmr_ui.h"

#include "contig_correct_type.h"

#include "contig_correct.h"

void contig_correct_build(int32_t index, char* seq_name, size_t seq_name_size, char* seq, size_t seq_size, void* user)
{
    //Recast the user into map.
    CONTIG_MAP* map = reinterpret_cast<CONTIG_MAP*>(user);
    //Record the sequence name and length.
    std::string contig_name(seq_name, seq_name_size);
    auto contig_finder = map->find(contig_name);
    if (contig_finder != map->end())
    {
        time_error(-1, "Corrupted FASTA: contig '%s' already existed.", contig_name.data());
    }
    //Insert the contig into map.
    map->insert(std::make_pair(contig_name, CONTIG_INFO {index, seq_size} ));
    free(seq_name);
    free(seq);
}
