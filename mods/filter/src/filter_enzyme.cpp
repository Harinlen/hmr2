#define _CRT_SECURE_NO_WARNINGS

#include <cassert>
#include <list>

#include "hmr_fastq.h"

#include "filter_enzyme.h"

void filter_enzyme_search_prepare(const char* enzyme, size_t enzyme_size, FILTER_ENZYME_SEARCH* search)
{
    //Save the enzyme and its size.
    search->enzyme = enzyme;
    search->enzyme_size = enzyme_size;
}

void filter_enzyme_search_submit(char* seq_name, size_t seq_name_length, char* seq, size_t seq_size, void* user)
{
    //Push the data to worker's pool.
    FILTER_ENZYME_USER *filter_user = reinterpret_cast<FILTER_ENZYME_USER *>(user);
    filter_user->workers->push_task(
        FILTER_ENZYME_ARGS 
        {
            seq_name, seq_name_length, seq, seq_size, 
            filter_user->search, filter_user->ranges
        });
}

void filter_enzyme_search(const FILTER_ENZYME_ARGS &args)
{
    FILTER_ENZYME_SEARCH* search = args.search;
    //Search the enzyme inside the sequence.
    std::vector<FILTER_ENZYME_RANGE> seq_ranges;
    {
        size_t offset = 0, max_offset = args.seq_size - search->enzyme_size;
        char* seq_pos = args.seq + offset, * pos = strstr(seq_pos, search->enzyme);
        std::list<size_t> enzyme_poses;
        while (offset < max_offset && pos != NULL)
        {
            //Save the enzyme position.
            offset = pos - args.seq;
            enzyme_poses.push_back(offset);
            //Set the next start position.
            seq_pos = args.seq + offset + search->enzyme_size;
            pos = strstr(seq_pos, search->enzyme);
        }
        //Extract the +-500bp of the position.
        seq_ranges.reserve(enzyme_poses.size());
        for (size_t enzyme_offset : enzyme_poses)
        {
            //Create the enzyme range.
            size_t start_pos = (enzyme_offset < 500) ? 0 : (enzyme_offset - 500),
                end_pos = (enzyme_offset + 500 < args.seq_size) ? (enzyme_offset + 500) : args.seq_size;
            size_t range_length = end_pos - start_pos;
            seq_ranges.push_back(FILTER_ENZYME_RANGE{ start_pos, end_pos, range_length });
        }
    }
    //Insert the data to ranges.
    {
        //Construct the CONTIG.
        FILTER_ENZYME_RANGES* ranges = args.ranges;
        //Lock the data.
        std::unique_lock<std::mutex> range_lock(*(ranges->access_mutex));
        //Push the data.
        ranges->contigs->push_back(HMR_CONTIG{args.seq_name, args.seq_name_length, args.seq_size});
        ranges->enzyme_ranges->push_back(FILTER_ENZYME_RANGE_DATA { seq_ranges, args.seq, args.seq_size} );
    }
}

void filter_read_matching_submit(const FASTQ_READ_PAIR& read_pair, void* user)
{
    FILTER_READ_PAIR_USER* reads_user = reinterpret_cast<FILTER_READ_PAIR_USER*>(user);
    //Push the pair to working thread.
    reads_user->workers->push_task(FILTER_READ_PAIR_ARGS {read_pair});
}

void filter_read_matching(const FILTER_READ_PAIR_ARGS& args)
{
    //Try to find whether this pair could be 90% similar to any enzyme range.
    ;
    printf("%zu %zu\n", args.read_pair.left.seq_length, args.read_pair.right.seq_length);
    hmr_fastq_read_free(args.read_pair.left);
    hmr_fastq_read_free(args.read_pair.right);
}
